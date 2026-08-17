//! Provider-owned backward-facing-step geometry, solve, and wall-shear metrics.
//!
//! The validation crate consumes this module as a thin benchmark adapter. The
//! geometry mask and the signed wall-shear extraction stay beside the SIMPLE
//! solver so a second consumer cannot silently implement a different step
//! boundary or reattachment definition.

use crate::scalar;
use crate::solvers::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, SolveResult};
use cfd_core::error::{Error, Result};
use cfd_core::geometry::StaggeredGrid2D;
use cfd_core::CfdScalar;
use eunomia::{FloatElement, RealField};

/// Physical dimensions of a backward-facing-step channel.
#[derive(Debug, Clone, Copy)]
pub struct BackwardFacingStepGeometry<T> {
    /// Height of the upstream solid step.
    pub step_height: T,
    /// Height of the downstream channel.
    pub channel_height: T,
    /// Fluid length before the step face.
    pub upstream_length: T,
    /// Fluid length after the step face.
    pub downstream_length: T,
    /// Uniform inlet velocity in the downstream fluid height.
    pub inlet_velocity: T,
}

impl<T: CfdScalar + Copy + FloatElement> BackwardFacingStepGeometry<T> {
    /// Construct and validate a backward-facing-step geometry.
    pub fn new(
        step_height: T,
        channel_height: T,
        upstream_length: T,
        downstream_length: T,
        inlet_velocity: T,
    ) -> Result<Self> {
        let zero = scalar::zero::<T>();
        let finite_positive = |value: T| value.is_finite() && value > zero;
        if !finite_positive(step_height)
            || !finite_positive(channel_height)
            || !finite_positive(upstream_length)
            || !finite_positive(downstream_length)
            || !finite_positive(inlet_velocity)
            || channel_height <= step_height
        {
            return Err(Error::InvalidConfiguration(
                "backward-facing-step dimensions and inlet velocity must be finite and positive; channel height must exceed step height".to_string(),
            ));
        }
        Ok(Self {
            step_height,
            channel_height,
            upstream_length,
            downstream_length,
            inlet_velocity,
        })
    }

    /// Total streamwise domain length.
    #[must_use]
    pub fn channel_length(self) -> T {
        self.upstream_length + self.downstream_length
    }
}

/// Inputs for one provider-owned backward-facing-step solve.
#[derive(Debug, Clone)]
pub struct BackwardFacingStepConfig<T: RealField + Copy> {
    /// Number of cells across the downstream channel height.
    pub resolution: usize,
    /// Reynolds number based on step height and inlet velocity.
    pub reynolds_number: T,
    /// SIMPLE pressure-velocity coupling configuration.
    pub simple: SIMPLEConfig<T>,
}

/// One signed lower-wall shear sample downstream of the step.
#[derive(Debug, Clone, Copy)]
pub struct WallShearSample<T> {
    /// Distance from the step face.
    pub x: T,
    /// Signed streamwise wall shear rate proxy `du/dy` at the lower wall.
    pub shear_rate: T,
}

/// Result of a backward-facing-step solve and its reattachment metric.
#[derive(Debug, Clone)]
pub struct BackwardFacingStepResult<T> {
    /// SIMPLE convergence result.
    pub solve: SolveResult<T>,
    /// Distance from the step face to the first negative-to-nonnegative shear crossing.
    pub reattachment_length: T,
    /// Signed lower-wall samples used to derive `reattachment_length`.
    pub wall_shear: Vec<WallShearSample<T>>,
}

/// Provider-owned backward-facing-step SIMPLE solver.
pub struct BackwardFacingStepSolver<T> {
    geometry: BackwardFacingStepGeometry<T>,
}

impl<T: CfdScalar + Copy + FloatElement> BackwardFacingStepSolver<T> {
    /// Create a solver for a validated backward-facing-step geometry.
    #[must_use]
    pub fn new(geometry: BackwardFacingStepGeometry<T>) -> Self {
        Self { geometry }
    }

    /// Solve the masked step domain and derive reattachment from signed wall shear.
    ///
    /// The mask marks the upstream region below the step as solid. The solver
    /// applies a parabolic velocity inlet on the fluid cells, zero normal
    /// velocity and no-slip walls on solid/physical boundaries, and a zero
    /// streamwise-gradient velocity with fixed-pressure outlet through the
    /// canonical SIMPLE boundary path. Reattachment is not inferred from a
    /// correlation: it is interpolated from the first downstream lower-wall
    /// `du/dy` sign change in the solved field.
    pub fn solve(
        &self,
        config: &BackwardFacingStepConfig<T>,
    ) -> Result<BackwardFacingStepResult<T>> {
        self.validate_config(config)?;
        let nx = config.resolution.checked_mul(3).ok_or_else(|| {
            Error::InvalidConfiguration("backward-facing-step grid dimensions overflow".to_string())
        })?;
        let ny = config.resolution;
        let grid = StaggeredGrid2D::new(
            nx,
            ny,
            self.geometry.channel_length(),
            self.geometry.channel_height,
        );
        let viscosity = config.reynolds_number.recip()
            * self.geometry.inlet_velocity
            * self.geometry.step_height;
        let mut solver = NavierStokesSolver2D::new(
            grid,
            BloodModel::Newtonian(viscosity),
            scalar::one::<T>(),
            config.simple.clone(),
        );
        self.apply_geometry_mask(&mut solver);
        let solve = solver.solve(self.geometry.inlet_velocity)?;
        let wall_shear = collect_downstream_wall_shear(&solver, self.geometry.upstream_length)?;
        let reattachment_length = interpolate_reattachment(&wall_shear)?;
        Ok(BackwardFacingStepResult {
            solve,
            reattachment_length,
            wall_shear,
        })
    }

    fn validate_config(&self, config: &BackwardFacingStepConfig<T>) -> Result<()> {
        let zero = scalar::zero::<T>();
        if config.resolution < 4
            || config.simple.max_iterations == 0
            || !config.reynolds_number.is_finite()
            || config.reynolds_number <= zero
            || !config.simple.tolerance.is_finite()
            || config.simple.tolerance <= zero
        {
            return Err(Error::InvalidConfiguration(
                "backward-facing-step resolution, SIMPLE iterations, tolerance, and Reynolds number must be positive and finite".to_string(),
            ));
        }
        Ok(())
    }

    fn apply_geometry_mask(&self, solver: &mut NavierStokesSolver2D<T>) {
        for i in 0..solver.grid.nx {
            let x = solver.grid.x_center(i);
            for j in 0..solver.grid.ny {
                let y = solver.grid.y_center(j);
                solver.field.mask[(i, j)] =
                    x >= self.geometry.upstream_length || y >= self.geometry.step_height;
            }
        }
    }
}

fn collect_downstream_wall_shear<T: CfdScalar + Copy + FloatElement>(
    solver: &NavierStokesSolver2D<T>,
    upstream_length: T,
) -> Result<Vec<WallShearSample<T>>> {
    let zero = scalar::zero::<T>();
    let half = <T as FloatElement>::from_f64(0.5);
    let mut samples = Vec::new();
    for i in 0..solver.grid.nx {
        if !solver.field.mask[(i, 0)] {
            continue;
        }
        let x = solver.grid.x_center(i) - upstream_length;
        if x < zero {
            continue;
        }
        let u_center = (solver.field.u[(i, 0)] + solver.field.u[(i + 1, 0)]) * half;
        let wall_distance = solver.grid.dy_at(0) * half;
        let shear_rate = u_center / wall_distance;
        if !x.is_finite() || !shear_rate.is_finite() {
            return Err(Error::InvalidInput(
                "backward-facing-step lower-wall shear contains a non-finite value".to_string(),
            ));
        }
        samples.push(WallShearSample { x, shear_rate });
    }
    if samples.len() < 2 {
        return Err(Error::InvalidInput(
            "backward-facing-step requires at least two downstream wall-shear samples".to_string(),
        ));
    }
    Ok(samples)
}

fn interpolate_reattachment<T: CfdScalar + Copy + FloatElement>(
    samples: &[WallShearSample<T>],
) -> Result<T> {
    let zero = scalar::zero::<T>();
    for pair in samples.windows(2) {
        let [previous, current] = pair else {
            continue;
        };
        if previous.shear_rate < zero && current.shear_rate >= zero {
            let fraction =
                (zero - previous.shear_rate) / (current.shear_rate - previous.shear_rate);
            return Ok(previous.x + fraction * (current.x - previous.x));
        }
    }
    Err(Error::InvalidInput(
        "backward-facing-step field has no downstream wall-shear reattachment crossing".to_string(),
    ))
}

#[cfg(test)]
mod tests {
    use super::{
        interpolate_reattachment, BackwardFacingStepGeometry, BackwardFacingStepSolver,
        WallShearSample,
    };
    use crate::solvers::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig};
    use cfd_core::geometry::StaggeredGrid2D;

    #[test]
    fn wall_shear_crossing_is_linearly_interpolated() {
        let samples = [
            WallShearSample {
                x: 1.0,
                shear_rate: -1.0,
            },
            WallShearSample {
                x: 2.0,
                shear_rate: 1.0,
            },
        ];
        let result = interpolate_reattachment(&samples).expect("crossing exists");
        assert!((result - 1.5).abs() < f64::EPSILON);
    }

    #[test]
    fn wall_shear_without_crossing_is_rejected() {
        let samples = [
            WallShearSample {
                x: 1.0,
                shear_rate: 1.0,
            },
            WallShearSample {
                x: 2.0,
                shear_rate: 2.0,
            },
        ];
        let error = interpolate_reattachment(&samples).expect_err("crossing is absent");
        assert_eq!(
            error.to_string(),
            "Invalid input: backward-facing-step field has no downstream wall-shear reattachment crossing"
        );
    }

    #[test]
    fn geometry_validation_rejects_step_taller_than_channel() {
        let error = BackwardFacingStepGeometry::new(2.0_f64, 1.0, 1.0, 1.0, 1.0)
            .expect_err("invalid expansion geometry");
        assert!(error
            .to_string()
            .contains("channel height must exceed step height"));
    }

    #[test]
    fn geometry_mask_contains_step_and_keeps_downstream_floor_fluid() {
        let geometry =
            BackwardFacingStepGeometry::new(1.0_f64, 2.0, 3.0, 7.0, 1.0).expect("valid geometry");
        let mut solver = NavierStokesSolver2D::new(
            StaggeredGrid2D::new(30, 8, 10.0, 2.0),
            BloodModel::Newtonian(1.0),
            1.0,
            SIMPLEConfig::default(),
        );
        BackwardFacingStepSolver::new(geometry).apply_geometry_mask(&mut solver);

        assert!(!solver.field.mask[(0, 0)]);
        assert!(solver.field.mask[(10, 0)]);
        assert!(solver.field.mask[(0, 7)]);
    }
}
