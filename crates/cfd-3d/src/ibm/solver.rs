#![allow(
    clippy::field_reassign_with_default,
    clippy::uninlined_format_args,
    clippy::explicit_iter_loop,
    clippy::cast_lossless
)]
//! # Immersed Boundary Method (IBM) Solver for 3D Complex Geometries
//!
//! This module implements the Immersed Boundary Method for simulating fluid flow
//! around complex geometries without body-fitted meshes.
//!
//! ## Mathematical Foundation
//!
//! ### Direct Forcing Theorem
//!
//! **Statement**: The direct forcing method enforces no-slip boundary conditions by
//! adding momentum source terms to the Navier-Stokes equations at immersed boundary points.
//!
//! **Mathematical Formulation**:
//!
//! ```math
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f
//! ```
//!
//! where the forcing term f is defined as:
//!
//! ```math
//! f(x,t) = ∫_Γ F(s,t) δ(x - X(s,t)) ds
//! ```
//!
//! **Discrete Implementation**: The forcing is applied as:
//!
//! ```math
//! f_i = (U_desired - u_i) / Δt    for points near boundary
//! ```
//!
//! **Convergence Properties**:
//! - **First-order accuracy** in boundary condition enforcement
//! - **Stability**: CFL condition must be satisfied for explicit forcing
//! - **Conservation**: Momentum conservation maintained through force balance
//!
//! **Assumptions**:
//! 1. **Linear forcing response** within each time step
//! 2. **Localized forcing** effects decay away from boundary
//! 3. **Small boundary displacement** compared to grid spacing
//!
//! **Literature**: Fadlun, E.A., Verzicco, R., Orlandi, P., Mohd-Yusof, J. (2000).
//! "Combined immersed-boundary finite-difference methods for three-dimensional
//! complex flow simulations". Journal of Computational Physics, 161(1), 35-60.
//!
//! ### Feedback Forcing Theorem
//!
//! **Statement**: The feedback forcing method uses a proportional-integral-derivative
//! (PID) controller to enforce boundary conditions through corrective forcing.
//!
//! **Mathematical Formulation**: The forcing term is computed as:
//!
//! ```math
//! f^{n+1} = K_p (U_desired - u^n) + K_i ∫(U_desired - u) dt + K_d d(U_desired - u)/dt
//! ```
//!
//! **PID Controller Parameters**:
//! - **Proportional Gain (K_p)**: Immediate correction based on current error
//! - **Integral Gain (K_i)**: Correction based on accumulated error over time
//! - **Derivative Gain (K_d)**: Anticipatory correction based on error trend
//!
//! **Stability Analysis**: The controller is stable if:
//!
//! ```math
//! |K_p| < 2/Δt,    K_i > 0,    K_d < Δt/2
//! ```
//!
//! **Convergence Properties**:
//! - **Exponential convergence** to desired boundary conditions
//! - **Robustness**: Less sensitive to grid resolution than direct forcing
//! - **Accuracy**: Higher-order boundary condition enforcement possible
//!
//! **Literature**: Goldstein, D., Handler, R., Sirovich, L. (1993).
//! "Modeling a no-slip flow boundary with an external force field".
//! Journal of Computational Physics, 105(2), 354-366.
//!
//! ### Interpolation Theory (Discrete Delta Functions)
//!
//! **Statement**: The discrete delta function provides a smooth interpolation kernel
//! for transferring quantities between Eulerian grid and Lagrangian boundary points.
//!
//! **Mathematical Formulation**: The discrete delta function D_h(x) satisfies:
//!
//! ```math
//! ∫_Ω D_h(x) dx = 1,    D_h(x) ≥ 0,    supp(D_h) ⊂ [-2h, 2h]^d
//! ```
//!
//! **Roma-Peskin Kernel**: The 4-point Roma kernel is defined as:
//!
//! ```math
//! D_h(r) = (1/h^d) ∏_{i=1}^d φ(r_i/h)
//! ```
//!
//! where φ is the 1D kernel:
//!
//! ```math
//! φ(r) = (1/8)(3 - 2|r| + √(1 + 4|r| - 4r²))    for |r| ≤ 1
//! φ(r) = (1/8)(5 - 2|r| - √(-7 + 12|r| - 4r²))  for 1 < |r| ≤ 2
//! φ(r) = 0                                       for |r| > 2
//! ```
//!
//! **Interpolation Properties**:
//! - **Conservation**: ∫ f δ_h(x - X) dx = f(X) (interpolation)
//! - **Adjoint Property**: ∫ f δ_h(x - X) dx = f(X) (spreading)
//! - **Accuracy**: Second-order accurate for smooth functions
//!
//! **Literature**: Roma, A.M., Peskin, C.S., Berger, M.J. (1999).
//! "An adaptive version of the immersed boundary method".
//! Journal of Computational Physics, 153(2), 509-534.
//!
//! Peskin, C.S. (2002). "The immersed boundary method".
//! Acta Numerica, 11, 479-517.
//!
//! ### Algorithm Implementation
//!
//! 1. **Boundary Representation**: Lagrangian points define immersed geometry
//! 2. **Force Computation**: Boundary forces computed from desired conditions
//! 3. **Force Spreading**: Forces spread to Eulerian grid using δ_h function
//! 4. **Flow Solution**: Navier-Stokes equations solved with forcing terms
//! 5. **Velocity Interpolation**: Velocities interpolated back to boundary points
//! 6. **Condition Enforcement**: Boundary conditions enforced via forcing update
//!
//! ### Convergence and Stability
//!
//! **Theorem (IBM Convergence)**: For sufficiently resolved boundaries (h << geometric features),
//! the IBM converges to the body-fitted solution with order O(h²) accuracy.
//!
//! **Stability Condition**: The forcing time scale must satisfy τ_force << Δt_flow
//! to maintain numerical stability.
//!
//! **Literature**: Mittal, R., Iaccarino, G. (2005). "Immersed boundary methods".
//! Annual Review of Fluid Mechanics, 37, 239-261.

use super::{
    config::IbmConfig,
    forcing::{DirectForcing, FeedbackForcing, ForcingMethod},
    interpolation::{DeltaFunction, InterpolationKernel},
    lagrangian::LagrangianPoint,
};
use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::FloatElement;
use eunomia::NumericElement;
use eunomia::RealField;
use leto::Vector3;

// Feedback control constants
const DEFAULT_PROPORTIONAL_GAIN: f64 = 10.0;
const DEFAULT_INTEGRAL_GAIN: f64 = 1.0;

fn validate_finite_positive_ibm<T: NumericElement>(name: &str, value: T) -> Result<()> {
    if !<T as NumericElement>::is_finite(value) || value <= <T as NumericElement>::ZERO {
        return Err(Error::InvalidConfiguration(format!(
            "IbmSolver: {name} must be finite and strictly positive, got {value:?}"
        )));
    }
    Ok(())
}

/// IBM solver for 3D flow around immersed boundaries
pub struct IbmSolver<T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy> {
    /// Configuration
    config: IbmConfig,
    /// Lagrangian points representing the immersed boundary
    lagrangian_points: Vec<LagrangianPoint<T>>,
    /// Interpolation kernel
    kernel: InterpolationKernel<T>,
    /// Forcing method
    forcing: Box<dyn ForcingMethod<T>>,
    /// Grid spacing
    dx: Vector3<T>,
    /// Grid dimensions
    grid_size: (usize, usize, usize),
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy> IbmSolver<T> {
    /// Create a new IBM solver
    ///
    /// Prefer [`IbmSolver::try_new`] for fallible construction; this
    /// constructor is retained for callers that already validated their
    /// inputs and need the infallible signature.
    pub fn new(config: IbmConfig, dx: Vector3<T>, grid_size: (usize, usize, usize)) -> Self {
        match Self::try_new(config, dx, grid_size) {
            Ok(solver) => solver,
            Err(error) => panic!("IbmSolver::new called with invalid inputs: {error}"),
        }
    }

    /// Create a new IBM solver, validating the physical inputs.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when:
    ///
    /// - `dx.x`, `dx.y`, or `dx.z` is non-finite or non-positive (the
    ///   discrete-delta interpolation divides by `dx`),
    /// - any of `grid_size.0/1/2` is zero (the stencil cannot iterate an
    ///   empty grid),
    /// - `config.smoothing_width` is non-finite or non-positive (the kernel
    ///   normalization is `1/h^d`).
    pub fn try_new(
        config: IbmConfig,
        dx: Vector3<T>,
        grid_size: (usize, usize, usize),
    ) -> Result<Self> {
        validate_finite_positive_ibm("dx.x", dx.x)?;
        validate_finite_positive_ibm("dx.y", dx.y)?;
        validate_finite_positive_ibm("dx.z", dx.z)?;
        if !<f64 as NumericElement>::is_finite(config.smoothing_width)
            || config.smoothing_width <= 0.0
        {
            return Err(Error::InvalidConfiguration(format!(
                "IbmSolver::try_new: config.smoothing_width must be finite and positive, got {}",
                config.smoothing_width
            )));
        }
        if grid_size.0 == 0 || grid_size.1 == 0 || grid_size.2 == 0 {
            return Err(Error::InvalidConfiguration(format!(
                "IbmSolver::try_new: grid_size must have a positive extent in each dimension, got {grid_size:?}"
            )));
        }

        let kernel = InterpolationKernel::new(
            DeltaFunction::RomaPeskin4,
            <T as FloatElement>::from_f64(config.smoothing_width),
        );

        let forcing: Box<dyn ForcingMethod<T>> = if config.use_direct_forcing {
            Box::new(DirectForcing::new())
        } else {
            Box::new(FeedbackForcing::new(
                <T as FloatElement>::from_f64(DEFAULT_PROPORTIONAL_GAIN),
                <T as FloatElement>::from_f64(DEFAULT_INTEGRAL_GAIN),
            ))
        };

        Ok(Self {
            config,
            lagrangian_points: Vec::new(),
            kernel,
            forcing,
            dx,
            grid_size,
        })
    }

    /// Configuration used to create this solver.
    pub fn config(&self) -> &IbmConfig {
        &self.config
    }

    /// Add a Lagrangian point, validating its position is inside the grid
    /// extent (a marker outside the grid silently produces zero contribution
    /// in the stencil iteration; we now reject this at insertion time).
    pub fn add_lagrangian_point(&mut self, point: LagrangianPoint<T>) {
        match self.try_add_lagrangian_point(point) {
            Ok(()) => {}
            Err(error) => panic!("IbmSolver::add_lagrangian_point rejected the point: {error}"),
        }
    }

    /// Add a Lagrangian point, returning an error if the position is outside
    /// the grid extent or non-finite.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when the point's position is
    /// non-finite, or when its position lies outside the grid extent (the
    /// `cfd-3d` AGENTS.md contract requires IBM markers to be at least
    /// `1.5Δx` from grid boundaries).
    pub fn try_add_lagrangian_point(&mut self, point: LagrangianPoint<T>) -> Result<()> {
        let p = &point.position;
        if !<T as NumericElement>::is_finite(p.x)
            || !<T as NumericElement>::is_finite(p.y)
            || !<T as NumericElement>::is_finite(p.z)
        {
            return Err(Error::InvalidConfiguration(format!(
                "IbmSolver::try_add_lagrangian_point: point position must be finite, got {p:?}"
            )));
        }
        let max_x = <T as FloatElement>::from_f64(self.grid_size.0 as f64) * self.dx.x;
        let max_y = <T as FloatElement>::from_f64(self.grid_size.1 as f64) * self.dx.y;
        let max_z = <T as FloatElement>::from_f64(self.grid_size.2 as f64) * self.dx.z;
        let lower = <T as NumericElement>::ZERO;
        if p.x < lower || p.x >= max_x || p.y < lower || p.y >= max_y || p.z < lower || p.z >= max_z
        {
            return Err(Error::InvalidConfiguration(format!(
                "IbmSolver::try_add_lagrangian_point: point position {p:?} must lie inside the grid extent [0, {max_x:?}) x [0, {max_y:?}) x [0, {max_z:?})"
            )));
        }
        self.lagrangian_points.push(point);
        Ok(())
    }

    /// Interpolate velocity from Eulerian to Lagrangian grid
    pub fn interpolate_velocity(&mut self, eulerian_velocity: &[Vector3<T>]) -> Result<()> {
        // Collect velocities first to avoid borrow conflict
        let velocities: Vec<_> = self
            .lagrangian_points
            .iter()
            .map(|point| self.interpolate_at_point(&point.position, eulerian_velocity))
            .collect::<Result<Vec<_>>>()?;

        // Then update the points
        for (point, velocity) in self.lagrangian_points.iter_mut().zip(velocities) {
            point.velocity = velocity;
        }
        Ok(())
    }

    /// Calculate forcing terms
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when `dt` is non-finite or
    /// non-positive (direct forcing divides by `dt`, feedback forcing uses
    /// `dt` in the derivative term).
    pub fn calculate_forcing(&mut self, desired_velocity: &[Vector3<T>], dt: T) -> Result<()> {
        validate_finite_positive_ibm("dt", dt)?;
        for (i, point) in self.lagrangian_points.iter_mut().enumerate() {
            point.force = self
                .forcing
                .calculate_force(&desired_velocity[i], &point.velocity, dt);
        }
        Ok(())
    }

    /// Spread forces from Lagrangian to Eulerian grid
    pub fn spread_forces(&self) -> Result<Vec<Vector3<T>>> {
        let total_points = self.grid_size.0 * self.grid_size.1 * self.grid_size.2;
        let mut eulerian_forces = vec![Vector3::zeros(); total_points];

        for point in &self.lagrangian_points {
            self.spread_from_point(point, &mut eulerian_forces)?;
        }

        Ok(eulerian_forces)
    }

    /// Interpolate at a specific point
    fn interpolate_at_point(
        &self,
        position: &Vector3<T>,
        field: &[Vector3<T>],
    ) -> Result<Vector3<T>> {
        let mut result = Vector3::zeros();
        let stencil = self.kernel.stencil_size();

        // Find grid indices and convert to integers
        let i_int = Self::floored_grid_index(position.x / self.dx.x, "x")?;
        let j_int = Self::floored_grid_index(position.y / self.dx.y, "y")?;
        let k_int = Self::floored_grid_index(position.z / self.dx.z, "z")?;

        let i_start = i_int - (stencil as isize / 2);
        let j_start = j_int - (stencil as isize / 2);
        let k_start = k_int - (stencil as isize / 2);

        for di in 0..stencil {
            for dj in 0..stencil {
                for dk in 0..stencil {
                    let ii = (i_start + di as isize) as usize;
                    let jj = (j_start + dj as isize) as usize;
                    let kk = (k_start + dk as isize) as usize;

                    if ii < self.grid_size.0 && jj < self.grid_size.1 && kk < self.grid_size.2 {
                        let idx =
                            ii + jj * self.grid_size.0 + kk * self.grid_size.0 * self.grid_size.1;

                        let rx = position.x / self.dx.x - scalar::from_usize::<T>(ii);
                        let ry = position.y / self.dx.y - scalar::from_usize::<T>(jj);
                        let rz = position.z / self.dx.z - scalar::from_usize::<T>(kk);

                        let weight =
                            self.kernel.delta(rx) * self.kernel.delta(ry) * self.kernel.delta(rz);
                        result += field[idx] * weight;
                    }
                }
            }
        }

        Ok(result)
    }

    /// Spread from a Lagrangian point to Eulerian grid
    fn spread_from_point(
        &self,
        point: &LagrangianPoint<T>,
        field: &mut [Vector3<T>],
    ) -> Result<()> {
        let stencil = self.kernel.stencil_size();

        // Find grid indices and convert to integers
        let i_int = Self::floored_grid_index(point.position.x / self.dx.x, "x")?;
        let j_int = Self::floored_grid_index(point.position.y / self.dx.y, "y")?;
        let k_int = Self::floored_grid_index(point.position.z / self.dx.z, "z")?;

        let i_start = i_int - (stencil as isize / 2);
        let j_start = j_int - (stencil as isize / 2);
        let k_start = k_int - (stencil as isize / 2);

        for di in 0..stencil {
            for dj in 0..stencil {
                for dk in 0..stencil {
                    let ii = (i_start + di as isize) as usize;
                    let jj = (j_start + dj as isize) as usize;
                    let kk = (k_start + dk as isize) as usize;

                    if ii < self.grid_size.0 && jj < self.grid_size.1 && kk < self.grid_size.2 {
                        let idx =
                            ii + jj * self.grid_size.0 + kk * self.grid_size.0 * self.grid_size.1;

                        let rx = point.position.x / self.dx.x - scalar::from_usize::<T>(ii);
                        let ry = point.position.y / self.dx.y - scalar::from_usize::<T>(jj);
                        let rz = point.position.z / self.dx.z - scalar::from_usize::<T>(kk);

                        let weight =
                            self.kernel.delta(rx) * self.kernel.delta(ry) * self.kernel.delta(rz);
                        field[idx] += point.force * weight * point.weight;
                    }
                }
            }
        }

        Ok(())
    }

    /// Update Lagrangian point positions
    pub fn update_positions(&mut self, dt: T) {
        for point in &mut self.lagrangian_points {
            point.update_position(dt);
        }
    }

    /// Get number of Lagrangian points
    pub fn num_points(&self) -> usize {
        self.lagrangian_points.len()
    }

    /// Get reference to Lagrangian points
    pub fn points(&self) -> &[LagrangianPoint<T>] {
        &self.lagrangian_points
    }

    /// Get mutable reference to Lagrangian points
    pub fn points_mut(&mut self) -> &mut [LagrangianPoint<T>] {
        &mut self.lagrangian_points
    }

    fn floored_grid_index(value: T, axis: &str) -> Result<isize> {
        let floored = scalar::floor::<T>(value);
        let index = <T as NumericElement>::to_f64(floored);
        if !index.is_finite() || index < isize::MIN as f64 || index > isize::MAX as f64 {
            return Err(Error::InvalidConfiguration(format!(
                "IBM grid {axis}-coordinate must be finite and within isize range after scaling"
            )));
        }
        Ok(index as isize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_solver() -> IbmSolver<f64> {
        let config = IbmConfig::default();
        let dx = Vector3::new(0.1, 0.1, 0.1);
        let grid_size = (10, 10, 10);
        IbmSolver::new(config, dx, grid_size)
    }

    #[test]
    fn creation_with_valid_config() {
        let solver = default_solver();
        assert_eq!(solver.num_points(), 0);
        assert!(solver.config().use_direct_forcing);
        assert_eq!(solver.grid_size, (10, 10, 10));
    }

    #[test]
    fn add_lagrangian_point_increments_count() {
        let mut solver = default_solver();
        let pt = LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 1.0);
        solver.add_lagrangian_point(pt);
        assert_eq!(solver.num_points(), 1);
    }

    #[test]
    fn spread_forces_produces_finite_values() {
        let mut solver = default_solver();

        // Place a point near the centre of the grid with a known force.
        let mut pt = LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 1.0);
        pt.force = Vector3::new(1.0, 0.0, 0.0);
        solver.add_lagrangian_point(pt);

        let forces = solver
            .spread_forces()
            .expect("spread_forces should succeed");

        assert_eq!(forces.len(), 10 * 10 * 10);

        // At least one grid cell should have received a non-zero, finite force.
        let has_nonzero = forces
            .iter()
            .any(|f| f.norm() > 0.0 && f.norm().is_finite());
        assert!(
            has_nonzero,
            "force spreading should produce at least one non-zero finite value on the grid"
        );

        // All values must be finite.
        for (i, f) in forces.iter().enumerate() {
            assert!(
                f.x.is_finite() && f.y.is_finite() && f.z.is_finite(),
                "non-finite force at grid index {}: {:?}",
                i,
                f
            );
        }
    }

    #[test]
    fn try_new_rejects_zero_dx() {
        let config = IbmConfig::default();
        let dx = Vector3::new(0.0, 0.1, 0.1);
        let result = IbmSolver::try_new(config, dx, (10, 10, 10));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("dx.x")),
            "expected InvalidConfiguration error for zero dx.x"
        );
    }

    #[test]
    fn try_new_rejects_non_finite_dx() {
        let config = IbmConfig::default();
        let dx = Vector3::new(f64::NAN, 0.1, 0.1);
        let result = IbmSolver::try_new(config, dx, (10, 10, 10));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("dx.x")),
            "expected InvalidConfiguration error for NaN dx.x"
        );
    }

    #[test]
    fn try_new_rejects_zero_grid_extent() {
        let config = IbmConfig::default();
        let dx = Vector3::new(0.1, 0.1, 0.1);
        let result = IbmSolver::try_new(config, dx, (0, 10, 10));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("grid_size")),
            "expected InvalidConfiguration error for zero grid_size.0"
        );
    }

    #[test]
    fn try_new_rejects_non_positive_smoothing_width() {
        let config = IbmConfig {
            smoothing_width: -1.5,
            ..IbmConfig::default()
        };
        let dx = Vector3::new(0.1, 0.1, 0.1);
        let result = IbmSolver::try_new(config, dx, (10, 10, 10));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("smoothing_width")),
            "expected InvalidConfiguration error for negative smoothing_width"
        );
    }

    #[test]
    fn try_new_accepts_physically_valid_inputs() {
        let config = IbmConfig::default();
        let dx = Vector3::new(0.1, 0.1, 0.1);
        let result = IbmSolver::try_new(config, dx, (10, 10, 10));
        assert!(
            result.is_ok(),
            "valid config, dx, and grid_size must succeed"
        );
    }

    #[test]
    fn try_add_lagrangian_point_rejects_out_of_grid_position() {
        let mut solver = default_solver();
        let out_of_grid = LagrangianPoint::new(Vector3::new(100.0, 0.5, 0.5), 1.0);
        let result = solver.try_add_lagrangian_point(out_of_grid);
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("inside the grid")),
            "expected InvalidConfiguration error for out-of-grid point"
        );
    }

    #[test]
    fn try_add_lagrangian_point_rejects_non_finite_position() {
        let mut solver = default_solver();
        let bad = LagrangianPoint::new(Vector3::new(f64::NAN, 0.5, 0.5), 1.0);
        let result = solver.try_add_lagrangian_point(bad);
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("finite")),
            "expected InvalidConfiguration error for non-finite point position"
        );
    }

    #[test]
    fn try_add_lagrangian_point_accepts_in_grid_position() {
        let mut solver = default_solver();
        let inside = LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 1.0);
        let result = solver.try_add_lagrangian_point(inside);
        assert!(result.is_ok(), "in-grid point insertion must succeed");
        assert_eq!(solver.num_points(), 1);
    }

    #[test]
    fn calculate_forcing_rejects_zero_dt() {
        let mut solver = default_solver();
        solver.add_lagrangian_point(LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 1.0));
        let desired = vec![Vector3::new(1.0, 0.0, 0.0)];
        let result = solver.calculate_forcing(&desired, 0.0);
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("dt")),
            "expected InvalidConfiguration error for zero dt"
        );
    }

    #[test]
    fn calculate_forcing_rejects_negative_dt() {
        let mut solver = default_solver();
        solver.add_lagrangian_point(LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 1.0));
        let desired = vec![Vector3::new(1.0, 0.0, 0.0)];
        let result = solver.calculate_forcing(&desired, -0.001);
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("dt")),
            "expected InvalidConfiguration error for negative dt"
        );
    }

    #[test]
    fn calculate_forcing_rejects_non_finite_dt() {
        let mut solver = default_solver();
        solver.add_lagrangian_point(LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 1.0));
        let desired = vec![Vector3::new(1.0, 0.0, 0.0)];
        let result = solver.calculate_forcing(&desired, f64::NAN);
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("dt")),
            "expected InvalidConfiguration error for NaN dt"
        );
    }
}
