//! Discretized Venturi solver and validation against analytical solutions.

use super::{BernoulliVenturi, VenturiFlowSolution, VenturiGeometry};
use crate::solvers::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Result as CfdResult;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Discretized Venturi Solver
// ============================================================================

/// 2D Venturi flow solver using Finite Volume Method (FVM)
pub struct VenturiSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    _geometry: VenturiGeometry<T>,
    solver: NavierStokesSolver2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> VenturiSolver2D<T> {
    /// Create a new discretized Venturi solver with uniform grid spacing.
    pub fn new(
        geometry: VenturiGeometry<T>,
        blood: BloodModel<T>,
        density: T,
        nx: usize,
        ny: usize,
    ) -> Self {
        let grid = StaggeredGrid2D::new(nx, ny, geometry.total_length(), geometry.w_inlet);
        let config = SIMPLEConfig::default();
        let mut solver = NavierStokesSolver2D::new(grid, blood, density, config);
        Self::populate_mask(&mut solver, &geometry, nx, ny);
        Self {
            _geometry: geometry,
            solver,
        }
    }

    /// Create a Venturi solver with non-uniform y-spacing that clusters cells
    /// near the channel centre (where the throat walls are).
    ///
    /// Uses a sine-based algebraic stretching:
    /// ```text
    /// y_j = ly · (η + β·sin(2πη) / (2π)),   η = j / ny
    /// ```
    ///
    /// The parameter `beta` ∈ (0, 1) controls the centre-clustering intensity:
    /// - `beta = 0.0` → uniform grid (same as [`new`])
    /// - `beta = 0.5` → 3× ratio (centre cells 3× finer than boundary)
    /// - `beta = 0.9` → 19× ratio (aggressive, for CR > 20)
    ///
    /// A useful heuristic: `beta = min(1 − 4·w_throat/w_inlet, 0.9)`.
    ///
    /// # Panics
    /// Panics if `beta >= 1.0` (cells would overlap) or `beta < 0.0`.
    pub fn new_stretched(
        geometry: VenturiGeometry<T>,
        blood: BloodModel<T>,
        density: T,
        nx: usize,
        ny: usize,
        beta: T,
    ) -> Self {
        Self::new_stretched_with_config(geometry, blood, density, nx, ny, beta, SIMPLEConfig::default())
    }

    /// Like [`new_stretched`] but accepts a custom [`SIMPLEConfig`], allowing
    /// callers to tune under-relaxation and iteration count for high-Re flows.
    pub fn new_stretched_with_config(
        geometry: VenturiGeometry<T>,
        blood: BloodModel<T>,
        density: T,
        nx: usize,
        ny: usize,
        beta: T,
        simple_config: SIMPLEConfig<T>,
    ) -> Self {
        assert!(
            beta >= T::zero() && beta < T::one(),
            "beta must be in [0, 1)"
        );

        let ly = geometry.w_inlet;
        let two_pi = T::from_f64(2.0 * std::f64::consts::PI).unwrap_or_else(num_traits::Zero::zero);
        let mut y_faces = Vec::with_capacity(ny + 1);
        let ny_t = T::from_usize(ny).unwrap_or_else(T::one);
        for j in 0..=ny {
            let eta = T::from_usize(j).unwrap_or_else(T::one) / ny_t;
            let s = eta + beta * Float::sin(two_pi * eta) / two_pi;
            y_faces.push(ly * s);
        }

        let lx = geometry.total_length();
        let grid = StaggeredGrid2D::new_stretched_y(nx, ny, lx, y_faces);
        let mut solver = NavierStokesSolver2D::new(grid, blood, density, simple_config);
        Self::populate_mask(&mut solver, &geometry, nx, ny);

        // Log resolution at throat for debugging
        let dy_min = (0..ny).map(|j| solver.grid.dy_at(j)).fold(ly, Float::min);
        let cr = geometry.w_inlet
            / Float::max(
                geometry.w_throat,
                T::from_f64(1e-12).unwrap_or_else(num_traits::Zero::zero),
            );
        if dy_min > geometry.w_throat / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero) {
            tracing::debug!(
                "[VenturiSolver2D] WARNING: dy_min ({:.2e}) > w_throat/2 ({:.2e}). \
                 CR={:.1}. Consider increasing ny or beta.",
                dy_min.to_f64().unwrap_or(0.0),
                (geometry.w_throat / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero))
                    .to_f64()
                    .unwrap_or(0.0),
                cr.to_f64().unwrap_or(0.0),
            );
        }

        Self {
            _geometry: geometry,
            solver,
        }
    }

    /// Populate the solver's fluid mask from the geometry.
    ///
    /// Uses `grid.x_center(i)` and `grid.y_center(j)` to handle both uniform
    /// and non-uniform grids correctly.
    fn populate_mask(
        solver: &mut NavierStokesSolver2D<T>,
        geometry: &VenturiGeometry<T>,
        nx: usize,
        ny: usize,
    ) {
        let half_h = geometry.w_inlet / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);
        for i in 0..nx {
            for j in 0..ny {
                let x = solver.grid.x_center(i);
                let y = solver.grid.y_center(j) - half_h;
                solver.field.mask[(i, j)] = geometry.contains(x, y);
            }
        }
    }

    /// Solve the Venturi flow for a given inlet velocity
    pub fn solve(&mut self, u_inlet: T) -> CfdResult<VenturiFlowSolution<T>> {
        let solve_result = self
            .solver
            .solve(u_inlet)
            .map_err(|e| cfd_core::error::Error::Solver(e.to_string()))?;

        // Extract metrics
        let nx = self.solver.grid.nx;
        let ny = self.solver.grid.ny;

        // Average inlet velocity from the first internal u-faces — zero-alloc fold.
        let (u_sum_in, count_inlet) = (0..ny)
            .filter(|&j| self.solver.field.mask[(0, j)])
            .fold((T::zero(), 0usize), |(s, n), j| {
                (s + self.solver.field.u[(1, j)], n + 1)
            });
        let u_inlet_sim = if count_inlet > 0 {
            u_sum_in / T::from_usize(count_inlet).unwrap_or_else(T::one)
        } else {
            T::zero()
        };

        // Find throat section: cell with maximum u-velocity — single flat_map fold.
        let (u_max, p_throat) = (0..nx)
            .flat_map(|i| (0..ny).map(move |j| (i, j)))
            .filter(|&(i, j)| self.solver.field.mask[(i, j)])
            .fold((T::zero(), T::zero()), |(u_best, p_best), (i, j)| {
                let u = self.solver.field.u[(i, j)];
                if u > u_best {
                    (u, self.solver.field.p[(i, j)])
                } else {
                    (u_best, p_best)
                }
            });

        // Area-averaged throat velocity: find the column at the throat midpoint
        // and average the x-velocity across all fluid cells in that column.
        //
        // Theorem (cross-section averaging): For comparison with 1D models that
        // predict mean velocity (Hagen-Poiseuille, Bernoulli continuity), the
        // area-averaged velocity ū = (1/A) ∫ u dA is the correct metric.
        // For fully-developed 2D laminar flow, ū = (2/3) u_max.
        let throat_x_mid = self._geometry.l_inlet + self._geometry.l_converge
            + self._geometry.l_throat / T::from_f64(2.0).unwrap_or_else(T::one);
        let i_throat = (0..nx)
            .min_by_key(|&i| {
                let dx = self.solver.grid.x_center(i) - throat_x_mid;
                // Convert to integer for Ord comparison (f64 doesn't implement Ord)
                (dx * dx * T::from_f64(1e12).unwrap_or_else(T::one))
                    .to_u64()
                    .unwrap_or(u64::MAX)
            })
            .unwrap_or(nx / 2);
        let (u_sum_throat, count_throat) = (0..ny)
            .filter(|&j| self.solver.field.mask[(i_throat, j)])
            .fold((T::zero(), 0usize), |(s, n), j| {
                (s + self.solver.field.u[(i_throat, j)], n + 1)
            });
        let u_throat_mean = if count_throat > 0 {
            u_sum_throat / T::from_usize(count_throat).unwrap_or_else(T::one)
        } else {
            u_max // fallback to max if no fluid cells found at throat
        };

        // Average outlet pressure from the last fluid column — zero-alloc fold.
        let (p_sum_out, count_outlet) = (0..ny)
            .filter(|&j| self.solver.field.mask[(nx - 1, j)])
            .fold((T::zero(), 0usize), |(s, n), j| {
                (s + self.solver.field.p[(nx - 1, j)], n + 1)
            });
        let p_outlet = if count_outlet > 0 {
            p_sum_out / T::from_usize(count_outlet).unwrap_or_else(T::one)
        } else {
            T::zero()
        };

        let p_inlet = self.solver.field.p[(0, ny / 2)]; // Reference inlet pressure

        let dp_throat = p_throat - p_inlet;
        let dp_recovery = p_outlet - p_inlet;

        let rho = self.solver.density;
        let q_dyn =
            T::from_f64(0.5).unwrap_or_else(num_traits::Zero::zero) * rho * u_inlet * u_inlet;
        // Guard: avoid division by near-zero dynamic pressure.  Use 1e-6 Pa (not 1.0)
        // so Cp remains physically meaningful even at very low inlet velocities.
        let q_dyn_safe = num_traits::Float::max(
            q_dyn,
            T::from_f64(1e-6).unwrap_or_else(num_traits::Zero::zero),
        );
        let cp_throat = dp_throat / q_dyn_safe;
        let cp_recovery = dp_recovery / q_dyn_safe;

        Ok(VenturiFlowSolution {
            u_inlet: u_inlet_sim,
            p_inlet,
            u_throat: u_max,
            u_throat_mean,
            p_throat,
            u_outlet: u_inlet_sim, // Assumes symmetric inlet/outlet
            p_outlet,
            dp_throat,
            dp_recovery,
            cp_throat,
            cp_recovery,
            converged: solve_result.converged,
        })
    }
}

// ============================================================================
// Validation Against Literature
// ============================================================================

/// Venturi validation against analytical and literature solutions
pub struct VenturiValidator<T: RealField + Copy> {
    geometry: VenturiGeometry<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> VenturiValidator<T> {
    /// Create new validator
    pub fn new(geometry: VenturiGeometry<T>) -> Self {
        Self { geometry }
    }

    /// Validate numerical solution against Bernoulli
    ///
    /// # Validation Criteria
    ///
    /// - Throat pressure: error < 5% (friction effects)
    /// - Outlet pressure: error < 10% (recovery losses)
    /// - Velocity magnitude: error < 1% (mass conservation)
    pub fn validate_against_bernoulli(
        &self,
        numerical: &VenturiFlowSolution<T>,
        u_inlet: T,
        p_inlet: T,
        rho: T,
    ) -> Result<VenturiValidationResult<T>, String> {
        let bernoulli = BernoulliVenturi::new(self.geometry.clone(), u_inlet, p_inlet, rho);

        let u_throat_analytical = bernoulli.velocity_throat();
        let p_throat_analytical = bernoulli.pressure_throat();

        // Calculate relative errors
        let u_throat_error =
            (numerical.u_throat - u_throat_analytical).abs() / u_throat_analytical.abs();
        let p_throat_error = (numerical.p_throat - p_throat_analytical).abs()
            / p_throat_analytical.abs().max(T::from_f64_or_one(1.0));

        let mut result = VenturiValidationResult {
            u_throat_error: Some(u_throat_error),
            p_throat_error: Some(p_throat_error),
            validation_passed: false,
            error_message: None,
        };

        // Check tolerances
        let tolerance_u = T::from_f64_or_one(0.01); // 1%
        let tolerance_p = T::from_f64_or_one(0.05); // 5%

        if u_throat_error < tolerance_u && p_throat_error < tolerance_p {
            result.validation_passed = true;
        } else {
            let mut msg = String::new();
            if u_throat_error >= tolerance_u {
                use std::fmt::Write;
                let _ = write!(
                    msg,
                    "Throat velocity error {:.2e} > 1%",
                    u_throat_error.to_f64().unwrap_or(f64::NAN)
                );
            }
            if p_throat_error >= tolerance_p {
                use std::fmt::Write;
                let _ = write!(
                    msg,
                    "Throat pressure error {:.2e} > 5%",
                    p_throat_error.to_f64().unwrap_or(f64::NAN)
                );
            }
            result.error_message = Some(msg);
        }

        Ok(result)
    }
}

/// Validation result for Venturi
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiValidationResult<T: RealField + Copy> {
    /// Relative error in throat velocity
    pub u_throat_error: Option<T>,
    /// Relative error in throat pressure
    pub p_throat_error: Option<T>,
    /// Validation passed
    pub validation_passed: bool,
    /// Error message if validation failed
    pub error_message: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_discretized_venturi_solver_convergence() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();
        let blood = BloodModel::Newtonian(1.0e-3);
        let density = 1060.0;

        let mut solver = VenturiSolver2D::new(geom, blood, density, 40, 20);
        let result = solver.solve(0.1); // 100 mm/s

        assert!(result.is_ok(), "Solver failed: {:?}", result.err());
        let sol = result.unwrap();
        println!("Venturi Solution: {:?}", sol);

        // Qualification checks
        assert!(
            sol.u_throat > sol.u_inlet,
            "Velocity should increase in throat ({:.4} > {:.4})",
            sol.u_throat,
            sol.u_inlet
        );
        assert!(
            sol.p_throat < sol.p_inlet,
            "Pressure should decrease in throat ({:.4} < {:.4})",
            sol.p_throat,
            sol.p_inlet
        );
        assert!(
            sol.cp_throat < 0.0,
            "Pressure coefficient at throat should be negative ({:.4})",
            sol.cp_throat
        );
    }

    #[test]
    fn test_stretched_venturi_solver_convergence() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();
        let blood = BloodModel::Newtonian(1.0e-3);
        let density = 1060.0;

        let mut solver = VenturiSolver2D::new_stretched(geom, blood, density, 40, 20, 0.5);
        let result = solver.solve(0.1);

        assert!(
            result.is_ok(),
            "Stretched solver failed: {:?}",
            result.err()
        );
        let sol = result.unwrap();

        assert!(
            sol.u_throat > sol.u_inlet,
            "Velocity should increase in throat ({:.4} > {:.4})",
            sol.u_throat,
            sol.u_inlet
        );
        assert!(
            sol.p_throat < sol.p_inlet,
            "Pressure should decrease in throat ({:.4} < {:.4})",
            sol.p_throat,
            sol.p_inlet
        );
    }

    #[test]
    fn test_stretched_grid_has_finer_center_cells() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();
        let blood = BloodModel::Newtonian(1.0e-3);
        let ny = 40;
        let solver = VenturiSolver2D::new_stretched(geom, blood, 1060.0, 40, ny, 0.8);
        // Centre cell should be smaller than boundary cell
        let dy_center = solver.solver.grid.dy_at(ny / 2);
        let dy_boundary = solver.solver.grid.dy_at(0);
        assert!(
            dy_center < dy_boundary,
            "Centre cells ({:.6}) should be finer than boundary cells ({:.6})",
            dy_center,
            dy_boundary
        );
    }
}
