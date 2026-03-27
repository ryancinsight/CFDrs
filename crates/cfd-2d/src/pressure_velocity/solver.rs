//! Main pressure-velocity coupling solver implementation
//!
//! # Theorem (SIMPLE Outer Loop Convergence)
//!
//! The coupled momentum–pressure iteration converges when the velocity and
//! pressure relaxation factors satisfy $0 < \alpha_u, \alpha_p < 1$ with
//! $\alpha_u + \alpha_p \le 1$. Each outer iteration reduces the global
//! continuity residual. See [`super`] module docs for the full proof.

use super::{PressureCorrectionSolver, PressureVelocityConfig, RhieChowInterpolation};
use crate::fields::SimulationFields;
use crate::grid::array2d::Array2D;
use crate::grid::StructuredGrid2D;
use crate::physics::{MomentumComponent, MomentumSolver};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use std::fmt::LowerExp;

/// SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) solver
/// Implementation follows Patankar (1980) "Numerical Heat Transfer and Fluid Flow"
pub struct PressureVelocitySolver<T: RealField + Copy> {
    /// Configuration
    config: PressureVelocityConfig<T>,
    /// Grid
    grid: StructuredGrid2D<T>,
    /// Momentum solver
    momentum_solver: MomentumSolver<T>,
    /// Pressure correction solver
    pressure_solver: PressureCorrectionSolver<T>,
    /// Rhie-Chow interpolation (optional)
    rhie_chow: Option<RhieChowInterpolation<T>>,
    /// Current velocity field
    u: Array2D<Vector2<T>>,
    /// Current pressure field
    p: Array2D<T>,
    /// Iteration counter
    iterations: usize,
}

impl<T: RealField + Copy + FromPrimitive + Copy + LowerExp + num_traits::ToPrimitive>
    PressureVelocitySolver<T>
{
    /// Create new pressure-velocity coupling solver
    pub fn new(
        grid: StructuredGrid2D<T>,
        config: PressureVelocityConfig<T>,
    ) -> cfd_core::error::Result<Self> {
        config.validate()?;

        let nx = grid.nx;
        let ny = grid.ny;

        let momentum_solver = MomentumSolver::new(&grid);
        let pressure_solver =
            PressureCorrectionSolver::new(grid.clone(), config.pressure_linear_solver)?;

        let rhie_chow = if config.use_rhie_chow {
            Some(RhieChowInterpolation::new(&grid))
        } else {
            None
        };

        Ok(Self {
            config,
            grid,
            momentum_solver,
            pressure_solver,
            rhie_chow,
            u: Array2D::new(nx, ny, Vector2::zeros()),
            p: Array2D::new(nx, ny, T::zero()),
            iterations: 0,
        })
    }

    /// Whether Rhie-Chow interpolation is enabled
    #[must_use]
    pub const fn has_rhie_chow(&self) -> bool {
        self.rhie_chow.is_some()
    }

    /// Set initial conditions
    pub fn set_initial_conditions(&mut self, u_init: Array2D<Vector2<T>>, p_init: Array2D<T>) {
        self.u = u_init;
        self.p = p_init;
        self.iterations = 0;
    }

    /// Perform one pressure-velocity coupling iteration
    pub fn step(
        &mut self,
        _bc: &cfd_core::physics::boundary::BoundaryCondition<T>,
        _nu: T,
        rho: T,
    ) -> cfd_core::error::Result<T> {
        // Step 1: Solve momentum equations for predicted velocity
        // Note: This needs proper integration with SimulationFields
        // Benchmark framework: Computational solver integration pending.
        // Current implementation creates fields structure for future solver hookup.
        let mut state_buffer = SimulationFields::new(self.grid.nx, self.grid.ny);
        // Copy current state to fields
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                state_buffer.set_velocity_at(i, j, &self.u[(i, j)]);
                if let Some(p) = state_buffer.p.at_mut(i, j) {
                    *p = self.p[(i, j)];
                }
            }
        }

        // Solve momentum equations (modifies state_buffer in place)
        self.momentum_solver
            .solve(MomentumComponent::U, &mut state_buffer, self.config.dt)?;
        self.momentum_solver
            .solve(MomentumComponent::V, &mut state_buffer, self.config.dt)?;

        // Extract predicted velocity field
        let mut u_star = Array2D::new(self.grid.nx, self.grid.ny, Vector2::zeros());
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                u_star[(i, j)] = Vector2::new(state_buffer.u.at(i, j), state_buffer.v.at(i, j));
            }
        }

        // Step 2: Solve pressure correction equation
        let p_correction =
            self.pressure_solver
                .solve_pressure_correction(&state_buffer, self.config.dt, rho)?;

        // Step 3: Correct velocity field
        let mut u_corrected = u_star;
        let (ap_u, _, ap_v, _) = self.momentum_solver.get_ap_coefficients();
        self.pressure_solver.correct_velocity(
            &mut u_corrected,
            &p_correction,
            ap_u,
            ap_v,
            rho,
            self.config.alpha_u,
            &state_buffer, // Pass buffer for mask access
        );

        // Step 4: Correct pressure field
        self.pressure_solver
            .correct_pressure(&mut self.p, &p_correction, self.config.alpha_p);

        // Step 5: Calculate residual for convergence check
        let residual = self.calculate_residual(&self.u, &u_corrected);

        // Update fields
        self.u = u_corrected;
        self.iterations += 1;

        Ok(residual)
    }

    /// Run solver until convergence or max iterations
    pub fn solve(
        &mut self,
        bc: &cfd_core::physics::boundary::BoundaryCondition<T>,
        nu: T,
        rho: T,
    ) -> cfd_core::error::Result<()> {
        let max_iter = self.config.base.convergence.max_iterations;
        let tolerance = self.config.base.convergence.tolerance;

        for _ in 0..max_iter {
            let residual = self.step(bc, nu, rho)?;

            if residual < tolerance {
                tracing::info!(
                    "Pressure-velocity coupling converged in {} iterations (residual: {:e})",
                    self.iterations,
                    residual
                );
                return Ok(());
            }
        }

        Err(cfd_core::error::Error::Convergence(
            cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { max: max_iter },
        ))
    }

    /// Calculate residual between two velocity fields
    fn calculate_residual(&self, previous: &Array2D<Vector2<T>>, current: &Array2D<Vector2<T>>) -> T {
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        let mut max_diff = T::zero();
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let diff = (current[(i, j)] - previous[(i, j)]).norm();
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }

        max_diff
    }

    /// Get current velocity field
    pub fn velocity(&self) -> &Array2D<Vector2<T>> {
        &self.u
    }

    /// Get current pressure field
    pub fn pressure(&self) -> &Array2D<T> {
        &self.p
    }

    /// Get iteration count
    pub fn iterations(&self) -> usize {
        self.iterations
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::StructuredGrid2D;

    fn make_solver(nx: usize, ny: usize) -> PressureVelocitySolver<f64> {
        let grid = StructuredGrid2D::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
        let config = PressureVelocityConfig::new().unwrap();
        PressureVelocitySolver::new(grid, config).unwrap()
    }

    #[test]
    fn solver_creation_with_valid_configuration() {
        let solver = make_solver(8, 8);
        assert_eq!(solver.iterations(), 0);
        assert!(!solver.velocity().as_slice().is_empty());
        assert!(!solver.pressure().as_slice().is_empty());
    }

    #[test]
    fn config_propagation_nx_ny() {
        let grid = StructuredGrid2D::new(10, 6, 0.0, 2.0, 0.0, 1.0).unwrap();
        let config = PressureVelocityConfig::new().unwrap();
        let solver = PressureVelocitySolver::new(grid, config).unwrap();

        // Velocity field dimensions should match grid
        assert_eq!(solver.velocity().rows(), 10);
        assert_eq!(solver.velocity().cols(), 6);
        assert_eq!(solver.pressure().rows(), 10);
        assert_eq!(solver.pressure().cols(), 6);
    }

    #[test]
    fn config_propagation_dx_dy() {
        let grid = StructuredGrid2D::new(5, 5, 0.0, 4.0, 0.0, 2.0).unwrap();
        let dx: f64 = grid.dx;
        let dy: f64 = grid.dy;
        let config = PressureVelocityConfig::new().unwrap();
        let solver = PressureVelocitySolver::new(grid, config).unwrap();

        assert!((solver.grid.dx - dx).abs() < 1e-15_f64);
        assert!((solver.grid.dy - dy).abs() < 1e-15_f64);
    }

    #[test]
    fn set_initial_conditions_resets_iteration_count() {
        let mut solver = make_solver(4, 4);
        let u_init = Array2D::new(4, 4, Vector2::new(1.0, 0.0));
        let p_init = Array2D::new(4, 4, 0.0);
        solver.set_initial_conditions(u_init, p_init);
        assert_eq!(solver.iterations(), 0);
    }
}
