//! Main pressure-velocity coupling solver implementation

use super::{PressureCorrectionSolver, PressureVelocityConfig, RhieChowInterpolation};
use crate::fields::SimulationFields;
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
    #[allow(dead_code)]
    rhie_chow: Option<RhieChowInterpolation<T>>,
    /// Current velocity field
    u: Vec<Vec<Vector2<T>>>,
    /// Current pressure field
    p: Vec<Vec<T>>,
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
            u: vec![vec![Vector2::zeros(); ny]; nx],
            p: vec![vec![T::zero(); ny]; nx],
            iterations: 0,
        })
    }

    /// Set initial conditions
    pub fn set_initial_conditions(&mut self, u_init: Vec<Vec<Vector2<T>>>, p_init: Vec<Vec<T>>) {
        self.u = u_init;
        self.p = p_init;
        self.iterations = 0;
    }

    /// Perform one pressure-velocity coupling iteration
    pub fn step(
        &mut self,
        _bc: &cfd_core::boundary::BoundaryCondition<T>,
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
                state_buffer.set_velocity_at(i, j, &self.u[i][j]);
                if let Some(p) = state_buffer.p.at_mut(i, j) {
                    *p = self.p[i][j];
                }
            }
        }

        // Solve momentum equations (modifies state_buffer in place)
        self.momentum_solver
            .solve(MomentumComponent::U, &mut state_buffer, self.config.dt)?;
        self.momentum_solver
            .solve(MomentumComponent::V, &mut state_buffer, self.config.dt)?;

        // Extract predicted velocity field
        let mut u_star = vec![vec![Vector2::zeros(); self.grid.ny]; self.grid.nx];
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                u_star[i][j] = Vector2::new(state_buffer.u.at(i, j), state_buffer.v.at(i, j));
            }
        }

        // Step 2: Solve pressure correction equation
        let p_correction =
            self.pressure_solver
                .solve_pressure_correction(&u_star, self.config.dt, rho)?;

        // Step 3: Correct velocity field
        let mut u_corrected = u_star;
        self.pressure_solver.correct_velocity(
            &mut u_corrected,
            &p_correction,
            self.config.dt,
            rho,
            self.config.alpha_u,
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
        bc: &cfd_core::boundary::BoundaryCondition<T>,
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
    fn calculate_residual(
        &self,
        previous: &Vec<Vec<Vector2<T>>>,
        current: &Vec<Vec<Vector2<T>>>,
    ) -> T {
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        let mut max_diff = T::zero();
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let diff = (current[i][j] - previous[i][j]).norm();
                if diff > max_diff {
                    max_diff = diff;
                }
            }
        }

        max_diff
    }

    /// Get current velocity field
    pub fn velocity(&self) -> &Vec<Vec<Vector2<T>>> {
        &self.u
    }

    /// Get current pressure field
    pub fn pressure(&self) -> &Vec<Vec<T>> {
        &self.p
    }

    /// Get iteration count
    pub fn iterations(&self) -> usize {
        self.iterations
    }
}
