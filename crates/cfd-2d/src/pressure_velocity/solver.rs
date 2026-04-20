//! Main pressure-velocity coupling solver implementation
//!
//! # Theorem (SIMPLE Outer Loop Convergence)
//!
//! For a fixed linearization with positive pressure-correction coefficients,
//! the SIMPLE outer iteration is contractive only when the relaxation factors
//! keep the fixed-point map stable. Under those assumptions the continuity
//! residual decreases monotonically in practice; outside them convergence is
//! empirical rather than guaranteed.

use super::{PressureCorrectionSolver, PressureVelocityConfig, RhieChowInterpolation};
use crate::fields::SimulationFields;
use crate::grid::array2d::Array2D;
use crate::grid::StructuredGrid2D;
use crate::physics::{MomentumComponent, MomentumSolver};
use cfd_core::physics::boundary::BoundaryCondition;
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
    /// Predicted velocity workspace reused across SIMPLE iterations.
    u_star_workspace: Array2D<Vector2<T>>,
    /// Pressure correction workspace reused across SIMPLE iterations.
    p_correction_workspace: Array2D<T>,
    /// Simulation field workspace reused across momentum and pressure solves.
    state_workspace: SimulationFields<T>,
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
            u_star_workspace: Array2D::new(nx, ny, Vector2::zeros()),
            p_correction_workspace: Array2D::new(nx, ny, T::zero()),
            state_workspace: SimulationFields::new(nx, ny),
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
        bc: &BoundaryCondition<T>,
        nu: T,
        rho: T,
    ) -> cfd_core::error::Result<T> {
        let dt = self.config.dt;
        let alpha_u = self.config.alpha_u;
        let alpha_p = self.config.alpha_p;
        let dynamic_viscosity = rho * nu;

        self.apply_uniform_boundary_condition(bc)?;

        let u_source = self.u.as_slice();
        let p_source = self.p.as_slice();

        // Step 1: Solve momentum equations for predicted velocity.
        // Reuse one field workspace to avoid per-iteration heap allocation.
        {
            let momentum_solver = &mut self.momentum_solver;
            let pressure_solver = &mut self.pressure_solver;
            let state_buffer = &mut self.state_workspace;
            let u_star_workspace = &mut self.u_star_workspace;
            let p_correction_workspace = &mut self.p_correction_workspace;

            state_buffer.density.as_mut_slice().fill(rho);
            state_buffer.viscosity.as_mut_slice().fill(dynamic_viscosity);
            state_buffer.force_u.as_mut_slice().fill(T::zero());
            state_buffer.force_v.as_mut_slice().fill(T::zero());
            state_buffer.mask.as_mut_slice().fill(true);

            for ((dst_u, dst_v), vel) in state_buffer
                .u
                .as_mut_slice()
                .iter_mut()
                .zip(state_buffer.v.as_mut_slice().iter_mut())
                .zip(u_source.iter())
            {
                *dst_u = vel.x;
                *dst_v = vel.y;
            }
            state_buffer.p.as_mut_slice().clone_from_slice(p_source);

            momentum_solver.solve(MomentumComponent::U, state_buffer, dt)?;
            momentum_solver.solve(MomentumComponent::V, state_buffer, dt)?;

            state_buffer.copy_velocity_star_to(u_star_workspace);
            state_buffer.promote_velocity_star_to_current();

            // Step 2: Solve pressure correction equation.
            p_correction_workspace.fill(T::zero());
            pressure_solver.solve_pressure_correction(
                state_buffer,
                dt,
                rho,
                true,
                p_correction_workspace,
            )?;
        }

        // Step 3: Correct velocity field
        std::mem::swap(&mut self.u, &mut self.u_star_workspace);
        let (ap_u, _, ap_v, _) = self.momentum_solver.get_ap_coefficients();
        self.pressure_solver.correct_velocity(
            &mut self.u,
            &self.p_correction_workspace,
            ap_u,
            ap_v,
            rho,
            alpha_u,
            &self.state_workspace, // Pass buffer for mask access
        );

        // Step 4: Correct pressure field
        self.pressure_solver
            .correct_pressure(&mut self.p, &self.p_correction_workspace, alpha_p);

        // Step 5: Calculate residual for convergence check
        let residual = self.calculate_residual(&self.u_star_workspace, &self.u);
        self.iterations += 1;

        Ok(residual)
    }

    fn apply_uniform_boundary_condition(
        &mut self,
        bc: &BoundaryCondition<T>,
    ) -> cfd_core::error::Result<()> {
        if matches!(bc, BoundaryCondition::Periodic { .. }) {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "PressureVelocitySolver requires a uniform non-periodic boundary condition"
                    .to_string(),
            ));
        }

        let needs_update = {
            let boundaries = self.momentum_solver.boundary_conditions();
            ["north", "south", "east", "west"]
                .iter()
                .any(|name| boundaries.get(*name) != Some(bc))
        };

        if needs_update {
            self.momentum_solver
                .set_boundary("north".to_string(), bc.clone());
            self.momentum_solver
                .set_boundary("south".to_string(), bc.clone());
            self.momentum_solver
                .set_boundary("east".to_string(), bc.clone());
            self.momentum_solver
                .set_boundary("west".to_string(), bc.clone());
        }

        Ok(())
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
    fn calculate_residual(
        &self,
        previous: &Array2D<Vector2<T>>,
        current: &Array2D<Vector2<T>>,
    ) -> T {
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
    use cfd_core::physics::boundary::{BoundaryCondition, WallType};
    use nalgebra::Vector3;

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

    #[test]
    fn step_reuses_internal_workspaces_across_iterations() {
        let mut solver = make_solver(4, 4);
        let u_init = Array2D::new(4, 4, Vector2::new(1.0, 0.0));
        let p_init = Array2D::new(4, 4, 0.0);
        solver.set_initial_conditions(u_init, p_init);

        let bc = BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        };

        let first = solver
            .step(&bc, 0.01, 1.0)
            .expect("first pressure-velocity iteration");
        let second = solver
            .step(&bc, 0.01, 1.0)
            .expect("second pressure-velocity iteration");

        assert!(first.is_finite(), "first residual must be finite");
        assert!(second.is_finite(), "second residual must be finite");
        assert_eq!(solver.iterations(), 2);

        for velocity in solver.velocity().as_slice() {
            assert!(
                velocity.x.is_finite() && velocity.y.is_finite(),
                "velocity field must remain finite after repeated steps"
            );
        }
    }

    #[test]
    fn step_responds_to_boundary_and_viscosity_inputs() {
        let mut no_slip_solver = make_solver(6, 6);
        let mut moving_wall_solver = make_solver(6, 6);
        let mut low_viscosity_solver = make_solver(6, 6);
        let mut high_viscosity_solver = make_solver(6, 6);

        let u_init = Array2D::new(6, 6, Vector2::new(1.0, 0.0));
        let p_init = Array2D::new(6, 6, 0.0);

        no_slip_solver.set_initial_conditions(u_init.clone(), p_init.clone());
        moving_wall_solver.set_initial_conditions(u_init.clone(), p_init.clone());
        low_viscosity_solver.set_initial_conditions(u_init.clone(), p_init.clone());
        high_viscosity_solver.set_initial_conditions(u_init, p_init);

        let no_slip = BoundaryCondition::Wall {
            wall_type: WallType::NoSlip,
        };
        let moving = BoundaryCondition::Wall {
            wall_type: WallType::Moving {
                velocity: Vector3::new(1.0, 0.0, 0.0),
            },
        };

        let no_slip_residual = no_slip_solver
            .step(&no_slip, 0.01, 1.0)
            .expect("no-slip pressure-velocity iteration");
        let moving_residual = moving_wall_solver
            .step(&moving, 0.01, 1.0)
            .expect("moving-wall pressure-velocity iteration");
        let low_viscosity_residual = low_viscosity_solver
            .step(&no_slip, 0.002, 1.0)
            .expect("low-viscosity pressure-velocity iteration");
        let high_viscosity_residual = high_viscosity_solver
            .step(&no_slip, 0.05, 1.0)
            .expect("high-viscosity pressure-velocity iteration");

        assert!(no_slip_residual.is_finite());
        assert!(moving_residual.is_finite());
        assert!(low_viscosity_residual.is_finite());
        assert!(high_viscosity_residual.is_finite());

        let center = (2, 2);
        let no_slip_center = no_slip_solver.velocity()[(center.0, center.1)];
        let moving_center = moving_wall_solver.velocity()[(center.0, center.1)];
        let low_viscosity_center = low_viscosity_solver.velocity()[(center.0, center.1)];
        let high_viscosity_center = high_viscosity_solver.velocity()[(center.0, center.1)];

        assert!(
            (moving_center.x - no_slip_center.x).abs() > 1e-12
                || (moving_center.y - no_slip_center.y).abs() > 1e-12,
            "boundary condition must affect the corrected velocity field"
        );
        assert!(
            (high_viscosity_center.x - low_viscosity_center.x).abs() > 1e-12
                || (high_viscosity_center.y - low_viscosity_center.y).abs() > 1e-12,
            "viscosity input must affect the corrected velocity field"
        );
    }
}
