//! Main LBM solver implementation.
//!
//! This module provides the core solver that integrates collision,
//! streaming, boundary conditions, and macroscopic computations.

use crate::grid::{Grid2D, StructuredGrid2D};
use crate::solvers::lbm::{
    boundary::BoundaryHandler,
    collision::{BgkCollision, CollisionOperator},
    lattice::{equilibrium, D2Q9},
    macroscopic::MacroscopicQuantities,
    streaming::StreamingOperator,
};
use cfd_core::boundary::BoundaryCondition;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Configuration for LBM solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LbmConfig<T: RealField + Copy> {
    /// Relaxation time (related to viscosity)
    pub tau: T,
    /// Maximum number of time steps
    pub max_steps: usize,
    /// Convergence tolerance for steady-state
    pub tolerance: T,
    /// Output frequency (steps between outputs)
    pub output_frequency: usize,
    /// Enable verbose output
    pub verbose: bool,
}

impl<T: RealField + Copy + FromPrimitive> Default for LbmConfig<T> {
    fn default() -> Self {
        Self {
            tau: T::from_f64(1.0).unwrap_or_else(T::zero),
            max_steps: 10000,
            tolerance: T::from_f64(1e-6).unwrap_or_else(T::zero),
            output_frequency: 100,
            verbose: false,
        }
    }
}

/// Lattice Boltzmann Method solver for 2D incompressible flows
pub struct LbmSolver<T: RealField + Copy> {
    /// Solver configuration
    config: LbmConfig<T>,
    /// Distribution functions
    f: Vec<Vec<[T; 9]>>,
    /// Distribution functions buffer for streaming
    f_buffer: Vec<Vec<[T; 9]>>,
    /// Macroscopic quantities
    macroscopic: MacroscopicQuantities<T>,
    /// Collision operator
    collision: Box<dyn CollisionOperator<T> + Send + Sync>,
    /// Boundary handler
    boundary_handler: BoundaryHandler<T>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    /// Current step count
    step_count: usize,
    /// Buffer for convergence checking (zero-copy optimization)
    previous_velocity: Vec<Vec<[T; 2]>>,
}

impl<T: RealField + Copy + FromPrimitive> LbmSolver<T>
where
    T: Send + Sync + std::fmt::LowerExp,
{
    /// Create a new LBM solver
    pub fn new(config: LbmConfig<T>, grid: &StructuredGrid2D<T>) -> Self {
        let nx = grid.nx();
        let ny = grid.ny();
        let dx = grid.dx;
        let dy = grid.dy;

        // Initialize distribution functions
        let f = vec![vec![[T::zero(); 9]; nx]; ny];
        let f_buffer = vec![vec![[T::zero(); 9]; nx]; ny];

        // Initialize macroscopic quantities
        let macroscopic = MacroscopicQuantities::new(nx, ny);

        // Create collision operator
        let collision = Box::new(BgkCollision::new(config.tau));

        // Create boundary handler
        let boundary_handler = BoundaryHandler::new();

        // Preallocate convergence buffer (zero-copy optimization)
        let previous_velocity = vec![vec![[T::zero(), T::zero()]; nx]; ny];

        Self {
            config,
            f,
            f_buffer,
            macroscopic,
            collision,
            boundary_handler,
            nx,
            ny,
            dx,
            dy,
            step_count: 0,
            previous_velocity,
        }
    }

    /// Compute equilibrium distribution for given density and velocity
    pub fn equilibrium_distribution(&self, density: T, velocity: Vector2<T>) -> Vec<T> {
        let u = [velocity.x, velocity.y];
        let mut feq = vec![T::zero(); 9];

        for q in 0..9 {
            let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
            let lattice_vel = D2Q9::VELOCITIES[q];
            feq[q] = equilibrium(density, &u, q, weight, lattice_vel);
        }

        feq
    }

    /// Compute macroscopic density and velocity at a grid point
    pub fn compute_macroscopic(&self, i: usize, j: usize) -> (T, Vector2<T>) {
        let density = self.macroscopic.density[j][i];
        let velocity = Vector2::new(
            self.macroscopic.velocity[j][i][0],
            self.macroscopic.velocity[j][i][1],
        );
        (density, velocity)
    }

    /// Initialize the solver with functions for density and velocity
    pub fn initialize<F1, F2>(&mut self, density_fn: F1, velocity_fn: F2) -> Result<()>
    where
        F1: Fn(T, T) -> T,
        F2: Fn(T, T) -> Vector2<T>,
    {
        // Initialize using the provided functions
        for j in 0..self.ny {
            for i in 0..self.nx {
                let x = T::from_usize(i).unwrap_or_else(T::zero) * self.dx;
                let y = T::from_usize(j).unwrap_or_else(T::zero) * self.dy;

                let density = density_fn(x, y);
                let velocity = velocity_fn(x, y);
                let u_init = [velocity.x, velocity.y];

                // Initialize distribution functions to equilibrium
                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
                    let lattice_vel = D2Q9::VELOCITIES[q];
                    self.f[j][i][q] = equilibrium(density, &u_init, q, weight, lattice_vel);
                }

                // Set macroscopic quantities
                self.macroscopic.density[j][i] = density;
                self.macroscopic.velocity[j][i] = u_init;
            }
        }

        self.step_count = 0;
        Ok(())
    }

    /// Perform one time step
    pub fn step(
        &mut self,
        boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        // Update macroscopic quantities
        self.macroscopic.update_from_distributions(&self.f);

        // Collision step
        self.collision.collide(
            &mut self.f,
            &self.macroscopic.density,
            &self.macroscopic.velocity,
        );

        // Streaming step
        StreamingOperator::stream(&self.f, &mut self.f_buffer);
        std::mem::swap(&mut self.f, &mut self.f_buffer);

        // Apply boundary conditions
        self.boundary_handler.apply_boundaries(
            &mut self.f,
            &mut self.macroscopic.density,
            &mut self.macroscopic.velocity,
            boundaries,
        );

        self.step_count += 1;
        Ok(())
    }

    /// Run the solver until convergence or max steps
    pub fn solve(
        &mut self,
        boundaries: HashMap<(usize, usize), BoundaryCondition<T>>,
        initial_density: T,
        initial_velocity: Vector2<T>,
    ) -> Result<()> {
        // Initialize with constant functions
        let density_fn = |_x: T, _y: T| initial_density;
        let velocity_fn = |_x: T, _y: T| initial_velocity;
        self.initialize(density_fn, velocity_fn)?;

        let mut converged = false;
        // Copy initial velocity to previous_velocity buffer (zero-copy: single allocation)
        self.copy_velocity_to_buffer();

        for step in 0..self.config.max_steps {
            // Perform time step
            self.step(&boundaries)?;

            // Check convergence
            if step % self.config.output_frequency == 0 {
                let max_change = self.compute_max_velocity_change();

                if self.config.verbose {
                    println!("Step {step}: max velocity change = {max_change:e}");
                }

                if max_change < self.config.tolerance {
                    converged = true;
                    if self.config.verbose {
                        println!("Converged after {step} steps");
                    }
                    break;
                }

                // Update buffer for next comparison (zero-copy: reuse allocation)
                self.copy_velocity_to_buffer();
            }
        }

        if !converged && self.config.verbose {
            println!(
                "Warning: Did not converge after {} steps",
                self.config.max_steps
            );
        }

        Ok(())
    }

    /// Copy current velocity to buffer for convergence checking (zero-copy optimization)
    fn copy_velocity_to_buffer(&mut self) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                self.previous_velocity[j][i] = self.macroscopic.velocity[j][i];
            }
        }
    }

    /// Compute maximum velocity change for convergence check (zero-copy optimization)
    fn compute_max_velocity_change(&self) -> T {
        let mut max_change = T::zero();

        for j in 0..self.ny {
            for i in 0..self.nx {
                let du =
                    (self.macroscopic.velocity[j][i][0] - self.previous_velocity[j][i][0]).abs();
                let dv =
                    (self.macroscopic.velocity[j][i][1] - self.previous_velocity[j][i][1]).abs();
                let change = du.max(dv);
                if change > max_change {
                    max_change = change;
                }
            }
        }

        max_change
    }

    /// Get velocity field
    pub fn velocity_field(&self) -> &Vec<Vec<[T; 2]>> {
        &self.macroscopic.velocity
    }

    /// Get density field
    pub fn density_field(&self) -> &Vec<Vec<T>> {
        &self.macroscopic.density
    }

    /// Get velocity at a specific point
    pub fn velocity_at(&self, i: usize, j: usize) -> Option<&[T; 2]> {
        if i < self.nx && j < self.ny {
            Some(&self.macroscopic.velocity[j][i])
        } else {
            None
        }
    }

    /// Get density at a specific point
    pub fn density_at(&self, i: usize, j: usize) -> Option<&T> {
        if i < self.nx && j < self.ny {
            Some(&self.macroscopic.density[j][i])
        } else {
            None
        }
    }

    /// Get macroscopic quantities
    pub fn get_macroscopic(&self) -> (&Vec<Vec<[T; 2]>>, &Vec<Vec<T>>) {
        (&self.macroscopic.velocity, &self.macroscopic.density)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::StructuredGrid2D;
    use approx::assert_relative_eq;

    #[test]
    fn test_equilibrium_distribution() -> Result<()> {
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let solver = LbmSolver::new(config, &grid);

        let rho = 1.0;
        let u = Vector2::new(0.1, 0.0);

        let feq = solver.equilibrium_distribution(rho, u);

        // Check that sum of distributions equals density
        let sum: f64 = feq.iter().sum();
        assert_relative_eq!(sum, rho, epsilon = 1e-10);

        Ok(())
    }

    #[test]
    fn test_initialization() -> Result<()> {
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let mut solver = LbmSolver::new(config, &grid);

        let initial_density = |_x: f64, _y: f64| 1.0;
        let initial_velocity = |_x: f64, _y: f64| Vector2::new(0.0, 0.0);

        solver.initialize(initial_density, initial_velocity)?;

        // Check that macroscopic properties match initial conditions
        let (rho, u) = solver.compute_macroscopic(5, 5);
        assert_relative_eq!(rho, 1.0, epsilon = 1e-10);
        assert_relative_eq!(u.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(u.y, 0.0, epsilon = 1e-10);

        Ok(())
    }
}
