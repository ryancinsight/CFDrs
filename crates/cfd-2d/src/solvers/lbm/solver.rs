//! Main LBM solver implementation.
//!
//! This module provides the core solver that integrates collision,
//! streaming, boundary conditions, and macroscopic computations.

use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use cfd_core::{Result, Error, BoundaryCondition};
use crate::grid::{Grid2D, StructuredGrid2D};
use crate::solvers::lbm::{
    lattice::{D2Q9, equilibrium},
    collision::{CollisionOperator, BgkCollision},
    streaming::StreamingOperator,
    boundary::BoundaryHandler,
    macroscopic::{MacroscopicQuantities, compute_density, compute_velocity},
};

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
    /// Temporary distribution functions for streaming
    f_temp: Vec<Vec<[T; 9]>>,
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
        let f_temp = vec![vec![[T::zero(); 9]; nx]; ny];
        
        // Initialize macroscopic quantities
        let macroscopic = MacroscopicQuantities::new(nx, ny);
        
        // Create collision operator
        let collision = Box::new(BgkCollision::new(config.tau));
        
        // Create boundary handler
        let boundary_handler = BoundaryHandler::new();
        
        Self {
            config,
            f,
            f_temp,
            macroscopic,
            collision,
            boundary_handler,
            nx,
            ny,
            dx,
            dy,
            step_count: 0,
        }
    }
    
    /// Initialize the solver with uniform flow
    pub fn initialize(&mut self, initial_density: T, initial_velocity: Vector2<T>) -> Result<()> {
        let u_init = [initial_velocity.x, initial_velocity.y];
        
        // Initialize distribution functions to equilibrium
        for j in 0..self.ny {
            for i in 0..self.nx {
                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
                    let lattice_vel = &D2Q9::VELOCITIES[q];
                    self.f[j][i][q] = equilibrium(
                        initial_density,
                        &u_init,
                        q,
                        weight,
                        lattice_vel,
                    );
                }
                
                // Set macroscopic quantities
                self.macroscopic.density[j][i] = initial_density;
                self.macroscopic.velocity[j][i] = u_init;
            }
        }
        
        self.step_count = 0;
        Ok(())
    }
    
    /// Perform one time step
    pub fn step(&mut self, boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>) -> Result<()> {
        // Update macroscopic quantities
        self.macroscopic.update_from_distributions(&self.f);
        
        // Collision step
        self.collision.collide(
            &mut self.f,
            &self.macroscopic.density,
            &self.macroscopic.velocity,
        );
        
        // Streaming step
        StreamingOperator::stream(&self.f, &mut self.f_temp);
        std::mem::swap(&mut self.f, &mut self.f_temp);
        
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
        // Initialize
        self.initialize(initial_density, initial_velocity)?;
        
        let mut converged = false;
        let mut previous_velocity = self.macroscopic.velocity.clone();
        
        for step in 0..self.config.max_steps {
            // Perform time step
            self.step(&boundaries)?;
            
            // Check convergence
            if step % self.config.output_frequency == 0 {
                let max_change = self.compute_max_velocity_change(&previous_velocity);
                
                if self.config.verbose {
                    println!("Step {}: max velocity change = {:e}", step, max_change);
                }
                
                if max_change < self.config.tolerance {
                    converged = true;
                    if self.config.verbose {
                        println!("Converged after {} steps", step);
                    }
                    break;
                }
                
                previous_velocity = self.macroscopic.velocity.clone();
            }
        }
        
        if !converged && self.config.verbose {
            println!("Warning: Did not converge after {} steps", self.config.max_steps);
        }
        
        Ok(())
    }
    
    /// Compute maximum velocity change for convergence check
    fn compute_max_velocity_change(&self, previous_velocity: &Vec<Vec<[T; 2]>>) -> T {
        let mut max_change = T::zero();
        
        for j in 0..self.ny {
            for i in 0..self.nx {
                let du = (self.macroscopic.velocity[j][i][0] - previous_velocity[j][i][0]).abs();
                let dv = (self.macroscopic.velocity[j][i][1] - previous_velocity[j][i][1]).abs();
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
    
    #[test]
    fn test_solver_creation() {
        let config = LbmConfig::<f64>::default();
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0).unwrap();
        let solver = LbmSolver::new(config, &grid);
        
        assert_eq!(solver.nx, 10);
        assert_eq!(solver.ny, 10);
        assert_eq!(solver.step_count, 0);
    }
    
    #[test]
    fn test_solver_initialization() {
        let config = LbmConfig::<f64>::default();
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0).unwrap();
        let mut solver = LbmSolver::new(config, &grid);
        
        let initial_density = 1.0;
        let initial_velocity = Vector2::new(0.1, 0.0);
        
        solver.initialize(initial_density, initial_velocity).unwrap();
        
        // Check that density and velocity are set correctly
        for j in 0..10 {
            for i in 0..10 {
                assert!((solver.macroscopic.density[j][i] - initial_density).abs() < 1e-10);
                assert!((solver.macroscopic.velocity[j][i][0] - 0.1).abs() < 1e-10);
                assert!(solver.macroscopic.velocity[j][i][1].abs() < 1e-10);
            }
        }
    }
}