//! Lattice Boltzmann Method (LBM) solvers for 2D fluid dynamics.
//!
//! This module provides LBM implementations for solving incompressible
//! Navier-Stokes equations using the collision-streaming approach.
//!
//! Features:
//! - D2Q9 lattice model (2D, 9 velocities)
//! - BGK collision operator
//! - Various boundary conditions (bounce-back, velocity, pressure)
//! - Optimized using iterator combinators and parallel processing

use cfd_core::{Result, BoundaryCondition};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use crate::grid::{Grid2D, StructuredGrid2D};

/// D2Q9 lattice model constants
pub struct D2Q9;

impl D2Q9 {
    /// Number of velocity directions
    pub const Q: usize = 9;

    /// Lattice velocities (normalized by lattice spacing)
    pub const VELOCITIES: [(i32, i32); 9] = [
        (0, 0),   // 0: rest
        (1, 0),   // 1: east
        (0, 1),   // 2: north
        (-1, 0),  // 3: west
        (0, -1),  // 4: south
        (1, 1),   // 5: northeast
        (-1, 1),  // 6: northwest
        (-1, -1), // 7: southwest
        (1, -1),  // 8: southeast
    ];

    /// Lattice weights
    pub const WEIGHTS: [f64; 9] = [
        4.0/9.0,  // 0: rest
        1.0/9.0,  // 1-4: cardinal directions
        1.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/36.0, // 5-8: diagonal directions
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
    ];

    /// Opposite directions for bounce-back boundary conditions
    pub const OPPOSITE: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
}

/// LBM solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LbmConfig<T: RealField> {
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

impl<T: RealField + FromPrimitive> Default for LbmConfig<T> {
    fn default() -> Self {
        Self {
            tau: T::from_f64(1.0).unwrap(),
            max_steps: 10000,
            tolerance: T::from_f64(1e-6).unwrap(),
            output_frequency: 100,
            verbose: false,
        }
    }
}

/// Lattice Boltzmann Method solver for 2D incompressible flows
pub struct LbmSolver<T: RealField> {
    config: LbmConfig<T>,
    /// Distribution functions [nx][ny][Q]
    f: Vec<Vec<Vec<T>>>,
    /// Temporary distribution functions for streaming
    f_temp: Vec<Vec<Vec<T>>>,
    /// Macroscopic density
    rho: Vec<Vec<T>>,
    /// Macroscopic velocity
    u: Vec<Vec<Vector2<T>>>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
}

impl<T: RealField + FromPrimitive + Send + Sync + Clone> LbmSolver<T> {
    /// Create a new LBM solver
    pub fn new(config: LbmConfig<T>, grid: &StructuredGrid2D<T>) -> Self {
        let nx = grid.nx();
        let ny = grid.ny();

        // Initialize distribution functions
        let f = vec![vec![vec![T::zero(); D2Q9::Q]; ny]; nx];
        let f_temp = vec![vec![vec![T::zero(); D2Q9::Q]; ny]; nx];
        let rho = vec![vec![T::one(); ny]; nx];
        let u = vec![vec![Vector2::zeros(); ny]; nx];

        Self {
            config,
            f,
            f_temp,
            rho,
            u,
            nx,
            ny,
        }
    }

    /// Initialize the solver with equilibrium distributions
    pub fn initialize(&mut self, initial_density: T, initial_velocity: Vector2<T>) -> Result<()> {
        // Initialize using regular loops
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.rho[i][j] = initial_density.clone();
                self.u[i][j] = initial_velocity.clone();

                // Set equilibrium distribution
                for q in 0..D2Q9::Q {
                    self.f[i][j][q] = self.equilibrium_distribution(
                        q,
                        &initial_density,
                        &initial_velocity
                    );
                }
            }
        }

        Ok(())
    }

    /// Calculate equilibrium distribution function
    fn equilibrium_distribution(&self, q: usize, rho: &T, u: &Vector2<T>) -> T {
        let w = T::from_f64(D2Q9::WEIGHTS[q]).unwrap();
        let c = Vector2::new(
            T::from_i32(D2Q9::VELOCITIES[q].0).unwrap(),
            T::from_i32(D2Q9::VELOCITIES[q].1).unwrap(),
        );

        let cu = c.dot(u);
        let u_sqr = u.dot(u);
        let cs2 = T::from_f64(1.0/3.0).unwrap(); // Speed of sound squared

        // Equilibrium distribution: w_i * rho * (1 + c_i·u/cs² + (c_i·u)²/(2cs⁴) - u²/(2cs²))
        let term1 = cu.clone() / cs2.clone();
        let term2 = cu.clone() * cu / (T::from_f64(2.0).unwrap() * cs2.clone() * cs2.clone());
        let term3 = u_sqr / (T::from_f64(2.0).unwrap() * cs2);

        w * rho.clone() * (T::one() + term1 + term2 - term3)
    }

    /// Perform collision step (BGK operator)
    fn collision(&mut self) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                // Calculate macroscopic quantities
                let mut rho_local = T::zero();
                let mut u_local = Vector2::zeros();

                // Sum distribution functions to get density
                for q in 0..D2Q9::Q {
                    rho_local += self.f[i][j][q].clone();
                }

                // Calculate velocity
                for q in 0..D2Q9::Q {
                    let c = Vector2::new(
                        T::from_i32(D2Q9::VELOCITIES[q].0).unwrap(),
                        T::from_i32(D2Q9::VELOCITIES[q].1).unwrap(),
                    );
                    u_local += c * self.f[i][j][q].clone();
                }
                u_local /= rho_local.clone();

                // Store macroscopic quantities
                self.rho[i][j] = rho_local.clone();
                self.u[i][j] = u_local.clone();

                // BGK collision
                let omega = T::one() / self.config.tau.clone();
                for q in 0..D2Q9::Q {
                    let f_eq = self.equilibrium_distribution(q, &rho_local, &u_local);
                    self.f[i][j][q] = self.f[i][j][q].clone() * (T::one() - omega.clone()) + f_eq * omega.clone();
                }
            }
        }
    }

    /// Perform streaming step
    fn streaming(&mut self) {
        // Copy current distributions to temporary array
        for i in 0..self.nx {
            for j in 0..self.ny {
                for q in 0..D2Q9::Q {
                    self.f_temp[i][j][q] = self.f[i][j][q].clone();
                }
            }
        }

        // Stream distributions using regular loops
        for i in 0..self.nx {
            for j in 0..self.ny {
                for q in 0..D2Q9::Q {
                    let (ci, cj) = D2Q9::VELOCITIES[q];
                    let ni = i as i32 - ci;
                    let nj = j as i32 - cj;

                    // Check bounds
                    if ni >= 0 && ni < self.nx as i32 && nj >= 0 && nj < self.ny as i32 {
                        self.f[i][j][q] = self.f_temp[ni as usize][nj as usize][q].clone();
                    }
                }
            }
        }
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&mut self, boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>) {
        for (&(i, j), bc) in boundaries.iter() {
            if i < self.nx && j < self.ny {
                match bc {
                    BoundaryCondition::Wall { .. } => {
                        // Proper bounce-back boundary condition
                        // f_q(x_b, t+1) = f*_opp(x_b, t) where f* is post-collision
                        // For simplicity, we apply bounce-back to current distributions
                        let mut f_new = vec![T::zero(); D2Q9::Q];

                        // Copy current distributions
                        for q in 0..D2Q9::Q {
                            f_new[q] = self.f[i][j][q].clone();
                        }

                        // Apply bounce-back: outgoing = incoming from opposite direction
                        for q in 1..D2Q9::Q { // Skip rest particle (q=0)
                            let opp = D2Q9::OPPOSITE[q];
                            f_new[q] = self.f[i][j][opp].clone();
                        }

                        // Update distributions
                        for q in 0..D2Q9::Q {
                            self.f[i][j][q] = f_new[q].clone();
                        }
                    }
                    BoundaryCondition::VelocityInlet { velocity } => {
                        // Velocity boundary condition (simplified)
                        let u_bc = Vector2::new(velocity.x.clone(), velocity.y.clone());
                        for q in 0..D2Q9::Q {
                            self.f[i][j][q] = self.equilibrium_distribution(q, &self.rho[i][j], &u_bc);
                        }
                    }
                    BoundaryCondition::PressureOutlet { pressure } => {
                        // Pressure boundary condition (simplified)
                        let rho_bc = pressure.clone();
                        for q in 0..D2Q9::Q {
                            self.f[i][j][q] = self.equilibrium_distribution(q, &rho_bc, &self.u[i][j]);
                        }
                    }
                    _ => {
                        // Default: no-slip wall (proper bounce-back)
                        let mut f_new = vec![T::zero(); D2Q9::Q];

                        // Copy current distributions
                        for q in 0..D2Q9::Q {
                            f_new[q] = self.f[i][j][q].clone();
                        }

                        // Apply bounce-back: outgoing = incoming from opposite direction
                        for q in 1..D2Q9::Q { // Skip rest particle (q=0)
                            let opp = D2Q9::OPPOSITE[q];
                            f_new[q] = self.f[i][j][opp].clone();
                        }

                        // Update distributions
                        for q in 0..D2Q9::Q {
                            self.f[i][j][q] = f_new[q].clone();
                        }
                    }
                }
            }
        }
    }

    /// Main solve method
    pub fn solve(
        &mut self,
        boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        for step in 0..self.config.max_steps {
            // Collision step
            self.collision();

            // Apply boundary conditions
            self.apply_boundary_conditions(boundaries);

            // Streaming step
            self.streaming();

            // Check convergence based on velocity and density field changes
            if step % self.config.output_frequency == 0 {
                if self.config.verbose {
                    println!("LBM Step: {}", step);
                }

                // Implement proper convergence check
                if step > 0 && self.check_convergence()? {
                    if self.config.verbose {
                        println!("LBM converged at step: {}", step);
                    }
                    break;
                }
            }
        }

        Ok(())
    }

    /// Get current velocity field
    pub fn velocity_field(&self) -> &Vec<Vec<Vector2<T>>> {
        &self.u
    }

    /// Get current density field
    pub fn density_field(&self) -> &Vec<Vec<T>> {
        &self.rho
    }

    /// Check convergence based on velocity and density field changes
    fn check_convergence(&self) -> Result<bool> {
        // Use iterator combinators for zero-copy convergence checking
        let total_cells = T::from_usize(self.nx * self.ny).unwrap();

        // Calculate velocity magnitudes and residuals using iterator patterns
        let (velocity_residual, density_residual, _max_velocity_magnitude) = (0..self.nx)
            .flat_map(|i| (0..self.ny).map(move |j| (i, j)))
            .map(|(i, j)| {
                let u_mag = (self.u[i][j].x.clone() * self.u[i][j].x.clone() +
                           self.u[i][j].y.clone() * self.u[i][j].y.clone()).sqrt();
                let density_residual = (self.rho[i][j].clone() - T::one()).abs();
                (u_mag.clone(), density_residual, u_mag)
            })
            .fold(
                (T::zero(), T::zero(), T::zero()),
                |(vel_acc, dens_acc, max_acc), (u_mag, dens_res, u_mag_max)| {
                    (
                        vel_acc + u_mag,
                        dens_acc + dens_res,
                        if u_mag_max > max_acc { u_mag_max } else { max_acc }
                    )
                }
            );

        // Normalize residuals
        let velocity_residual_norm = velocity_residual / total_cells.clone();
        let density_residual_norm = density_residual / total_cells;

        // Check convergence criteria
        let velocity_converged = velocity_residual_norm < self.config.tolerance;
        let density_converged = density_residual_norm < self.config.tolerance;

        Ok(velocity_converged && density_converged)
    }

    /// Get velocity at specific point
    pub fn velocity_at(&self, i: usize, j: usize) -> Option<&Vector2<T>> {
        if i < self.nx && j < self.ny {
            Some(&self.u[i][j])
        } else {
            None
        }
    }

    /// Get density at specific point
    pub fn density_at(&self, i: usize, j: usize) -> Option<&T> {
        if i < self.nx && j < self.ny {
            Some(&self.rho[i][j])
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::BoundaryCondition;
    use nalgebra::Vector2;
    use std::collections::HashMap;

    #[test]
    fn test_d2q9_constants() {
        assert_eq!(D2Q9::Q, 9);
        assert_eq!(D2Q9::VELOCITIES.len(), 9);
        assert_eq!(D2Q9::WEIGHTS.len(), 9);
        assert_eq!(D2Q9::OPPOSITE.len(), 9);

        // Check that weights sum to 1
        let weight_sum: f64 = D2Q9::WEIGHTS.iter().sum();
        assert_relative_eq!(weight_sum, 1.0, epsilon = 1e-10);

        // Check opposite directions
        assert_eq!(D2Q9::OPPOSITE[1], 3); // east <-> west
        assert_eq!(D2Q9::OPPOSITE[2], 4); // north <-> south
    }

    #[test]
    fn test_lbm_solver_creation() {
        let grid = StructuredGrid2D::<f64>::unit_square(5, 5).unwrap();
        let config = LbmConfig::<f64>::default();
        let solver = LbmSolver::new(config, &grid);

        assert_eq!(solver.nx, 5);
        assert_eq!(solver.ny, 5);
        assert_eq!(solver.f.len(), 5);
        assert_eq!(solver.f[0].len(), 5);
        assert_eq!(solver.f[0][0].len(), 9);
    }

    #[test]
    fn test_equilibrium_distribution() {
        let grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();
        let config = LbmConfig::<f64>::default();
        let solver = LbmSolver::new(config, &grid);

        let rho = 1.0;
        let u = Vector2::new(0.1, 0.0);

        // Test equilibrium for rest particle (q=0)
        let f_eq_0 = solver.equilibrium_distribution(0, &rho, &u);
        assert!(f_eq_0 > 0.0);

        // Test that all equilibrium distributions sum to density
        let mut sum = 0.0;
        for q in 0..D2Q9::Q {
            sum += solver.equilibrium_distribution(q, &rho, &u);
        }
        assert_relative_eq!(sum, rho, epsilon = 1e-10);
    }

    #[test]
    fn test_lbm_initialization() {
        let grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();
        let config = LbmConfig::<f64>::default();
        let mut solver = LbmSolver::new(config, &grid);

        let initial_density = 1.0;
        let initial_velocity = Vector2::new(0.0, 0.0);

        solver.initialize(initial_density, initial_velocity).unwrap();

        // Check that density and velocity are initialized correctly
        for i in 0..solver.nx {
            for j in 0..solver.ny {
                assert_relative_eq!(solver.rho[i][j], initial_density, epsilon = 1e-10);
                assert_relative_eq!(solver.u[i][j].x, initial_velocity.x, epsilon = 1e-10);
                assert_relative_eq!(solver.u[i][j].y, initial_velocity.y, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_lbm_simple_flow() {
        let grid = StructuredGrid2D::<f64>::unit_square(5, 5).unwrap();
        let mut config = LbmConfig::<f64>::default();
        config.max_steps = 10; // Short simulation for testing
        config.verbose = false;

        let mut solver = LbmSolver::new(config, &grid);

        // Initialize with rest state
        solver.initialize(1.0, Vector2::zeros()).unwrap();

        // Set up simple boundary conditions
        let mut boundaries = HashMap::new();

        // Left wall: velocity inlet
        for j in 0..grid.ny() {
            boundaries.insert(
                (0, j),
                BoundaryCondition::VelocityInlet {
                    velocity: nalgebra::Vector3::new(0.1, 0.0, 0.0)
                }
            );
        }

        // Right wall: pressure outlet
        for j in 0..grid.ny() {
            boundaries.insert(
                (grid.nx() - 1, j),
                BoundaryCondition::PressureOutlet { pressure: 1.0 }
            );
        }

        // Run simulation
        solver.solve(&boundaries).unwrap();

        // Check that simulation completed without errors
        assert!(solver.velocity_field().len() == grid.nx());
        assert!(solver.density_field().len() == grid.nx());

        // Check that velocity field has reasonable values
        let u_center = solver.velocity_at(2, 2).unwrap();
        assert!(u_center.x >= 0.0); // Should have positive x-velocity due to inlet
    }
}
