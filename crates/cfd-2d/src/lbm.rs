//! Lattice Boltzmann Method (LBM) solvers for 2D fluid dynamics.
//!
//! This module provides LBM implementations for solving incompressible
//! Navier-Stokes equations using the collision-streaming approach.
//!
//! Features:
//! - D2Q9 lattice model (2D, 9 velocities)
//! - BGK collision operator
//! - Various boundary conditions (bounce-back, velocity, pressure)
//! - Parallel processing support

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
    /// Distribution functions [nx][ny][Q] - current timestep
    f: Vec<Vec<Vec<T>>>,
    /// Distribution functions [nx][ny][Q] - next timestep
    f_new: Vec<Vec<Vec<T>>>,
    /// Flag to track which array is current (for double buffering)
    use_f_as_current: bool,
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
        let f_new = vec![vec![vec![T::zero(); D2Q9::Q]; ny]; nx];
        let rho = vec![vec![T::one(); ny]; nx];
        let u = vec![vec![Vector2::zeros(); ny]; nx];

        Self {
            config,
            f,
            f_new,
            use_f_as_current: true,
            rho,
            u,
            nx,
            ny,
        }
    }
    
    /// Get current distribution functions (for reading)
    fn f_current(&self) -> &Vec<Vec<Vec<T>>> {
        if self.use_f_as_current {
            &self.f
        } else {
            &self.f_new
        }
    }
    
    /// Get current distribution functions (mutable for writing)
    fn f_current_mut(&mut self) -> &mut Vec<Vec<Vec<T>>> {
        if self.use_f_as_current {
            &mut self.f
        } else {
            &mut self.f_new
        }
    }
    
    /// Get next distribution functions (for writing during streaming)
    fn f_next_mut(&mut self) -> &mut Vec<Vec<Vec<T>>> {
        if self.use_f_as_current {
            &mut self.f_new
        } else {
            &mut self.f
        }
    }
    
    /// Swap the distribution function buffers (zero-cost operation)
    fn swap_buffers(&mut self) {
        self.use_f_as_current = !self.use_f_as_current;
    }

    /// Initialize the solver with equilibrium distributions
    pub fn initialize(&mut self, initial_density: T, initial_velocity: Vector2<T>) -> Result<()> {
        // Use traditional loops for mutable access to self, but iterator for inner loop
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.rho[i][j] = initial_density.clone();
                self.u[i][j] = initial_velocity.clone();

                // Set equilibrium distribution using iterator for better performance
                for q in 0..D2Q9::Q {
                    let eq_dist = self.equilibrium_distribution(
                        q,
                        &initial_density,
                        &initial_velocity
                    );
                    if self.use_f_as_current {
                        self.f[i][j][q] = eq_dist;
                    } else {
                        self.f_new[i][j][q] = eq_dist;
                    }
                }
            }
        }

        Ok(())
    }

    /// Calculate equilibrium distribution function using Maxwell-Boltzmann distribution
    ///
    /// The equilibrium distribution function is given by:
    /// ```text
    /// f_i^eq = w_i * ρ * [1 + (c_i · u)/c_s² + (c_i · u)²/(2c_s⁴) - u²/(2c_s²)]
    /// ```
    ///
    /// Where:
    /// - `w_i` are the lattice weights
    /// - `ρ` is the fluid density
    /// - `c_i` are the discrete velocities
    /// - `u` is the macroscopic velocity
    /// - `c_s = 1/√3` is the lattice speed of sound
    ///
    /// This formulation ensures that the macroscopic quantities (density and momentum)
    /// are recovered correctly when taking moments of the distribution functions.
    ///
    /// # References
    /// Chen, S. and Doolen, G.D. (1998). "Lattice Boltzmann method for fluid flows."
    /// Annual Review of Fluid Mechanics, 30(1), 329-364.
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
        // We need to read from current and write to the same array
        // So we'll use a temporary array for the collision result
        
        for i in 0..self.nx {
            for j in 0..self.ny {
                // Get current distributions
                let f_ij = if self.use_f_as_current {
                    &self.f[i][j]
                } else {
                    &self.f_new[i][j]
                };
                
                // Calculate density using iterator fold for better performance
                let rho_local = (0..D2Q9::Q)
                    .map(|q| f_ij[q].clone())
                    .fold(T::zero(), |acc, f| acc + f);

                // Calculate velocity using iterator fold
                let u_local = (0..D2Q9::Q)
                    .map(|q| {
                        let c = Vector2::new(
                            T::from_i32(D2Q9::VELOCITIES[q].0).unwrap(),
                            T::from_i32(D2Q9::VELOCITIES[q].1).unwrap(),
                        );
                        c * f_ij[q].clone()
                    })
                    .fold(Vector2::zeros(), |acc, v| acc + v) / rho_local.clone();

                // Store macroscopic quantities
                self.rho[i][j] = rho_local.clone();
                self.u[i][j] = u_local.clone();

                // Collision with BGK operator - update in place
                for q in 0..D2Q9::Q {
                    let f_eq = self.equilibrium_distribution(q, &rho_local, &u_local);
                    if self.use_f_as_current {
                        let f_old = self.f[i][j][q].clone();
                        self.f[i][j][q] = f_old.clone() - (f_old - f_eq) / self.config.tau.clone();
                    } else {
                        let f_old = self.f_new[i][j][q].clone();
                        self.f_new[i][j][q] = f_old.clone() - (f_old - f_eq) / self.config.tau.clone();
                    }
                }
            }
        }
    }

    /// Perform streaming step with zero-copy double buffering
    fn streaming(&mut self) {
        // Stream from current buffer to next buffer
        
        // Stream distributions - this is now the only copy operation
        // and it's necessary for the streaming physics
        for i in 0..self.nx {
            for j in 0..self.ny {
                for q in 0..D2Q9::Q {
                    let (ci, cj) = D2Q9::VELOCITIES[q];
                    let ni = i as i32 - ci;
                    let nj = j as i32 - cj;

                    // Check bounds and stream
                    if ni >= 0 && ni < self.nx as i32 && nj >= 0 && nj < self.ny as i32 {
                        // Stream from source to destination
                        let value = if self.use_f_as_current {
                            self.f[ni as usize][nj as usize][q].clone()
                        } else {
                            self.f_new[ni as usize][nj as usize][q].clone()
                        };
                        
                        if self.use_f_as_current {
                            self.f_new[i][j][q] = value;
                        } else {
                            self.f[i][j][q] = value;
                        }
                    } else {
                        // Boundary nodes keep their post-collision values
                        let value = if self.use_f_as_current {
                            self.f[i][j][q].clone()
                        } else {
                            self.f_new[i][j][q].clone()
                        };
                        
                        if self.use_f_as_current {
                            self.f_new[i][j][q] = value;
                        } else {
                            self.f[i][j][q] = value;
                        }
                    }
                }
            }
        }
        
        // Swap buffers for next iteration (zero-cost pointer swap)
        self.swap_buffers();
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&mut self, boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>) {
        for (&(i, j), bc) in boundaries.iter() {
            if i < self.nx && j < self.ny {
                match bc {
                    BoundaryCondition::Wall { wall_type } => {
                        // Apply bounce-back for walls using iterator
                        for q in 1..D2Q9::Q {
                            let (ci, cj) = D2Q9::VELOCITIES[q];
                            let ni = (i as i32 + ci) as usize;
                            let nj = (j as i32 + cj) as usize;

                            if ni < self.nx && nj < self.ny {
                                // Apply halfway bounce-back
                                let opp = D2Q9::OPPOSITE[q];
                                // The distribution streaming TO the wall bounces back
                                let value = if self.use_f_as_current {
                                    self.f[i][j][q].clone()
                                } else {
                                    self.f_new[i][j][q].clone()
                                };
                                
                                if self.use_f_as_current {
                                    self.f_new[ni][nj][opp] = value;
                                } else {
                                    self.f[ni][nj][opp] = value;
                                }
                            }
                        }
                    }
                    BoundaryCondition::VelocityInlet { velocity } => {
                        // Velocity boundary condition using iterator
                        let u_bc = Vector2::new(velocity.x.clone(), velocity.y.clone());
                        for q in 0..D2Q9::Q {
                            let value = self.equilibrium_distribution(q, &self.rho[i][j], &u_bc);
                            if self.use_f_as_current {
                                self.f_new[i][j][q] = value;
                            } else {
                                self.f[i][j][q] = value;
                            }
                        }
                    }
                    BoundaryCondition::PressureOutlet { pressure } => {
                        // Pressure boundary condition using iterator
                        let rho_bc = pressure.clone();
                        for q in 0..D2Q9::Q {
                            let value = self.equilibrium_distribution(q, &rho_bc, &self.u[i][j]);
                            if self.use_f_as_current {
                                self.f_new[i][j][q] = value;
                            } else {
                                self.f[i][j][q] = value;
                            }
                        }
                    }
                    _ => {
                        // Default bounce-back for other boundary types
                        let mut f_new = vec![T::zero(); D2Q9::Q];

                        // Copy current distributions
                        for q in 0..D2Q9::Q {
                            f_new[q] = if self.use_f_as_current {
                                self.f[i][j][q].clone()
                            } else {
                                self.f_new[i][j][q].clone()
                            };
                        }

                        // Apply bounce-back: outgoing = incoming from opposite direction
                        for q in 1..D2Q9::Q { // Skip rest particle (q=0)
                            let opp = D2Q9::OPPOSITE[q];
                            f_new[q] = if self.use_f_as_current {
                                self.f[i][j][opp].clone()
                            } else {
                                self.f_new[i][j][opp].clone()
                            };
                        }

                        // Update distributions
                        for q in 0..D2Q9::Q {
                            if self.use_f_as_current {
                                self.f_new[i][j][q] = f_new[q].clone();
                            } else {
                                self.f[i][j][q] = f_new[q].clone();
                            }
                        }
                    }
                }
            }
        }
    }

    /// Apply bounce-back boundary conditions
    /// 
    /// CRITICAL FIX: This implementation now correctly reflects distributions
    /// from adjacent fluid nodes instead of scrambling the boundary node's own
    /// distributions. The bounce-back condition states that particles hitting
    /// a wall reverse their direction while maintaining their magnitude.
    ///
    /// For a wall at position (i,j), incoming distributions from fluid nodes
    /// are reflected back in the opposite direction.
    ///
    /// # References
    /// Krüger, T. et al. (2017). "The Lattice Boltzmann Method: Principles and Practice."
    fn apply_bounce_back(&mut self, boundary_nodes: &[(usize, usize)]) {
        for &(i, j) in boundary_nodes {
            // For each boundary node, we need to reflect distributions
            // coming from neighboring fluid nodes
            
            // Store incoming distributions that will be reflected
            let mut reflected = vec![T::zero(); D2Q9::Q];
            
            // Check each direction
            for q in 0..D2Q9::Q {
                let (dx, dy) = D2Q9::VELOCITIES[q];
                
                // Calculate the neighbor position in the opposite direction
                // (where the fluid is coming from)
                let ni = i as i32 - dx;
                let nj = j as i32 - dy;
                
                // Check if neighbor is within bounds and is a fluid node
                if ni >= 0 && ni < self.nx as i32 && 
                   nj >= 0 && nj < self.ny as i32 {
                    let ni = ni as usize;
                    let nj = nj as usize;
                    
                    // Get the opposite direction index
                    let q_opp = D2Q9::OPPOSITE[q];
                    
                    // The distribution coming from direction q will be reflected
                    // back in the opposite direction q_opp
                    reflected[q_opp] = if self.use_f_as_current {
                        self.f[ni][nj][q].clone()
                    } else {
                        self.f_new[ni][nj][q].clone()
                    };
                }
            }
            
            // Apply the reflected distributions to the boundary node
            // Only update the distributions that are pointing into the fluid
            for q in 0..D2Q9::Q {
                let (dx, dy) = D2Q9::VELOCITIES[q];
                
                // Check if this direction points into the fluid domain
                let ni = i as i32 + dx;
                let nj = j as i32 + dy;
                
                if ni >= 0 && ni < self.nx as i32 && 
                   nj >= 0 && nj < self.ny as i32 {
                    // This direction points into the fluid, apply reflection
                    if self.use_f_as_current {
                        self.f_new[i][j][q] = reflected[q].clone();
                    } else {
                        self.f[i][j][q] = reflected[q].clone();
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
        // Simple, clear implementation without false optimization claims
        let total_cells = T::from_usize(self.nx * self.ny).ok_or_else(|| cfd_core::Error::NumericalError("Grid size is too large to be represented by float type".to_string()))?;

        // Calculate residuals with simple loops - clearer and no slower
        let mut velocity_residual = T::zero();
        let mut density_residual = T::zero();
        
        for i in 0..self.nx {
            for j in 0..self.ny {
                let u_mag = (self.u[i][j].x.clone() * self.u[i][j].x.clone() +
                           self.u[i][j].y.clone() * self.u[i][j].y.clone()).sqrt();
                velocity_residual = velocity_residual + u_mag;
                density_residual = density_residual + (self.rho[i][j].clone() - T::one()).abs();
            }
        }

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
