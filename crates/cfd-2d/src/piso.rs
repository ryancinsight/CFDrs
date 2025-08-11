//! PISO (Pressure-Implicit with Splitting of Operators) algorithm implementation
//!
//! The PISO algorithm is an extension of SIMPLE with additional corrector steps
//! for improved accuracy in transient calculations.

use crate::grid::StructuredGrid2D;
use cfd_core::{Result, Error};
use nalgebra::{Vector2, ComplexField, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Named constants for PISO algorithm
const DEFAULT_MAX_CORRECTORS: usize = 2;
const DEFAULT_VELOCITY_RELAXATION: f64 = 1.0; // No under-relaxation in PISO
const DEFAULT_PRESSURE_RELAXATION: f64 = 1.0;
const MOMENTUM_COEFFICIENT_FACTOR: f64 = 4.0;
const GRADIENT_FACTOR: f64 = 2.0;

/// PISO solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PisoConfig<T: RealField> {
    /// Base solver configuration
    pub base: cfd_core::SolverConfig<T>,
    /// Time step for transient simulation
    pub time_step: T,
    /// Number of corrector steps (typically 2)
    pub num_correctors: usize,
    /// Convergence tolerance for velocity
    pub velocity_tolerance: T,
    /// Convergence tolerance for pressure
    pub pressure_tolerance: T,
    /// Non-orthogonal correction iterations
    pub non_orthogonal_correctors: usize,
}

impl<T: RealField + FromPrimitive> Default for PisoConfig<T> {
    fn default() -> Self {
        let base = cfd_core::SolverConfig::builder()
            .max_iterations(50) // PISO needs fewer iterations than SIMPLE
            .tolerance(T::from_f64(1e-6).unwrap())
            .build();

        Self {
            base,
            time_step: T::from_f64(0.01).unwrap(), // Default time step
            num_correctors: DEFAULT_MAX_CORRECTORS,
            velocity_tolerance: T::from_f64(1e-6).unwrap(),
            pressure_tolerance: T::from_f64(1e-6).unwrap(),
            non_orthogonal_correctors: 1,
        }
    }
}

/// PISO (Pressure-Implicit with Splitting of Operators) solver
///
/// ## Algorithm Steps:
/// 1. **Momentum Predictor**: Solve momentum equations with current pressure
/// 2. **Pressure Corrector 1**: First pressure correction from continuity
/// 3. **Velocity Corrector 1**: Update velocities with first pressure correction
/// 4. **Pressure Corrector 2**: Second pressure correction for improved accuracy
/// 5. **Velocity Corrector 2**: Final velocity update
///
/// ## Key Differences from SIMPLE:
/// - No under-relaxation needed (fully implicit)
/// - Multiple pressure corrections improve transient accuracy
/// - Better suited for time-dependent problems
///
/// ## References:
/// Issa, R.I. (1986). "Solution of the implicitly discretised fluid flow equations
/// by operator-splitting." Journal of Computational Physics, 62(1), 40-65.
pub struct PisoSolver<T: RealField> {
    config: PisoConfig<T>,
    /// Velocity field [nx][ny]
    u: Vec<Vec<Vector2<T>>>,
    /// Pressure field [nx][ny]
    p: Vec<Vec<T>>,
    /// First pressure correction [nx][ny]
    p_prime: Vec<Vec<T>>,
    /// Second pressure correction [nx][ny]
    p_double_prime: Vec<Vec<T>>,
    /// Velocity corrections
    u_prime: Vec<Vec<Vector2<T>>>,
    u_double_prime: Vec<Vec<Vector2<T>>>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    /// Fluid properties
    density: T,
    viscosity: T,
    /// Momentum equation coefficients
    a_p: Vec<Vec<T>>,
    a_e: Vec<Vec<T>>,
    a_w: Vec<Vec<T>>,
    a_n: Vec<Vec<T>>,
    a_s: Vec<Vec<T>>,
}

impl<T: RealField + FromPrimitive + Send + Sync> PisoSolver<T> {
    /// Create a new PISO solver
    pub fn new(
        config: PisoConfig<T>,
        grid: &StructuredGrid2D<T>,
        density: T,
        viscosity: T,
    ) -> Self {
        let nx = grid.nx;
        let ny = grid.ny;

        Self {
            config,
            u: vec![vec![Vector2::zeros(); ny]; nx],
            p: vec![vec![T::zero(); ny]; nx],
            p_prime: vec![vec![T::zero(); ny]; nx],
            p_double_prime: vec![vec![T::zero(); ny]; nx],
            u_prime: vec![vec![Vector2::zeros(); ny]; nx],
            u_double_prime: vec![vec![Vector2::zeros(); ny]; nx],
            nx,
            ny,
            dx: grid.dx.clone(),
            dy: grid.dy.clone(),
            density,
            viscosity,
            a_p: vec![vec![T::zero(); ny]; nx],
            a_e: vec![vec![T::zero(); ny]; nx],
            a_w: vec![vec![T::zero(); ny]; nx],
            a_n: vec![vec![T::zero(); ny]; nx],
            a_s: vec![vec![T::zero(); ny]; nx],
        }
    }

    /// Initialize solver with initial conditions
    pub fn initialize(&mut self, initial_velocity: Vector2<T>, initial_pressure: T) -> Result<()> {
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.u[i][j] = initial_velocity.clone();
                self.p[i][j] = initial_pressure.clone();
            }
        }
        self.compute_momentum_coefficients();
        Ok(())
    }

    /// Compute momentum equation coefficients
    fn compute_momentum_coefficients(&mut self) {
        let dt = self.config.time_step.clone();
        let dx2 = self.dx.clone() * self.dx.clone();
        let dy2 = self.dy.clone() * self.dy.clone();
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                // Diffusion coefficients
                let diff_x = self.viscosity.clone() / dx2.clone();
                let diff_y = self.viscosity.clone() / dy2.clone();
                
                // Convection coefficients (using central differencing)
                let u_e = if i+1 < self.nx {
                    (self.u[i][j].x.clone() + self.u[i+1][j].x.clone()) / T::from_f64(GRADIENT_FACTOR).unwrap()
                } else {
                    self.u[i][j].x.clone() / T::from_f64(GRADIENT_FACTOR).unwrap()
                };
                let u_w = (self.u[i-1][j].x.clone() + self.u[i][j].x.clone()) / T::from_f64(GRADIENT_FACTOR).unwrap();
                let v_n = (self.u[i][j].y.clone() + self.u[i][j+1].y.clone()) / T::from_f64(GRADIENT_FACTOR).unwrap();
                let v_s = (self.u[i][j-1].y.clone() + self.u[i][j].y.clone()) / T::from_f64(GRADIENT_FACTOR).unwrap();
                
                let conv_e = self.density.clone() * u_e / self.dx.clone();
                let conv_w = self.density.clone() * u_w / self.dx.clone();
                let conv_n = self.density.clone() * v_n / self.dy.clone();
                let conv_s = self.density.clone() * v_s / self.dy.clone();
                
                // Hybrid scheme coefficients
                let zero = T::zero();
                self.a_e[i][j] = diff_x.clone() + if -conv_e.clone() > zero { -conv_e.clone() } else { zero.clone() };
                self.a_w[i][j] = diff_x.clone() + if conv_w.clone() > zero { conv_w.clone() } else { zero.clone() };
                self.a_n[i][j] = diff_y.clone() + if -conv_n.clone() > zero { -conv_n.clone() } else { zero.clone() };
                self.a_s[i][j] = diff_y.clone() + if conv_s.clone() > zero { conv_s } else { zero.clone() };
                
                // Central coefficient (includes transient term)
                self.a_p[i][j] = self.a_e[i][j].clone() + self.a_w[i][j].clone() 
                    + self.a_n[i][j].clone() + self.a_s[i][j].clone()
                    + self.density.clone() * self.dx.clone() * self.dy.clone() / dt.clone();
            }
        }
    }

    /// Solve momentum predictor step
    fn solve_momentum_predictor(&mut self) -> Result<()> {
        // Store old velocities
        let u_old = self.u.clone();
        
        // Solve momentum equations with current pressure field
        for _ in 0..10 { // Simple iteration
            for i in 1..self.nx-1 {
                for j in 1..self.ny-1 {
                    // X-momentum
                    let sum_nb_x = self.a_e[i][j].clone() * self.u[i+1][j].x.clone()
                        + self.a_w[i][j].clone() * self.u[i-1][j].x.clone()
                        + self.a_n[i][j].clone() * self.u[i][j+1].x.clone()
                        + self.a_s[i][j].clone() * self.u[i][j-1].x.clone();
                    
                    let pressure_gradient_x = (self.p[i-1][j].clone() - self.p[i+1][j].clone()) 
                        * self.dy.clone() / T::from_f64(GRADIENT_FACTOR).unwrap();
                    
                    let source_x = self.density.clone() * self.dx.clone() * self.dy.clone() 
                        * u_old[i][j].x.clone() / self.config.time_step.clone();
                    
                    self.u[i][j].x = (sum_nb_x + pressure_gradient_x + source_x) / self.a_p[i][j].clone();
                    
                    // Y-momentum
                    let sum_nb_y = self.a_e[i][j].clone() * self.u[i+1][j].y.clone()
                        + self.a_w[i][j].clone() * self.u[i-1][j].y.clone()
                        + self.a_n[i][j].clone() * self.u[i][j+1].y.clone()
                        + self.a_s[i][j].clone() * self.u[i][j-1].y.clone();
                    
                    let pressure_gradient_y = (self.p[i][j-1].clone() - self.p[i][j+1].clone()) 
                        * self.dx.clone() / T::from_f64(GRADIENT_FACTOR).unwrap();
                    
                    let source_y = self.density.clone() * self.dx.clone() * self.dy.clone() 
                        * u_old[i][j].y.clone() / self.config.time_step.clone();
                    
                    self.u[i][j].y = (sum_nb_y + pressure_gradient_y + source_y) / self.a_p[i][j].clone();
                }
            }
        }
        
        Ok(())
    }

    /// Solve pressure correction equation
    fn solve_pressure_correction(&mut self, use_double_prime: bool) -> Result<()> {
        // Choose which correction field to use
        let correction_field = if use_double_prime {
            &mut self.p_double_prime
        } else {
            &mut self.p_prime
        };
        
        // Reset correction field
        for i in 0..self.nx {
            for j in 0..self.ny {
                correction_field[i][j] = T::zero();
            }
        }
        
        // Solve pressure correction equation iteratively
        for _ in 0..20 {
            for i in 1..self.nx-1 {
                for j in 1..self.ny-1 {
                    // Calculate mass imbalance (continuity error)
                    let dudx = (self.u[i+1][j].x.clone() - self.u[i-1][j].x.clone()) 
                        / (T::from_f64(GRADIENT_FACTOR).unwrap() * self.dx.clone());
                    let dvdy = (self.u[i][j+1].y.clone() - self.u[i][j-1].y.clone()) 
                        / (T::from_f64(GRADIENT_FACTOR).unwrap() * self.dy.clone());
                    let mass_imbalance = self.density.clone() * (dudx + dvdy);
                    
                    // Pressure correction coefficients
                    let d_e = self.dy.clone() / self.a_p[i][j].clone();
                    let d_w = self.dy.clone() / self.a_p[i][j].clone();
                    let d_n = self.dx.clone() / self.a_p[i][j].clone();
                    let d_s = self.dx.clone() / self.a_p[i][j].clone();
                    
                    let a_p_prime = d_e.clone() / self.dx.clone() + d_w.clone() / self.dx.clone()
                        + d_n.clone() / self.dy.clone() + d_s.clone() / self.dy.clone();
                    
                    // Solve for pressure correction
                    let sum_nb = (d_e * correction_field[i+1][j].clone() 
                        + d_w * correction_field[i-1][j].clone()) / self.dx.clone()
                        + (d_n * correction_field[i][j+1].clone() 
                        + d_s * correction_field[i][j-1].clone()) / self.dy.clone();
                    
                    correction_field[i][j] = (sum_nb - mass_imbalance) / a_p_prime;
                }
            }
        }
        
        Ok(())
    }

    /// Apply velocity correction
    fn apply_velocity_correction(&mut self, use_double_prime: bool) {
        let pressure_correction = if use_double_prime {
            &self.p_double_prime
        } else {
            &self.p_prime
        };
        
        let velocity_correction = if use_double_prime {
            &mut self.u_double_prime
        } else {
            &mut self.u_prime
        };
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                // Velocity correction from pressure gradient
                let dp_dx = (pressure_correction[i+1][j].clone() - pressure_correction[i-1][j].clone())
                    / (T::from_f64(GRADIENT_FACTOR).unwrap() * self.dx.clone());
                let dp_dy = (pressure_correction[i][j+1].clone() - pressure_correction[i][j-1].clone())
                    / (T::from_f64(GRADIENT_FACTOR).unwrap() * self.dy.clone());
                
                let d_coeff = self.dx.clone() * self.dy.clone() / self.a_p[i][j].clone();
                
                velocity_correction[i][j].x = -d_coeff.clone() * dp_dx;
                velocity_correction[i][j].y = -d_coeff * dp_dy;
                
                // Update velocity
                self.u[i][j] += velocity_correction[i][j].clone();
            }
        }
    }

    /// Execute one PISO iteration
    pub fn step(&mut self) -> Result<()> {
        // Step 1: Momentum predictor
        self.solve_momentum_predictor()?;
        
        // Step 2: First pressure correction
        self.solve_pressure_correction(false)?;
        
        // Step 3: First velocity correction
        self.apply_velocity_correction(false);
        
        // Update pressure with first correction
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.p[i][j] += self.p_prime[i][j].clone();
            }
        }
        
        // Step 4: Second pressure correction (PISO specific)
        if self.config.num_correctors >= 2 {
            self.solve_pressure_correction(true)?;
            
            // Step 5: Second velocity correction
            self.apply_velocity_correction(true);
            
            // Update pressure with second correction
            for i in 0..self.nx {
                for j in 0..self.ny {
                    self.p[i][j] += self.p_double_prime[i][j].clone();
                }
            }
        }
        
        // Apply boundary conditions
        self.apply_boundary_conditions();
        
        Ok(())
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&mut self) {
        // Left and right walls (no-slip)
        for j in 0..self.ny {
            self.u[0][j] = Vector2::zeros();
            self.u[self.nx-1][j] = Vector2::zeros();
        }
        
        // Top and bottom walls (no-slip)
        for i in 0..self.nx {
            self.u[i][0] = Vector2::zeros();
            self.u[i][self.ny-1] = Vector2::zeros();
        }
        
        // Pressure boundary conditions (zero gradient)
        for j in 0..self.ny {
            self.p[0][j] = self.p[1][j].clone();
            self.p[self.nx-1][j] = self.p[self.nx-2][j].clone();
        }
        for i in 0..self.nx {
            self.p[i][0] = self.p[i][1].clone();
            self.p[i][self.ny-1] = self.p[i][self.ny-2].clone();
        }
    }

    /// Check convergence
    pub fn check_convergence(&self) -> bool {
        let mut max_continuity_error = T::zero();
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                let dudx = (self.u[i+1][j].x.clone() - self.u[i-1][j].x.clone()) 
                    / (T::from_f64(GRADIENT_FACTOR).unwrap() * self.dx.clone());
                let dvdy = (self.u[i][j+1].y.clone() - self.u[i][j-1].y.clone()) 
                    / (T::from_f64(GRADIENT_FACTOR).unwrap() * self.dy.clone());
                let continuity_error = (dudx + dvdy).abs();
                
                if continuity_error > max_continuity_error {
                    max_continuity_error = continuity_error;
                }
            }
        }
        
        max_continuity_error < self.config.pressure_tolerance
    }

    /// Get velocity field
    pub fn velocity_field(&self) -> &Vec<Vec<Vector2<T>>> {
        &self.u
    }

    /// Get pressure field
    pub fn pressure_field(&self) -> &Vec<Vec<T>> {
        &self.p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_piso_creation() {
        let grid = StructuredGrid2D::<f64>::unit_square(5, 5).unwrap();
        let config = PisoConfig::default();
        let solver = PisoSolver::new(config, &grid, 1.0, 0.001);
        
        assert_eq!(solver.nx, 5);
        assert_eq!(solver.ny, 5);
    }

    #[test]
    fn test_piso_initialization() {
        let grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();
        let config = PisoConfig::default();
        let mut solver = PisoSolver::new(config, &grid, 1.0, 0.001);
        
        solver.initialize(Vector2::new(1.0, 0.5), 101325.0).unwrap();
        
        assert_relative_eq!(solver.u[1][1].x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(solver.u[1][1].y, 0.5, epsilon = 1e-10);
        assert_relative_eq!(solver.p[1][1], 101325.0, epsilon = 1e-10);
    }
}