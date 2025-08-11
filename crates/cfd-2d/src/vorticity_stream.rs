//! Vorticity-Stream Function formulation for 2D incompressible flow
//!
//! This solver uses the vorticity-stream function approach which automatically
//! satisfies continuity and reduces the number of variables.

use crate::grid::StructuredGrid2D;
use cfd_core::Result;
use nalgebra::{Vector2, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants
const DEFAULT_MAX_ITERATIONS: usize = 1000;
const DEFAULT_TOLERANCE: f64 = 1e-6;
const GRADIENT_FACTOR: f64 = 2.0;
const _LAPLACIAN_FACTOR: f64 = 4.0;
const SOR_OPTIMAL_FACTOR: f64 = 1.85; // Optimal for Poisson on square grid

/// Vorticity-Stream function solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VorticityStreamConfig<T: RealField> {
    /// Base solver configuration
    pub base: cfd_core::SolverConfig<T>,
    /// Time step for transient simulation
    pub time_step: T,
    /// Convergence tolerance for stream function
    pub stream_tolerance: T,
    /// Convergence tolerance for vorticity
    pub vorticity_tolerance: T,
    /// SOR relaxation parameter
    pub sor_omega: T,
}

impl<T: RealField + FromPrimitive> Default for VorticityStreamConfig<T> {
    fn default() -> Self {
        let base = cfd_core::SolverConfig::builder()
            .max_iterations(DEFAULT_MAX_ITERATIONS)
            .tolerance(T::from_f64(DEFAULT_TOLERANCE).unwrap())
            .build();

        Self {
            base,
            time_step: T::from_f64(0.001).unwrap(), // Default time step
            stream_tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap(),
            vorticity_tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap(),
            sor_omega: T::from_f64(SOR_OPTIMAL_FACTOR).unwrap(),
        }
    }
}

/// Vorticity-Stream Function solver for 2D incompressible flow
///
/// ## Governing Equations:
/// 
/// The Navier-Stokes equations in vorticity-stream function form:
/// - Stream function: ∇²ψ = -ω
/// - Vorticity transport: ∂ω/∂t + u·∇ω = ν∇²ω
/// - Velocity recovery: u = ∂ψ/∂y, v = -∂ψ/∂x
///
/// ## Advantages:
/// - Automatically satisfies continuity equation
/// - No pressure-velocity coupling issues
/// - Reduced number of variables (2 instead of 3)
///
/// ## References:
/// Anderson, J.D. (1995). "Computational Fluid Dynamics: The Basics with Applications."
/// McGraw-Hill.
pub struct VorticityStreamSolver<T: RealField> {
    config: VorticityStreamConfig<T>,
    /// Stream function field [nx][ny]
    psi: Vec<Vec<T>>,
    /// Vorticity field [nx][ny]
    omega: Vec<Vec<T>>,
    /// Velocity field (derived from stream function)
    u: Vec<Vec<Vector2<T>>>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    /// Reynolds number (Re = U*L/ν)
    reynolds: T,
}

impl<T: RealField + FromPrimitive + Send + Sync> VorticityStreamSolver<T> {
    /// Create a new vorticity-stream solver
    pub fn new(
        config: VorticityStreamConfig<T>,
        grid: &StructuredGrid2D<T>,
        reynolds: T,
    ) -> Self {
        let nx = grid.nx;
        let ny = grid.ny;

        Self {
            config,
            psi: vec![vec![T::zero(); ny]; nx],
            omega: vec![vec![T::zero(); ny]; nx],
            u: vec![vec![Vector2::zeros(); ny]; nx],
            nx,
            ny,
            dx: grid.dx.clone(),
            dy: grid.dy.clone(),
            reynolds,
        }
    }

    /// Initialize with lid-driven cavity conditions
    pub fn initialize_lid_driven_cavity(&mut self, lid_velocity: T) -> Result<()> {
        // Initialize stream function to zero
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.psi[i][j] = T::zero();
                self.omega[i][j] = T::zero();
            }
        }
        
        // Set top lid velocity boundary condition
        for i in 0..self.nx {
            self.u[i][self.ny-1].x = lid_velocity.clone();
            self.u[i][self.ny-1].y = T::zero();
        }
        
        Ok(())
    }

    /// Solve stream function Poisson equation: ∇²ψ = -ω
    fn solve_stream_function(&mut self) -> Result<()> {
        let omega_sor = self.config.sor_omega.clone();
        let dx2 = self.dx.clone() * self.dx.clone();
        let dy2 = self.dy.clone() * self.dy.clone();
        let denominator = T::from_f64(GRADIENT_FACTOR).unwrap() * (T::one() / dx2.clone() + T::one() / dy2.clone());
        
        // SOR iteration for Poisson equation
        for _ in 0..100 {
            let mut max_change = T::zero();
            
            for i in 1..self.nx-1 {
                for j in 1..self.ny-1 {
                    let psi_old = self.psi[i][j].clone();
                    
                    // Five-point stencil for Laplacian
                    let psi_new = ((self.psi[i+1][j].clone() + self.psi[i-1][j].clone()) / dx2.clone()
                        + (self.psi[i][j+1].clone() + self.psi[i][j-1].clone()) / dy2.clone()
                        + self.omega[i][j].clone()) / denominator.clone();
                    
                    // SOR update
                    self.psi[i][j] = (T::one() - omega_sor.clone()) * psi_old.clone() 
                        + omega_sor.clone() * psi_new;
                    
                    let change = (self.psi[i][j].clone() - psi_old).abs();
                    if change > max_change {
                        max_change = change;
                    }
                }
            }
            
            if max_change < self.config.stream_tolerance {
                break;
            }
        }
        
        Ok(())
    }

    /// Solve vorticity transport equation
    fn solve_vorticity_transport(&mut self) -> Result<()> {
        let dt = self.config.time_step.clone();
        let dx2 = self.dx.clone() * self.dx.clone();
        let dy2 = self.dy.clone() * self.dy.clone();
        let _two_dx = T::from_f64(GRADIENT_FACTOR).unwrap() * self.dx.clone();
        let _two_dy = T::from_f64(GRADIENT_FACTOR).unwrap() * self.dy.clone();
        
        // Store old vorticity
        let omega_old = self.omega.clone();
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                // Convection terms (upwind differencing for stability)
                let u = self.u[i][j].x.clone();
                let v = self.u[i][j].y.clone();
                
                let dwdx = if u >= T::zero() {
                    (self.omega[i][j].clone() - self.omega[i-1][j].clone()) / self.dx.clone()
                } else {
                    (self.omega[i+1][j].clone() - self.omega[i][j].clone()) / self.dx.clone()
                };
                
                let dwdy = if v >= T::zero() {
                    (self.omega[i][j].clone() - self.omega[i][j-1].clone()) / self.dy.clone()
                } else {
                    (self.omega[i][j+1].clone() - self.omega[i][j].clone()) / self.dy.clone()
                };
                
                let convection = u * dwdx + v * dwdy;
                
                // Diffusion term (central differencing)
                let d2wdx2 = (omega_old[i+1][j].clone() - T::from_f64(GRADIENT_FACTOR).unwrap() * omega_old[i][j].clone() 
                    + omega_old[i-1][j].clone()) / dx2.clone();
                let d2wdy2 = (omega_old[i][j+1].clone() - T::from_f64(GRADIENT_FACTOR).unwrap() * omega_old[i][j].clone() 
                    + omega_old[i][j-1].clone()) / dy2.clone();
                
                let diffusion = (d2wdx2 + d2wdy2) / self.reynolds.clone();
                
                // Time integration (explicit Euler)
                self.omega[i][j] = omega_old[i][j].clone() + dt.clone() * (diffusion - convection);
            }
        }
        
        Ok(())
    }

    /// Update velocity field from stream function
    fn update_velocity(&mut self) {
        let two_dx = T::from_f64(GRADIENT_FACTOR).unwrap() * self.dx.clone();
        let two_dy = T::from_f64(GRADIENT_FACTOR).unwrap() * self.dy.clone();
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                // u = ∂ψ/∂y
                self.u[i][j].x = (self.psi[i][j+1].clone() - self.psi[i][j-1].clone()) / two_dy.clone();
                // v = -∂ψ/∂x
                self.u[i][j].y = -(self.psi[i+1][j].clone() - self.psi[i-1][j].clone()) / two_dx.clone();
            }
        }
    }

    /// Update boundary vorticity using Thom's formula
    fn update_boundary_vorticity(&mut self) {
        let dx2 = self.dx.clone() * self.dx.clone();
        let dy2 = self.dy.clone() * self.dy.clone();
        let two = T::from_f64(GRADIENT_FACTOR).unwrap();
        
        // Bottom wall (y = 0)
        for i in 1..self.nx-1 {
            self.omega[i][0] = -two.clone() * self.psi[i][1].clone() / dy2.clone();
        }
        
        // Top wall (y = L) - moving lid
        for i in 1..self.nx-1 {
            let lid_velocity = self.u[i][self.ny-1].x.clone();
            self.omega[i][self.ny-1] = -two.clone() * self.psi[i][self.ny-2].clone() / dy2.clone()
                - two.clone() * lid_velocity * self.dy.clone() / dy2.clone();
        }
        
        // Left wall (x = 0)
        for j in 1..self.ny-1 {
            self.omega[0][j] = -two.clone() * self.psi[1][j].clone() / dx2.clone();
        }
        
        // Right wall (x = L)
        for j in 1..self.ny-1 {
            self.omega[self.nx-1][j] = -two.clone() * self.psi[self.nx-2][j].clone() / dx2.clone();
        }
    }

    /// Execute one time step
    pub fn step(&mut self) -> Result<()> {
        // Step 1: Update boundary vorticity
        self.update_boundary_vorticity();
        
        // Step 2: Solve vorticity transport equation
        self.solve_vorticity_transport()?;
        
        // Step 3: Solve stream function Poisson equation
        self.solve_stream_function()?;
        
        // Step 4: Update velocity field
        self.update_velocity();
        
        Ok(())
    }

    /// Check convergence based on change in stream function
    pub fn check_convergence(&self, psi_old: &Vec<Vec<T>>) -> bool {
        let mut max_change = T::zero();
        
        for i in 0..self.nx {
            for j in 0..self.ny {
                let change = (self.psi[i][j].clone() - psi_old[i][j].clone()).abs();
                if change > max_change {
                    max_change = change;
                }
            }
        }
        
        max_change < self.config.stream_tolerance
    }

    /// Get stream function field
    pub fn stream_function(&self) -> &Vec<Vec<T>> {
        &self.psi
    }

    /// Get vorticity field
    pub fn vorticity(&self) -> &Vec<Vec<T>> {
        &self.omega
    }

    /// Get velocity field
    pub fn velocity_field(&self) -> &Vec<Vec<Vector2<T>>> {
        &self.u
    }

    /// Calculate stream function at center for validation
    pub fn stream_at_center(&self) -> T {
        self.psi[self.nx/2][self.ny/2].clone()
    }

    /// Calculate vorticity at center for validation
    pub fn vorticity_at_center(&self) -> T {
        self.omega[self.nx/2][self.ny/2].clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vorticity_stream_creation() {
        let grid = StructuredGrid2D::<f64>::unit_square(10, 10).unwrap();
        let config = VorticityStreamConfig::default();
        let solver = VorticityStreamSolver::new(config, &grid, 100.0);
        
        assert_eq!(solver.nx, 10);
        assert_eq!(solver.ny, 10);
        assert_relative_eq!(solver.reynolds, 100.0, epsilon = 1e-10);
    }

    #[test]
    fn test_lid_driven_cavity_initialization() {
        let grid = StructuredGrid2D::<f64>::unit_square(5, 5).unwrap();
        let config = VorticityStreamConfig::default();
        let mut solver = VorticityStreamSolver::new(config, &grid, 100.0);
        
        solver.initialize_lid_driven_cavity(1.0).unwrap();
        
        // Check lid velocity
        assert_relative_eq!(solver.u[2][4].x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(solver.u[2][4].y, 0.0, epsilon = 1e-10);
        
        // Check initial stream function and vorticity
        assert_relative_eq!(solver.psi[2][2], 0.0, epsilon = 1e-10);
        assert_relative_eq!(solver.omega[2][2], 0.0, epsilon = 1e-10);
    }
}