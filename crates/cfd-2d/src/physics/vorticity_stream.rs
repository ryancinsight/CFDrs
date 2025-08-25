//! Vorticity-Stream Function formulation for 2D incompressible flow
//!
//! This solver uses the vorticity-stream function approach which automatically
//! satisfies continuity and reduces the number of variables.

use crate::grid::StructuredGrid2D;
use cfd_core::error::Result;
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
pub struct VorticityStreamConfig<T: RealField + Copy> {
    /// Base solver configuration
    pub base: cfd_core::solver::SolverConfig<T>,
    /// Time step for transient simulation
    pub time_step: T,
    /// Convergence tolerance for stream function
    pub stream_tolerance: T,
    /// Convergence tolerance for vorticity
    pub vorticity_tolerance: T,
    /// SOR relaxation parameter
    pub sor_omega: T,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for VorticityStreamConfig<T> {
    fn default() -> Self {
        let base = cfd_core::solver::SolverConfig::builder()
            .max_iterations(DEFAULT_MAX_ITERATIONS)
            .tolerance(T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::zero()))
            .build();

        Self {
            base,
            time_step: T::from_f64(0.001).unwrap_or_else(|| T::zero()), // Default time step
            stream_tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::zero()),
            vorticity_tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::zero()),
            sor_omega: T::from_f64(SOR_OPTIMAL_FACTOR).unwrap_or_else(|| T::zero()),
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
/// Anderson, J.D. (1995). "Computational Fluid Dynamics: The Cores with Applications."
/// McGraw-Hill.
pub struct VorticityStreamSolver<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive + Send + Sync> VorticityStreamSolver<T> {
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
            dx: grid.dx,
            dy: grid.dy,
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
            self.u[i][self.ny-1].x = lid_velocity;
            self.u[i][self.ny-1].y = T::zero();
        }
        
        Ok(())
    }

    /// Solve stream function Poisson equation: ∇²ψ = -ω
    fn solve_stream_function(&mut self) -> Result<()> {
        let omega_sor = self.config.sor_omega;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let denominator = T::from_f64(GRADIENT_FACTOR).unwrap_or_else(|| T::zero()) * (T::one() / dx2 + T::one() / dy2);
        
        // SOR iteration for Poisson equation
        for _ in 0..100 {
            let mut max_change = T::zero();
            
            for i in 1..self.nx-1 {
                for j in 1..self.ny-1 {
                    let psi_old = self.psi[i][j];
                    
                    // Five-point stencil for Laplacian
                    let psi_new = ((self.psi[i+1][j] + self.psi[i-1][j]) / dx2
                        + (self.psi[i][j+1] + self.psi[i][j-1]) / dy2
                        + self.omega[i][j]) / denominator;
                    
                    // SOR update
                    self.psi[i][j] = (T::one() - omega_sor) * psi_old 
                        + omega_sor * psi_new;
                    
                    let change = (self.psi[i][j] - psi_old).abs();
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
        let dt = self.config.time_step;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        
        // Store old vorticity
        let omega_old = self.omega.clone();
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                // Convection terms (upwind differencing for stability)
                let u = self.u[i][j].x;
                let v = self.u[i][j].y;
                
                let dwdx = if u >= T::zero() {
                    (self.omega[i][j] - self.omega[i-1][j]) / self.dx
                } else {
                    (self.omega[i+1][j] - self.omega[i][j]) / self.dx
                };
                
                let dwdy = if v >= T::zero() {
                    (self.omega[i][j] - self.omega[i][j-1]) / self.dy
                } else {
                    (self.omega[i][j+1] - self.omega[i][j]) / self.dy
                };
                
                let convection = u * dwdx + v * dwdy;
                
                // Diffusion term (central differencing)
                let d2wdx2 = (omega_old[i+1][j] - T::from_f64(GRADIENT_FACTOR).unwrap_or_else(|| T::zero()) * omega_old[i][j] 
                    + omega_old[i-1][j]) / dx2;
                let d2wdy2 = (omega_old[i][j+1] - T::from_f64(GRADIENT_FACTOR).unwrap_or_else(|| T::zero()) * omega_old[i][j] 
                    + omega_old[i][j-1]) / dy2;
                
                let diffusion = (d2wdx2 + d2wdy2) / self.reynolds;
                
                // Time integration (explicit Euler)
                self.omega[i][j] = omega_old[i][j] + dt * (diffusion - convection);
            }
        }
        
        Ok(())
    }

    /// Update velocity field from stream function
    fn update_velocity(&mut self) {
        let two_dx = T::from_f64(GRADIENT_FACTOR).unwrap_or_else(|| T::zero()) * self.dx;
        let two_dy = T::from_f64(GRADIENT_FACTOR).unwrap_or_else(|| T::zero()) * self.dy;
        
        for i in 1..self.nx-1 {
            for j in 1..self.ny-1 {
                // u = ∂ψ/∂y
                self.u[i][j].x = (self.psi[i][j+1] - self.psi[i][j-1]) / two_dy;
                // v = -∂ψ/∂x
                self.u[i][j].y = -(self.psi[i+1][j] - self.psi[i-1][j]) / two_dx;
            }
        }
    }

    /// Update boundary vorticity using Thom's formula
    fn update_boundary_vorticity(&mut self) {
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let two = T::from_f64(GRADIENT_FACTOR).unwrap_or_else(|| T::zero());
        
        // Bottom wall (y = 0)
        for i in 1..self.nx-1 {
            self.omega[i][0] = -two * self.psi[i][1] / dy2;
        }
        
        // Top wall (y = L) - moving lid
        for i in 1..self.nx-1 {
            let lid_velocity = self.u[i][self.ny-1].x;
            self.omega[i][self.ny-1] = -two * self.psi[i][self.ny-2] / dy2
                - two * lid_velocity * self.dy / dy2;
        }
        
        // Left wall (x = 0)
        for j in 1..self.ny-1 {
            self.omega[0][j] = -two * self.psi[1][j] / dx2;
        }
        
        // Right wall (x = L)
        for j in 1..self.ny-1 {
            self.omega[self.nx-1][j] = -two * self.psi[self.nx-2][j] / dx2;
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
                let change = (self.psi[i][j] - psi_old[i][j]).abs();
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
        self.psi[self.nx/2][self.ny/2]
    }

    /// Calculate vorticity at center for validation
    pub fn vorticity_at_center(&self) -> T {
        self.omega[self.nx/2][self.ny/2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vorticity_stream_creation() {
        let grid = StructuredGrid2D::<f64>::unit_square(10, 10).expect("CRITICAL: Add proper error handling");
        let config = VorticityStreamConfig::default();
        let solver = VorticityStreamSolver::new(config, &grid, 100.0);
        
        assert_eq!(solver.nx, 10);
        assert_eq!(solver.ny, 10);
        assert_relative_eq!(solver.reynolds, 100.0, epsilon = 1e-10);
    }

    #[test]
    fn test_lid_driven_cavity_initialization() {
        let grid = StructuredGrid2D::<f64>::unit_square(5, 5).expect("CRITICAL: Add proper error handling");
        let config = VorticityStreamConfig::default();
        let mut solver = VorticityStreamSolver::new(config, &grid, 100.0);
        
        solver.initialize_lid_driven_cavity(1.0).expect("CRITICAL: Add proper error handling");
        
        // Check lid velocity
        assert_relative_eq!(solver.u[2][4].x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(solver.u[2][4].y, 0.0, epsilon = 1e-10);
        
        // Check initial stream function and vorticity
        assert_relative_eq!(solver.psi[2][2], 0.0, epsilon = 1e-10);
        assert_relative_eq!(solver.omega[2][2], 0.0, epsilon = 1e-10);
    }
}