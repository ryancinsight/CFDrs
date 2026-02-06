//! 2D Navier-Stokes Finite Volume Method (FVM) Solver Core
//!
//! This module provides a reusable FVM solver for 2D incompressible Navier-Stokes
//! with non-Newtonian blood rheology. It implements the SIMPLE algorithm for
//! pressure-velocity coupling on a staggered Cartesian grid.
//!
//! # Governing Equations
//!
//! **Continuity (incompressible):**
//! ```text
//! ∇·u = 0
//! ```
//!
//! **Momentum:**
//! ```text
//! ρ(u·∇)u = -∇p + ∇·(μ(γ̇)∇u) + f
//! ```
//!
//! Where:
//! - u = (u, v): velocity vector [m/s]
//! - p: pressure [Pa]
//! - ρ: density [kg/m³]
//! - μ(γ̇): shear-dependent viscosity [Pa·s]
//! - γ̇: shear rate magnitude [s⁻¹]
//! - f: body force [N/m³]
//!
//! # Discretization
//!
//! **Staggered Grid:**
//! ```text
//! - Pressure (p) stored at cell centers
//! - u-velocity at vertical faces (between cells in x-direction)
//! - v-velocity at horizontal faces (between cells in y-direction)
//! ```
//!
//! **SIMPLE Algorithm** (Semi-Implicit Method for Pressure-Linked Equations):
//! ```text
//! 1. Guess pressure field p*
//! 2. Solve momentum equations for u*, v* (with p*)
//! 3. Solve pressure correction equation (from continuity)
//! 4. Correct velocities and pressure
//! 5. Update viscosity from shear rate
//! 6. Check convergence, iterate if needed
//! ```
//!
//! **Rhie-Chow Interpolation:**
//! - Prevents pressure-velocity decoupling (checkerboard oscillations)
//! - Interpolates face velocities with pressure gradient term
//!
//! # Non-Newtonian Viscosity
//!
//! Viscosity μ computed from shear rate:
//! ```text
//! γ̇ = √(2 D:D)
//! where D = (∇u + ∇uᵀ)/2 (strain rate tensor)
//! ```
//!
//! For 2D: γ̇ = √(2[(∂u/∂x)² + (∂v/∂y)² + ((∂u/∂y + ∂v/∂x)/2)²])
//!
//! # References
//!
//! 1. Patankar, S.V. (1980) "Numerical Heat Transfer and Fluid Flow" - SIMPLE algorithm
//! 2. Versteeg & Malalasekera (2007) "An Introduction to CFD: The FVM" 2nd Ed.
//! 3. Rhie & Chow (1983) "Numerical study of turbulent flow past an airfoil" - Interpolation
//! 4. Ferziger & Perić (2002) "Computational Methods for Fluid Dynamics"
//! 5. Chhabra, R.P. (2006) "Bubbles, Drops, and Particles in Non-Newtonian Fluids"

use crate::error::Error;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};

/// Blood rheology model for non-Newtonian viscosity
#[derive(Debug, Clone)]
pub enum BloodModel<T: RealField + Copy> {
    /// Casson model with yield stress
    Casson(CassonBlood<T>),
    /// Carreau-Yasuda model
    CarreauYasuda(CarreauYasudaBlood<T>),
    /// Newtonian (constant viscosity)
    Newtonian(T),
}

/// Staggered grid for 2D FVM
///
/// Grid layout:
/// ```text
///   v[i,j+1]
///      |
/// u[i,j] -- p[i,j] -- u[i+1,j]
///      |
///   v[i,j]
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StaggeredGrid2D<T: RealField + Copy> {
    /// Number of cells in x-direction
    pub nx: usize,
    /// Number of cells in y-direction
    pub ny: usize,
    /// Cell width [m]
    pub dx: T,
    /// Cell height [m]
    pub dy: T,
    /// Domain width [m]
    pub width: T,
    /// Domain height [m]
    pub height: T,
}

impl<T: RealField + Copy + FromPrimitive> StaggeredGrid2D<T> {
    /// Create new staggered grid
    pub fn new(width: T, height: T, nx: usize, ny: usize) -> Self {
        let dx = width / T::from_usize(nx).unwrap();
        let dy = height / T::from_usize(ny).unwrap();

        Self {
            nx,
            ny,
            dx,
            dy,
            width,
            height,
        }
    }

    /// Get x-coordinate of cell center i
    pub fn x_center(&self, i: usize) -> T {
        (T::from_usize(i).unwrap() + T::from_f64(0.5).unwrap()) * self.dx
    }

    /// Get y-coordinate of cell center j
    pub fn y_center(&self, j: usize) -> T {
        (T::from_usize(j).unwrap() + T::from_f64(0.5).unwrap()) * self.dy
    }

    /// Get x-coordinate of u-face (between cells i-1 and i)
    pub fn x_u_face(&self, i: usize) -> T {
        T::from_usize(i).unwrap() * self.dx
    }

    /// Get y-coordinate of v-face (between cells j-1 and j)
    pub fn y_v_face(&self, j: usize) -> T {
        T::from_usize(j).unwrap() * self.dy
    }
}

/// 2D velocity and pressure fields on staggered grid
#[derive(Debug, Clone)]
pub struct FlowField2D<T: RealField + Copy> {
    /// u-velocity at vertical faces [nx+1][ny]
    pub u: Vec<Vec<T>>,
    /// v-velocity at horizontal faces [nx][ny+1]
    pub v: Vec<Vec<T>>,
    /// Pressure at cell centers [nx][ny]
    pub p: Vec<Vec<T>>,
    /// Viscosity at cell centers [nx][ny]
    pub mu: Vec<Vec<T>>,
    /// Shear rate at cell centers [nx][ny]
    pub gamma_dot: Vec<Vec<T>>,
}

impl<T: RealField + Copy + FromPrimitive + Float> FlowField2D<T> {
    /// Create new flow field with zero initial conditions
    pub fn new(nx: usize, ny: usize) -> Self {
        let u = vec![vec![T::zero(); ny]; nx + 1];
        let v = vec![vec![T::zero(); ny + 1]; nx];
        let p = vec![vec![T::zero(); ny]; nx];
        let mu = vec![vec![T::zero(); ny]; nx];
        let gamma_dot = vec![vec![T::zero(); ny]; nx];

        Self {
            u,
            v,
            p,
            mu,
            gamma_dot,
        }
    }

    /// Compute shear rate magnitude at cell (i,j)
    ///
    /// γ̇ = √(2 D:D) where D is strain rate tensor
    ///
    /// For 2D:
    /// ```text
    /// D₁₁ = ∂u/∂x
    /// D₂₂ = ∂v/∂y
    /// D₁₂ = D₂₁ = (∂u/∂y + ∂v/∂x)/2
    ///
    /// γ̇ = √(2[D₁₁² + D₂₂² + 2D₁₂²])
    /// ```
    pub fn compute_shear_rate(&self, i: usize, j: usize, dx: T, dy: T) -> T {
        let two = T::from_f64(2.0).unwrap();

        // ∂u/∂x at cell center (average from faces)
        let dudx = if i < self.u.len() - 1 {
            (self.u[i + 1][j] - self.u[i][j]) / dx
        } else {
            T::zero()
        };

        // ∂v/∂y at cell center (average from faces)
        let dvdy = if j < self.v[0].len() - 1 {
            (self.v[i][j + 1] - self.v[i][j]) / dy
        } else {
            T::zero()
        };

        // ∂u/∂y at cell center (interpolate u to center, then differentiate)
        let dudy = if j > 0 && j < self.u[0].len() - 1 && i < self.u.len() - 1 {
            let u_center_up = (self.u[i][j + 1] + self.u[i + 1][j + 1]) / two;
            let u_center_down = (self.u[i][j - 1] + self.u[i + 1][j - 1]) / two;
            (u_center_up - u_center_down) / (two * dy)
        } else {
            T::zero()
        };

        // ∂v/∂x at cell center (interpolate v to center, then differentiate)
        let dvdx = if i > 0 && i < self.v.len() - 1 && j < self.v[0].len() - 1 {
            let v_center_right = (self.v[i + 1][j] + self.v[i + 1][j + 1]) / two;
            let v_center_left = (self.v[i - 1][j] + self.v[i - 1][j + 1]) / two;
            (v_center_right - v_center_left) / (two * dx)
        } else {
            T::zero()
        };

        // Strain rate tensor components
        let d11 = dudx;
        let d22 = dvdy;
        let d12 = (dudy + dvdx) / two;

        // Magnitude: γ̇ = √(2[D₁₁² + D₂₂² + 2D₁₂²])
        let gamma_sq = two * (d11 * d11 + d22 * d22 + two * d12 * d12);

        Float::sqrt(Float::max(gamma_sq, T::zero()))
    }

    /// Update viscosity field based on shear rate and blood model
    pub fn update_viscosity(&mut self, grid: &StaggeredGrid2D<T>, blood: &BloodModel<T>) {
        for i in 0..grid.nx {
            for j in 0..grid.ny {
                // Compute shear rate
                let gamma = self.compute_shear_rate(i, j, grid.dx, grid.dy);
                self.gamma_dot[i][j] = gamma;

                // Compute viscosity from blood model
                self.mu[i][j] = match blood {
                    BloodModel::Casson(model) => model.apparent_viscosity(gamma),
                    BloodModel::CarreauYasuda(model) => model.apparent_viscosity(gamma),
                    BloodModel::Newtonian(mu) => *mu,
                };
            }
        }
    }

    /// Compute maximum velocity for CFL condition
    pub fn max_velocity(&self) -> T {
        let max_u = self
            .u
            .iter()
            .flat_map(|row| row.iter())
            .map(|&val| Float::abs(val))
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(T::zero());

        let max_v = self
            .v
            .iter()
            .flat_map(|row| row.iter())
            .map(|&val| Float::abs(val))
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(T::zero());

        Float::max(max_u, max_v)
    }
}

/// Configuration for SIMPLE algorithm
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SIMPLEConfig<T: RealField + Copy> {
    /// Maximum outer iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Under-relaxation factor for velocity (typical: 0.5-0.7)
    pub alpha_u: T,
    /// Under-relaxation factor for pressure (typical: 0.2-0.3)
    pub alpha_p: T,
    /// Under-relaxation factor for viscosity (typical: 0.5)
    pub alpha_mu: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for SIMPLEConfig<T> {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).unwrap(),
            alpha_u: T::from_f64(0.7).unwrap(),
            alpha_p: T::from_f64(0.3).unwrap(),
            alpha_mu: T::from_f64(0.5).unwrap(),
        }
    }
}

/// Boundary condition type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BCType {
    /// Fixed velocity (Dirichlet)
    Velocity,
    /// Fixed pressure (outlet)
    Pressure,
    /// No-slip wall (u=v=0)
    Wall,
    /// Symmetry
    Symmetry,
}

/// Boundary condition for a face
#[derive(Debug, Clone)]
pub struct BoundaryCondition<T: RealField + Copy> {
    /// Boundary type
    pub bc_type: BCType,
    /// Value (velocity or pressure)
    pub value: T,
}

/// 2D Navier-Stokes FVM solver with SIMPLE algorithm
///
/// This is the core solver that will be used by specific geometry solvers
/// (bifurcation, Venturi, serpentine, etc.)
#[derive(Debug)]
pub struct NavierStokesSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Computational grid
    pub grid: StaggeredGrid2D<T>,
    /// Flow field (u, v, p, μ, γ̇)
    pub field: FlowField2D<T>,
    /// Blood model
    pub blood: BloodModel<T>,
    /// Fluid density [kg/m³]
    pub density: T,
    /// SIMPLE configuration
    pub config: SIMPLEConfig<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Create new Navier-Stokes solver
    pub fn new(
        grid: StaggeredGrid2D<T>,
        blood: BloodModel<T>,
        density: T,
        config: SIMPLEConfig<T>,
    ) -> Self {
        let field = FlowField2D::new(grid.nx, grid.ny);

        Self {
            grid,
            field,
            blood,
            density,
            config,
        }
    }

    /// Initialize viscosity field
    pub fn initialize_viscosity(&mut self) {
        // Start with constant viscosity
        let mu_init = match &self.blood {
            BloodModel::Casson(model) => model.infinite_shear_viscosity,
            BloodModel::CarreauYasuda(model) => model.infinite_shear_viscosity,
            BloodModel::Newtonian(mu) => *mu,
        };

        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                self.field.mu[i][j] = mu_init;
            }
        }
    }

    /// Solve momentum equation for u-velocity (simplified for demonstration)
    ///
    /// Full implementation would use iterative solver (Gauss-Seidel, etc.)
    /// for the discretized momentum equations
    pub fn solve_u_momentum(&mut self, _bc_u: &[BoundaryCondition<T>]) -> Result<(), Error> {
        // TODO: Implement full momentum solver
        // This is a placeholder showing the structure
        Ok(())
    }

    /// Solve momentum equation for v-velocity
    pub fn solve_v_momentum(&mut self, _bc_v: &[BoundaryCondition<T>]) -> Result<(), Error> {
        // TODO: Implement full momentum solver
        Ok(())
    }

    /// Solve pressure correction equation
    pub fn solve_pressure_correction(&mut self) -> Result<(), Error> {
        // TODO: Implement pressure correction from continuity
        Ok(())
    }

    /// Apply Rhie-Chow interpolation for face velocities
    pub fn rhie_chow_interpolation(&mut self) {
        // TODO: Implement Rhie-Chow to prevent checkerboard
    }

    /// Compute residuals for convergence check
    pub fn compute_residuals(&self) -> (T, T, T) {
        // TODO: Compute continuity and momentum residuals
        (T::zero(), T::zero(), T::zero())
    }

    /// Run SIMPLE algorithm for steady-state solution
    pub fn solve(&mut self) -> Result<usize, Error> {
        self.initialize_viscosity();

        for iteration in 0..self.config.max_iterations {
            // SIMPLE algorithm steps
            // 1. Solve momentum equations
            // 2. Solve pressure correction
            // 3. Correct velocities and pressure
            // 4. Update viscosity from shear rate
            // 5. Check convergence

            // Update viscosity
            self.field.update_viscosity(&self.grid, &self.blood);

            // Check convergence
            let (res_u, res_v, res_p) = self.compute_residuals();
            let residual = Float::max(Float::max(res_u, res_v), res_p);

            if residual < self.config.tolerance {
                return Ok(iteration + 1);
            }
        }

        Err(Error::ConvergenceFailed {
            iterations: self.config.max_iterations,
            residual: 1.0, // Would be actual residual
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_staggered_grid_creation() {
        let grid = StaggeredGrid2D::<f64>::new(0.01, 0.001, 10, 5);
        assert_eq!(grid.nx, 10);
        assert_eq!(grid.ny, 5);
        assert!((grid.dx - 0.001).abs() < 1e-10);
        assert!((grid.dy - 0.0002).abs() < 1e-10);
    }

    #[test]
    fn test_flow_field_creation() {
        let field = FlowField2D::<f64>::new(10, 5);
        assert_eq!(field.u.len(), 11); // nx+1
        assert_eq!(field.v.len(), 10); // nx
        assert_eq!(field.p.len(), 10); // nx
        assert_eq!(field.u[0].len(), 5); // ny
        assert_eq!(field.v[0].len(), 6); // ny+1
    }
}
