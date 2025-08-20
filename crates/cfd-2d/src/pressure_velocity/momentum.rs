//! Momentum equation solver for STANDARD algorithm

use nalgebra::{RealField, Vector2, DVector};
use num_traits::FromPrimitive;
use cfd_math::{LinearSolver, ConjugateGradient, BiCGSTAB};
use crate::grid::StructuredGrid2D;
use crate::schemes::{SpatialScheme, FiniteDifference};
use super::coefficients::CellCoefficients;
use super::config::PressureVelocityConfig;

/// Momentum equation solver
pub struct MomentumSolver<T: RealField> {
    /// Grid
    grid: StructuredGrid2D<T>,
    /// Configuration
    config: PressureVelocityConfig<T>,
}

impl<T: RealField + FromPrimitive> MomentumSolver<T> {
    /// Create new momentum solver
    pub fn new(grid: StructuredGrid2D<T>, config: PressureVelocityConfig<T>) -> Self {
        Self { grid, config }
    }
    
    /// Solve momentum equations for predicted velocity
    pub fn solve_momentum(
        &self,
        u: &Vec<Vec<Vector2<T>>>,
        p: &Vec<Vec<T>>,
        bc: &cfd_core::boundary::BoundaryCondition<T>,
        nu: T,
    ) -> cfd_core::error::Result<Vec<Vec<Vector2<T>>>> {
        let nx = self.grid.nx();
        let ny = self.grid.ny();
        let dx = self.grid.dx();
        let dy = self.grid.dy();
        
        let mut u_star = u.clone();
        
        // Build coefficient matrix for each velocity component
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let coeffs = self.compute_momentum_coefficients(
                    i, j, u, nu, dx, dy
                );
                
                // Apply pressure gradient
                let dp_dx = (p[i+1][j] - p[i-1][j]) / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
                let dp_dy = (p[i][j+1] - p[i][j-1]) / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);
                
                // Solve for velocity components
                if self.config.implicit_momentum {
                    // Implicit solver - more stable
                    u_star[i][j] = self.solve_implicit_momentum(
                        coeffs, Vector2::new(dp_dx, dp_dy), u[i][j]
                    )?;
                } else {
                    // Explicit update
                    u_star[i][j] = self.solve_explicit_momentum(
                        coeffs, Vector2::new(dp_dx, dp_dy), u[i][j]
                    );
                }
            }
        }
        
        // Apply boundary conditions
        self.apply_boundary_conditions(&mut u_star, bc)?;
        
        Ok(u_star)
    }
    
    /// Compute momentum equation coefficients
    fn compute_momentum_coefficients(
        &self,
        i: usize,
        j: usize,
        u: &Vec<Vec<Vector2<T>>>,
        nu: T,
        dx: T,
        dy: T,
    ) -> CellCoefficients<T> {
        let mut coeffs = CellCoefficients::zero();
        
        // Diffusion coefficients
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        coeffs.ae = nu / dx2;
        coeffs.aw = nu / dx2;
        coeffs.an = nu / dy2;
        coeffs.as_ = nu / dy2;
        
        // Convection coefficients based on scheme
        match self.config.convection_scheme {
            SpatialScheme::CentralDifference => {
                // Central differencing (may cause oscillations)
                let ue = (u[i][j].x + u[i+1][j].x) / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                let uw = (u[i-1][j].x + u[i][j].x) / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                let vn = (u[i][j].y + u[i][j+1].y) / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                let vs = (u[i][j-1].y + u[i][j].y) / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                
                coeffs.ae = coeffs.ae - ue / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
                coeffs.aw = coeffs.aw + uw / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
                coeffs.an = coeffs.an - vn / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);
                coeffs.as_ = coeffs.as_ + vs / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);
            }
            SpatialScheme::FirstOrderUpwind => {
                // First-order upwind (stable but diffusive)
                let ue = u[i][j].x.max(T::zero());
                let uw = u[i-1][j].x.min(T::zero());
                let vn = u[i][j].y.max(T::zero());
                let vs = u[i][j-1].y.min(T::zero());
                
                coeffs.ae = coeffs.ae + ue.abs() / dx;
                coeffs.aw = coeffs.aw + uw.abs() / dx;
                coeffs.an = coeffs.an + vn.abs() / dy;
                coeffs.as_ = coeffs.as_ + vs.abs() / dy;
            }
            _ => {
                // Default to upwind for stability
                coeffs.ae = coeffs.ae + u[i][j].x.abs() / dx;
                coeffs.aw = coeffs.aw + u[i][j].x.abs() / dx;
                coeffs.an = coeffs.an + u[i][j].y.abs() / dy;
                coeffs.as_ = coeffs.as_ + u[i][j].y.abs() / dy;
            }
        }
        
        // Calculate central coefficient
        coeffs.calculate_ap();
        
        // Apply under-relaxation
        coeffs.apply_relaxation(self.config.alpha_u);
        
        coeffs
    }
    
    /// Solve implicit momentum equation
    fn solve_implicit_momentum(
        &self,
        coeffs: CellCoefficients<T>,
        pressure_grad: Vector2<T>,
        u_old: Vector2<T>,
    ) -> cfd_core::error::Result<Vector2<T>> {
        // For now, use explicit update
        // Full implicit solver would require matrix assembly
        Ok(self.solve_explicit_momentum(coeffs, pressure_grad, u_old))
    }
    
    /// Solve explicit momentum equation
    fn solve_explicit_momentum(
        &self,
        coeffs: CellCoefficients<T>,
        pressure_grad: Vector2<T>,
        u_old: Vector2<T>,
    ) -> Vector2<T> {
        let source = coeffs.su - pressure_grad;
        let u_new = (source + coeffs.d * u_old) / coeffs.ap;
        
        // Apply relaxation
        u_old * (T::one() - self.config.alpha_u) + u_new * self.config.alpha_u
    }
    
    /// Apply boundary conditions
    fn apply_boundary_conditions(
        &self,
        u: &mut Vec<Vec<Vector2<T>>>,
        bc: &cfd_core::boundary::BoundaryCondition<T>,
    ) -> cfd_core::error::Result<()> {
        let nx = self.grid.nx();
        let ny = self.grid.ny();
        
        // Apply boundary conditions based on type
        match bc {
            cfd_core::boundary::BoundaryCondition::Dirichlet(value) => {
                // Fixed velocity at boundaries
                for i in 0..nx {
                    u[i][0] = Vector2::new(*value, T::zero());
                    u[i][ny-1] = Vector2::new(*value, T::zero());
                }
                for j in 0..ny {
                    u[0][j] = Vector2::new(*value, T::zero());
                    u[nx-1][j] = Vector2::new(*value, T::zero());
                }
            }
            cfd_core::boundary::BoundaryCondition::Neumann(_) => {
                // Zero gradient at boundaries
                for i in 0..nx {
                    u[i][0] = u[i][1];
                    u[i][ny-1] = u[i][ny-2];
                }
                for j in 0..ny {
                    u[0][j] = u[1][j];
                    u[nx-1][j] = u[nx-2][j];
                }
            }
            _ => {
                // Default to no-slip walls
                for i in 0..nx {
                    u[i][0] = Vector2::zeros();
                    u[i][ny-1] = Vector2::zeros();
                }
                for j in 0..ny {
                    u[0][j] = Vector2::zeros();
                    u[nx-1][j] = Vector2::zeros();
                }
            }
        }
        
        Ok(())
    }
}