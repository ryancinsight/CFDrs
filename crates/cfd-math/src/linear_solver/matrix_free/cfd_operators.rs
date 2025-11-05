//! CFD-specific linear operators for matrix-free computations.
//!
//! This module provides implementations of common CFD operators that can
//! be used with matrix-free solvers. These operators represent discretized
//! PDEs without explicitly storing the coefficient matrices.

use super::operator::LinearOperator;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;

/// 2D Laplacian operator using finite differences.
///
/// This operator represents the discrete Laplacian -∇² on a uniform grid
/// without storing the sparse matrix. It's commonly used for pressure
/// projection in incompressible Navier-Stokes equations.
pub struct LaplacianOperator2D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
}

impl<T: RealField + Copy> LaplacianOperator2D<T> {
    /// Create a new 2D Laplacian operator.
    ///
    /// # Arguments
    ///
    /// * `nx` - Number of grid points in x-direction
    /// * `ny` - Number of grid points in y-direction
    /// * `dx` - Grid spacing in x-direction
    /// * `dy` - Grid spacing in y-direction
    pub fn new(nx: usize, ny: usize, dx: T, dy: T) -> Self {
        Self { nx, ny, dx, dy }
    }
}

impl<T: RealField + Copy + From<f64>> LinearOperator<T> for LaplacianOperator2D<T> {
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        if x.len() != self.size() || y.len() != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);

        // Interior points: -∇²x = (x[i-1] + x[i+1] - 2*x[i])/dx² + (x[j-1] + x[j+1] - 2*x[j])/dy²
        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x[idx];

                // Second derivatives
                let d2x_dx2 = (x[idx - 1] + x[idx + 1] - T::from_f64(2.0).unwrap() * center) * dx2_inv;
                let d2x_dy2 = (x[idx - self.nx] + x[idx + self.nx] - T::from_f64(2.0).unwrap() * center) * dy2_inv;

                y[idx] = -(d2x_dx2 + d2x_dy2); // Negative Laplacian
            }
        }

        // Boundary conditions (homogeneous Neumann for pressure)
        self.apply_boundary_conditions(y);

        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }

    fn is_symmetric(&self) -> bool {
        true // Laplacian is symmetric
    }

    fn is_positive_definite(&self) -> Option<bool> {
        Some(true) // Negative Laplacian is positive definite for pressure
    }
}

impl<T: RealField + Copy + From<f64>> LaplacianOperator2D<T> {
    /// Apply homogeneous Neumann boundary conditions.
    fn apply_boundary_conditions(&self, y: &mut [T]) {
        // Left and right boundaries: zero normal derivative
        for j in 0..self.ny {
            let left_idx = j * self.nx;
            let right_idx = j * self.nx + (self.nx - 1);

            // For Neumann BC, we set the boundary value to equal the adjacent interior value
            if self.nx > 1 {
                y[left_idx] = y[left_idx + 1];
                y[right_idx] = y[right_idx - 1];
            }
        }

        // Top and bottom boundaries: zero normal derivative
        for i in 0..self.nx {
            let bottom_idx = i;
            let top_idx = (self.ny - 1) * self.nx + i;

            // For Neumann BC, we set the boundary value to equal the adjacent interior value
            if self.ny > 1 {
                y[bottom_idx] = y[bottom_idx + self.nx];
                y[top_idx] = y[top_idx - self.nx];
            }
        }
    }
}

/// Simplified momentum operator for 1D convection-diffusion.
///
/// This operator represents the 1D momentum equation:
/// ∂u/∂t + u∂u/∂x = ν∂²u/∂x² - ∂p/∂x
/// Without explicit matrix storage.
pub struct MomentumOperator1D<T: RealField + Copy> {
    n: usize,
    dx: T,
    viscosity: T,
    convection_velocity: T,
}

impl<T: RealField + Copy> MomentumOperator1D<T> {
    /// Create a new 1D momentum operator.
    pub fn new(n: usize, dx: T, viscosity: T, convection_velocity: T) -> Self {
        Self { n, dx, viscosity, convection_velocity }
    }
}

impl<T: RealField + Copy + From<f64>> LinearOperator<T> for MomentumOperator1D<T> {
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        if x.len() != self.size() || y.len() != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dx_inv = T::one() / self.dx;

        // Interior points: ν*d²u/dx² - u*du/dx
        for i in 1..(self.n - 1) {
            let center = x[i];

            // Diffusion term: ν*d²u/dx²
            let diffusion = self.viscosity * (x[i - 1] + x[i + 1] - T::from_f64(2.0).unwrap() * center) * dx2_inv;

            // Convection term: -u*du/dx (upwind approximation)
            let du_dx = (center - x[i - 1]) * dx_inv; // Backward difference
            let convection = -self.convection_velocity * center * du_dx;

            y[i] = diffusion + convection;
        }

        // Boundary conditions (simplified)
        y[0] = x[0]; // Left boundary: no-slip or inlet
        y[self.n - 1] = x[self.n - 1]; // Right boundary: outlet

        Ok(())
    }

    fn size(&self) -> usize {
        self.n
    }

    fn is_symmetric(&self) -> bool {
        false // Convection makes it non-symmetric
    }
}

/// 3D Poisson operator for pressure projection in Navier-Stokes.
///
pub struct PoissonOperator3D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
}

impl<T: RealField + Copy> PoissonOperator3D<T> {
    /// Create a new 3D Poisson operator.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: T, dy: T, dz: T) -> Self {
        Self { nx, ny, nz, dx, dy, dz }
    }
}

impl<T: RealField + Copy + From<f64>> LinearOperator<T> for PoissonOperator3D<T> {
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        if x.len() != self.size() || y.len() != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let dz2_inv = T::one() / (self.dz * self.dz);

        // Interior points: -∇²x = d²x/dx² + d²x/dy² + d²x/dz²
        for k in 1..(self.nz - 1) {
            for j in 1..(self.ny - 1) {
                for i in 1..(self.nx - 1) {
                    let idx = k * self.ny * self.nx + j * self.nx + i;

                    // Second derivatives in all directions
                    let d2x_dx2 = (x[idx - 1] + x[idx + 1] - T::from_f64(2.0).unwrap() * x[idx]) * dx2_inv;
                    let d2x_dy2 = (x[idx - self.nx] + x[idx + self.nx] - T::from_f64(2.0).unwrap() * x[idx]) * dy2_inv;
                    let d2x_dz2 = (x[idx - self.nx * self.ny] + x[idx + self.nx * self.ny] - T::from_f64(2.0).unwrap() * x[idx]) * dz2_inv;

                    // Negative Laplacian
                    y[idx] = -(d2x_dx2 + d2x_dy2 + d2x_dz2);
                }
            }
        }

        // Apply homogeneous Neumann boundary conditions
        self.apply_boundary_conditions(y);
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    fn is_symmetric(&self) -> bool {
        true
    }

    fn is_positive_definite(&self) -> Option<bool> {
        Some(true)
    }
}

impl<T: RealField + Copy + From<f64>> PoissonOperator3D<T> {
    fn apply_boundary_conditions(&self, y: &mut [T]) {
        // Apply Neumann BC on all faces
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let idx = k * self.ny * self.nx + j * self.nx + i;

                    // Check if on boundary
                    let on_x_boundary = i == 0 || i == self.nx - 1;
                    let on_y_boundary = j == 0 || j == self.ny - 1;
                    let on_z_boundary = k == 0 || k == self.nz - 1;

                    // Skip if not on boundary
                    if !on_x_boundary && !on_y_boundary && !on_z_boundary {
                        continue;
                    }

                    // For Neumann BC, set boundary value equal to adjacent interior value
                    if on_x_boundary && i == 0 && self.nx > 1 {
                        y[idx] = y[idx + 1];
                    } else if on_x_boundary && i == self.nx - 1 && self.nx > 1 {
                        y[idx] = y[idx - 1];
                    } else if on_y_boundary && j == 0 && self.ny > 1 {
                        y[idx] = y[idx + self.nx];
                    } else if on_y_boundary && j == self.ny - 1 && self.ny > 1 {
                        y[idx] = y[idx - self.nx];
                    } else if on_z_boundary && k == 0 && self.nz > 1 {
                        y[idx] = y[idx + self.nx * self.ny];
                    } else if on_z_boundary && k == self.nz - 1 && self.nz > 1 {
                        y[idx] = y[idx - self.nx * self.ny];
                    }
                }
            }
        }
    }
}

/// 2D Momentum operator for Navier-Stokes equations.
///
pub struct MomentumOperator2D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    viscosity: T,
    /// Velocity field for nonlinear convection (nx*ny vector)
    velocity_u: Vec<T>,
    velocity_v: Vec<T>,
}

impl<T: RealField + Copy> MomentumOperator2D<T> {
    /// Create a new 2D momentum operator.
    pub fn new(
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        viscosity: T,
        velocity_u: Vec<T>,
        velocity_v: Vec<T>,
    ) -> Self {
        assert_eq!(velocity_u.len(), nx * ny);
        assert_eq!(velocity_v.len(), nx * ny);
        Self {
            nx, ny, dx, dy, viscosity,
            velocity_u, velocity_v
        }
    }
}

impl<T: RealField + Copy + From<f64>> LinearOperator<T> for MomentumOperator2D<T> {
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        if x.len() != self.size() || y.len() != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let dx_inv = T::one() / self.dx;
        let dy_inv = T::one() / self.dy;

        // Interior points: -ν∇²u + (u·∇)u + ∇p/ρ (simplified linearization)
        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x[idx];

                // Laplacian term: ν∇²u
                let laplacian = self.viscosity * (
                    (x[idx - 1] + x[idx + 1] - T::from_f64(2.0).unwrap() * center) * dx2_inv +
                    (x[idx - self.nx] + x[idx + self.nx] - T::from_f64(2.0).unwrap() * center) * dy2_inv
                );

                // Simplified convection term: u·∇u (central difference approximation)
                let ui = self.velocity_u[idx];
                let vi = self.velocity_v[idx];
                let du_dx = (x[idx + 1] - x[idx - 1]) * T::from_f64(0.5).unwrap() * dx_inv;
                let du_dy = (x[idx + self.nx] - x[idx - self.nx]) * T::from_f64(0.5).unwrap() * dy_inv;
                let convection = ui * du_dx + vi * du_dy;

                // Combine terms: result = -ν∇²u + convection
                y[idx] = -laplacian + convection;
            }
        }

        // Boundary conditions (simplified no-slip)
        self.apply_boundary_conditions(y);
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }

    fn is_symmetric(&self) -> bool {
        false // Convection makes it non-symmetric
    }
}

impl<T: RealField + Copy + From<f64>> MomentumOperator2D<T> {
    fn apply_boundary_conditions(&self, y: &mut [T]) {
        // Left and right boundaries (no-slip)
        for j in 0..self.ny {
            let idx_left = j * self.nx;
            let idx_right = j * self.nx + (self.nx - 1);
            y[idx_left] = T::zero();
            y[idx_right] = T::zero();
        }

        // Top and bottom boundaries (no-slip)
        for i in 0..self.nx {
            let idx_bottom = i;
            let idx_top = (self.ny - 1) * self.nx + i;
            y[idx_bottom] = T::zero();
            y[idx_top] = T::zero();
        }
    }
}

/// Energy equation operator for thermal transport.
///
pub struct EnergyOperator2D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    thermal_diffusivity: T,
    /// Velocity field for advection (nx*ny vector)
    velocity_u: Vec<T>,
    velocity_v: Vec<T>,
}

impl<T: RealField + Copy> EnergyOperator2D<T> {
    /// Create a new 2D energy operator.
    pub fn new(
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        thermal_diffusivity: T,
        velocity_u: Vec<T>,
        velocity_v: Vec<T>,
    ) -> Self {
        assert_eq!(velocity_u.len(), nx * ny);
        assert_eq!(velocity_v.len(), nx * ny);
        Self {
            nx, ny, dx, dy, thermal_diffusivity,
            velocity_u, velocity_v
        }
    }
}

impl<T: RealField + Copy + From<f64>> LinearOperator<T> for EnergyOperator2D<T> {
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        if x.len() != self.size() || y.len() != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let dx_inv = T::one() / self.dx;
        let dy_inv = T::one() / self.dy;

        // Interior points: α∇²T - (u·∇)T (energy equation without source terms)
        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x[idx];

                // Diffusion: α∇²T
                let diffusion = self.thermal_diffusivity * (
                    (x[idx - 1] + x[idx + 1] - T::from_f64(2.0).unwrap() * center) * dx2_inv +
                    (x[idx - self.nx] + x[idx + self.nx] - T::from_f64(2.0).unwrap() * center) * dy2_inv
                );

                // Advection: -(u·∇)T
                let ui = self.velocity_u[idx];
                let vi = self.velocity_v[idx];
                let d_t_dx = (x[idx + 1] - x[idx - 1]) * T::from_f64(0.5).unwrap() * dx_inv;
                let d_t_dy = (x[idx + self.nx] - x[idx - self.nx]) * T::from_f64(0.5).unwrap() * dy_inv;
                let advection = ui * d_t_dx + vi * d_t_dy;

                y[idx] = diffusion - advection;
            }
        }

        // Boundary conditions (simplified)
        self.apply_boundary_conditions(y);
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }

    fn is_symmetric(&self) -> bool {
        false // Advection makes it non-symmetric
    }
}

impl<T: RealField + Copy + From<f64>> EnergyOperator2D<T> {
    fn apply_boundary_conditions(&self, y: &mut [T]) {
        // Left boundary: fixed temperature T=0
        for j in 0..self.ny {
            y[j * self.nx] = T::zero();
        }

        // Right boundary: outflow
        for j in 0..self.ny {
            let idx = j * self.nx + (self.nx - 1);
            if self.nx > 1 {
                y[idx] = y[idx - 1];
            }
        }

        // Top and bottom: insulated (zero normal derivative)
        for i in 0..self.nx {
            let idx_bottom = i;
            let idx_top = (self.ny - 1) * self.nx + i;
            if self.ny > 1 {
                y[idx_bottom] = y[idx_bottom + self.nx];
                y[idx_top] = y[idx_top - self.nx];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_laplacian_operator_2d() {
        let op = LaplacianOperator2D::new(5, 5, 1.0, 1.0);
        assert_eq!(op.size(), 25);
        assert!(op.is_symmetric());
        assert_eq!(op.is_positive_definite(), Some(true));

        // Test with constant field (should give zero)
        let x = vec![1.0; 25];
        let mut y = vec![0.0; 25];
        op.apply(&x, &mut y).unwrap();

        // Interior points should be approximately zero
        for j in 1..4 {
            for i in 1..4 {
                let idx = j * 5 + i;
                assert!((y[idx] as f64).abs() < 1e-10, "Interior point {} should be near zero, got {}", idx, y[idx]);
            }
        }
    }

    #[test]
    fn test_momentum_operator_1d() {
        let op = MomentumOperator1D::new(10, 0.1, 0.01, 1.0);
        assert_eq!(op.size(), 10);
        assert!(!op.is_symmetric());

        // Test basic functionality
        let x = vec![1.0; 10];
        let mut y = vec![0.0; 10];
        op.apply(&x, &mut y).unwrap();

        // Should not panic and produce some result
        assert_eq!(y.len(), 10);
    }

    #[test]
    fn test_poisson_operator_3d() {
        let op = PoissonOperator3D::new(3, 3, 3, 1.0, 1.0, 1.0);
        assert_eq!(op.size(), 27);
        assert!(op.is_symmetric());
        assert_eq!(op.is_positive_definite(), Some(true));

        // Test with constant field
        let x = vec![2.0; 27];
        let mut y = vec![0.0; 27];
        op.apply(&x, &mut y).unwrap();

        // Interior points should be zero for constant field
        let center_idx = 1 * 9 + 1 * 3 + 1; // (k=1,j=1,i=1) in 3x3x3 grid
        assert!((y[center_idx] as f64).abs() < 1e-10, "Center should be near zero, got {}", y[center_idx]);
    }

    #[test]
    fn test_energy_operator_2d() {
        let velocity_u = vec![1.0; 25]; // Uniform x-velocity
        let velocity_v = vec![0.0; 25]; // Zero y-velocity
        let op = EnergyOperator2D::new(5, 5, 1.0, 1.0, 0.1, velocity_u, velocity_v);
        assert_eq!(op.size(), 25);
        assert!(!op.is_symmetric());

        let x = (0..25).map(|i| i as f64 * 0.1).collect::<Vec<_>>();
        let mut y = vec![0.0; 25];
        op.apply(&x, &mut y).unwrap();

        // Should produce some non-zero result
        assert!(y.iter().any(|&val| val.abs() > 0.0));
    }
}
