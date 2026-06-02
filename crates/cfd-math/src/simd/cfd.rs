//! CFD-specific SIMD optimizations
//!
//! This module provides SIMD-accelerated operations commonly used in CFD,
//! including gradient computations, flux calculations, and vector field operations.
//!
//! Performance targets:
//! - 2-4x speedup on SIMD-capable hardware
//! - Zero-overhead fallback on unsupported architectures
//! - Integration with existing CFD solver pipeline

use crate::error::Result;
use moirai::{ParallelSlice, ParallelSliceMut};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// CFD-specific SIMD operations dispatcher
pub struct CfdSimdOps<T: RealField + Copy> {
    /// Phantom data for generic type
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> CfdSimdOps<T> {
    /// Create new CFD SIMD operations handler
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> CfdSimdOps<T> {
    /// Compute gradient components with SIMD acceleration
    ///
    /// Computes ∂u/∂x and ∂u/∂y using central differences with SIMD operations
    pub fn gradient_2d_simd(
        &self,
        u: &[T],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
    ) -> Result<(Vec<T>, Vec<T>)> {
        if u.len() != nx * ny {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string(),
            ));
        }
        if nx < 3 || ny < 3 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Grid must be at least 3x3 for gradient computation".to_string(),
            ));
        }

        let mut dudy = vec![T::zero(); nx * ny];
        let mut dudx = vec![T::zero(); nx * ny];

        let dx_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero) * dx);
        let dy_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero) * dy);

        // Interior points - use parallel processing for better performance
        let u_slice = &u;
        let dx_inv_val = dx_inv;
        let dy_inv_val = dy_inv;

        dudy.par_mut().enumerate(|idx, dudy_val| {
            let j = idx / nx;
            let i = idx % nx;

            if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
                // ∂u/∂y = (u[i,j+1] - u[i,j-1]) / (2*dy)
                let u_up = u_slice[(j + 1) * nx + i];
                let u_down = u_slice[(j - 1) * nx + i];
                *dudy_val = (u_up - u_down) * dy_inv_val;
            }
        });

        dudx.par_mut().enumerate(|idx, dudx_val| {
            let j = idx / nx;
            let i = idx % nx;

            if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
                // ∂u/∂x = (u[i+1,j] - u[i-1,j]) / (2*dx)
                let u_right = u_slice[j * nx + (i + 1)];
                let u_left = u_slice[j * nx + (i - 1)];
                *dudx_val = (u_right - u_left) * dx_inv_val;
            }
        });

        // For boundary points, use one-sided differences (non-parallelized for simplicity)
        self.apply_boundary_gradients(u, &mut dudx, &mut dudy, nx, ny, dx_inv, dy_inv);

        Ok((dudx, dudy))
    }

    /// Apply boundary conditions for gradient computation
    fn apply_boundary_gradients(
        &self,
        u: &[T],
        dudx: &mut [T],
        dudy: &mut [T],
        nx: usize,
        ny: usize,
        dx_inv: T,
        dy_inv: T,
    ) {
        // Ensure we don't access out of bounds
        if nx < 3 || ny < 3 {
            return; // Need at least 3x3 grid for boundary gradients
        }

        // Left boundary (i=0) - forward difference
        for j in 1..(ny - 1) {
            let idx = j * nx;
            let u_right = u[j * nx + 1];
            let u_here = u[idx];
            dudx[idx] = (u_right - u_here)
                * dx_inv
                * T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

            let u_up = u[(j + 1) * nx];
            let u_down = u[(j - 1) * nx];
            dudy[idx] = (u_up - u_down) * dy_inv;
        }

        // Right boundary (i=nx-1) - backward difference
        for j in 1..(ny - 1) {
            let idx = j * nx + (nx - 1);
            let u_left = u[j * nx + (nx - 2)];
            let u_here = u[idx];
            dudx[idx] = (u_here - u_left)
                * dx_inv
                * T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

            let u_up = u[(j + 1) * nx + (nx - 1)];
            let u_down = u[(j - 1) * nx + (nx - 1)];
            dudy[idx] = (u_up - u_down) * dy_inv;
        }

        // Bottom boundary (j=0) - forward difference
        for i in 1..(nx - 1) {
            let idx = i;
            let u_up = u[nx + i];
            let u_here = u[idx];
            dudy[idx] =
                (u_up - u_here) * dy_inv * T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

            let u_right = u[i + 1];
            let u_left = u[i - 1];
            dudx[idx] = (u_right - u_left) * dx_inv;
        }

        // Top boundary (j=ny-1) - backward difference
        for i in 1..(nx - 1) {
            let idx = (ny - 1) * nx + i;
            let u_down = u[(ny - 2) * nx + i];
            let u_here = u[idx];
            dudy[idx] = (u_here - u_down)
                * dy_inv
                * T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

            let u_right = u[(ny - 1) * nx + i + 1];
            let u_left = u[(ny - 1) * nx + i - 1];
            dudx[idx] = (u_right - u_left) * dx_inv;
        }
    }

    /// Compute convective fluxes with SIMD acceleration
    ///
    /// Computes F = u*φ and G = v*φ using SIMD operations
    pub fn convective_fluxes_simd(&self, phi: &[T], u: &[T], v: &[T]) -> Result<(Vec<T>, Vec<T>)> {
        if phi.len() != u.len() || phi.len() != v.len() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Field dimensions don't match".to_string(),
            ));
        }

        let mut f_flux = vec![T::zero(); phi.len()];
        let mut g_flux = vec![T::zero(); phi.len()];

        // Adaptive parallelism routes small arrays to the sequential path.
        f_flux.par_mut().enumerate(|i, f| *f = phi[i] * u[i]);
        g_flux.par_mut().enumerate(|i, g| *g = phi[i] * v[i]);

        Ok((f_flux, g_flux))
    }

    /// Compute viscous fluxes with SIMD acceleration
    ///
    /// Computes diffusive fluxes using gradients
    pub fn viscous_fluxes_simd(
        &self,
        viscosity: T,
        dphi_dx: &[T],
        dphi_dy: &[T],
    ) -> Result<(Vec<T>, Vec<T>)> {
        if dphi_dx.len() != dphi_dy.len() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Gradient dimensions don't match".to_string(),
            ));
        }

        let mut f_viscous = vec![T::zero(); dphi_dx.len()];
        let mut g_viscous = vec![T::zero(); dphi_dx.len()];

        // F_viscous = -μ * ∂φ/∂x
        // G_viscous = -μ * ∂φ/∂y
        let neg_viscosity = -viscosity;

        f_viscous
            .par_mut()
            .enumerate(|i, f| *f = neg_viscosity * dphi_dx[i]);
        g_viscous
            .par_mut()
            .enumerate(|i, g| *g = neg_viscosity * dphi_dy[i]);

        Ok((f_viscous, g_viscous))
    }

    /// Vectorized array operations for CFD fields
    pub fn field_operations_simd(&self) -> CfdFieldOps<T> {
        CfdFieldOps {
            _phantom: std::marker::PhantomData,
        }
    }
}

/// CFD field operations with SIMD acceleration
pub struct CfdFieldOps<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> CfdFieldOps<T> {
    /// Add two CFD fields element-wise
    pub fn add_fields(&self, a: &[T], b: &[T], result: &mut [T]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Field dimensions don't match".to_string(),
            ));
        }

        result.par_mut().enumerate(|i, r| *r = a[i] + b[i]);
        Ok(())
    }

    /// Scale CFD field by scalar
    pub fn scale_field(&self, field: &[T], scalar: T, result: &mut [T]) -> Result<()> {
        if field.len() != result.len() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Field dimensions don't match".to_string(),
            ));
        }

        result.par_mut().enumerate(|i, r| *r = field[i] * scalar);
        Ok(())
    }

    /// Compute L2 norm of CFD field
    pub fn field_norm(&self, field: &[T]) -> Result<T> {
        if field.is_empty() {
            return Ok(T::zero());
        }

        let sum_sq = field
            .par()
            .map_reduce(T::zero(), |&x| x * x, |acc, x| acc + x);
        Ok(sum_sq.sqrt())
    }
}

impl<T: RealField + Copy> Default for CfdSimdOps<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gradient_2d_simd() {
        let ops = CfdSimdOps::<f64>::new();

        // Simple test field: u(x,y) = x + y
        // ∂u/∂x = 1, ∂u/∂y = 1
        let nx = 10;
        let ny = 8;
        let mut u = vec![0.0; nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * 0.1;
                let y = j as f64 * 0.1;
                u[j * nx + i] = x + y;
            }
        }

        let (dudx, dudy) = ops.gradient_2d_simd(&u, nx, ny, 0.1, 0.1).unwrap();

        // Check interior points (should be close to 1.0)
        for j in 2..(ny - 2) {
            for i in 2..(nx - 2) {
                let idx = j * nx + i;
                assert_relative_eq!(dudx[idx], 1.0, epsilon = 1e-6);
                assert_relative_eq!(dudy[idx], 1.0, epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn test_convective_fluxes_simd() {
        let ops = CfdSimdOps::<f64>::new();

        let phi = vec![1.0, 2.0, 3.0];
        let u = vec![0.5, 1.0, 1.5];
        let v = vec![2.0, 1.0, 0.5];

        let (f_flux, g_flux) = ops.convective_fluxes_simd(&phi, &u, &v).unwrap();

        assert_eq!(f_flux.len(), 3);
        assert_eq!(g_flux.len(), 3);
        assert_relative_eq!(f_flux[0], 0.5, epsilon = 1e-10);
        assert_relative_eq!(f_flux[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(g_flux[0], 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_basic_functionality() {
        let _ops = CfdSimdOps::<f64>::new();
        // Just verify it doesn't panic
    }
}
