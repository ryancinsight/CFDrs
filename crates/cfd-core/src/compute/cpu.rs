//! CPU compute backend
//!
//! This module provides rigorously implemented CPU kernels with mathematically
//! documented behavior. The 2D upwind advection kernel implements the standard
//! first-order upwind discretization for the constant-velocity linear advection
//! equation, with support for periodic and Dirichlet-zero boundary conditions.
//!
//! References:
//! - R. J. LeVeque, "Finite Volume Methods for Hyperbolic Problems", Cambridge University Press.
//! - C. Hirsch, "Numerical Computation of Internal and External Flows", Vol. 2.

use super::traits::{BoundaryCondition2D, ComputeBackend, ComputeBuffer, ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::marker::PhantomData;

/// Safe conversion from f64 with fallback
fn safe_f64_to_t<T: RealField + Copy + FromPrimitive>(value: f64, fallback: T) -> T {
    <T as FromPrimitive>::from_f64(value).unwrap_or(fallback)
}

/// CPU buffer implementation
pub struct CpuBuffer<T: RealField + Copy> {
    data: Vec<T>,
}

impl<T: RealField + Copy + FromPrimitive> CpuBuffer<T> {
    /// Create a new CPU buffer
    #[must_use]
    pub fn new(size: usize) -> Self {
        Self {
            data: vec![T::zero(); size],
        }
    }

    /// Create buffer with initial data
    #[must_use]
    pub fn from_data(data: Vec<T>) -> Self {
        Self { data }
    }
}

impl<T: RealField + Copy + FromPrimitive> ComputeBuffer<T> for CpuBuffer<T> {
    fn size(&self) -> usize {
        self.data.len()
    }

    fn read(&self) -> Result<Vec<T>> {
        Ok(self.data.clone())
    }

    fn write(&mut self, data: &[T]) -> Result<()> {
        self.data.copy_from_slice(data);
        Ok(())
    }

    fn map(&self) -> Option<&[T]> {
        Some(&self.data)
    }

    fn map_mut(&mut self) -> Option<&mut [T]> {
        Some(&mut self.data)
    }
}

impl<T: RealField + Copy + FromPrimitive> std::fmt::Debug for CpuBuffer<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CpuBuffer")
            .field("size", &self.data.len())
            .finish()
    }
}

/// 2D first-order upwind advection kernel (constant velocity)
///
/// Implements the discrete update
/// `ϕ^{n+1}_{i,j} = ϕ^n_{i,j} - Δt ( u D_x ϕ^n_{i,j} + v D_y ϕ^n_{i,j} )`,
/// where `D_x` and `D_y` are upwind differences chosen by the sign of `u` and `v`.
///
/// Theorem (Exactness for linear fields under constant velocity):
/// For `ϕ(x,y) = a x + b y` and constant `(u,v)`, the upwind scheme produces the
/// exact one-step update `ϕ^{n+1} = ϕ^n - Δt (u a + v b)` on the discrete grid,
/// assuming uniform spacing and consistent boundary handling.
///
/// Stability: CFL condition `max(|u| Δt/Δx, |v| Δt/Δy) ≤ 1` is required for
/// monotonicity preservation and avoidance of spurious oscillations.
pub struct CpuAdvectionKernel<T: RealField + Copy + FromPrimitive> {
    _phantom: PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive> Default for CpuAdvectionKernel<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive> CpuAdvectionKernel<T> {
    /// Creates a new CPU advection kernel
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ComputeKernel<T> for CpuAdvectionKernel<T> {
    fn name(&self) -> &'static str {
        "CPU Advection (Upwind)"
    }

    fn execute(&self, input: &[T], output: &mut [T], params: KernelParams) -> Result<()> {
        let (nx, ny, nz) = params.domain_params.grid_dims;
        assert_eq!(nz, 1, "CpuAdvectionKernel is 2D; expected nz == 1");
        let (dx, dy, _dz) = params.domain_params.grid_spacing;
        let dt = params.domain_params.dt;
        let (u, v, _w) = params.domain_params.velocity;
        let bc = params.domain_params.boundary;

        let dx_t = safe_f64_to_t(dx, T::one());
        let dy_t = safe_f64_to_t(dy, T::one());
        let dt_t = safe_f64_to_t(dt, T::one());
        let u_t = safe_f64_to_t(u, T::zero());
        let v_t = safe_f64_to_t(v, T::zero());

        assert_eq!(input.len(), nx * ny, "input length must equal nx*ny");
        assert_eq!(output.len(), nx * ny, "output length must equal nx*ny");

        // Helper to fetch values with boundary handling
        let get = |ii: isize, jj: isize| -> T {
            match bc {
                BoundaryCondition2D::Periodic => {
                    let i_wrap = ((ii % nx as isize) + nx as isize) % nx as isize;
                    let j_wrap = ((jj % ny as isize) + ny as isize) % ny as isize;
                    input[(j_wrap as usize) * nx + (i_wrap as usize)]
                }
                BoundaryCondition2D::DirichletZero => {
                    if ii < 0 || ii >= nx as isize || jj < 0 || jj >= ny as isize {
                        T::zero()
                    } else {
                        input[(jj as usize) * nx + (ii as usize)]
                    }
                }
            }
        };

        for j in 0..ny as isize {
            for i in 0..nx as isize {
                let idx = (j as usize) * nx + (i as usize);

                // Upwind differences based on velocity sign
                let grad_x = if u_t > T::zero() {
                    (get(i, j) - get(i - 1, j)) / dx_t
                } else {
                    (get(i + 1, j) - get(i, j)) / dx_t
                };

                let grad_y = if v_t > T::zero() {
                    (get(i, j) - get(i, j - 1)) / dy_t
                } else {
                    (get(i, j + 1) - get(i, j)) / dy_t
                };

                output[idx] = get(i, j) - dt_t * (u_t * grad_x + v_t * grad_y);
            }
        }

        Ok(())
    }

    fn complexity(&self, size: usize) -> usize {
        // Approximately 10 FLOPs per grid point
        size * 10
    }

    fn supports_backend(&self, backend: &ComputeBackend) -> bool {
        matches!(backend, ComputeBackend::Cpu | ComputeBackend::Hybrid)
    }
}

impl<T: RealField + Copy + FromPrimitive> std::fmt::Debug for CpuAdvectionKernel<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CpuAdvectionKernel").finish()
    }
}
