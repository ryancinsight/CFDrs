//! CPU compute backend

use super::traits::{ComputeBackend, ComputeBuffer, ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use std::marker::PhantomData;

/// Safe conversion from f64 with fallback
fn safe_f64_to_t<T: RealField + Copy>(value: f64, fallback: T) -> T {
    T::from_f64(value).unwrap_or(fallback)
}

/// CPU buffer implementation
pub struct CpuBuffer<T: RealField + Copy> {
    data: Vec<T>,
}

impl<T: RealField + Copy> CpuBuffer<T> {
    /// Create a new CPU buffer
    #[must_use] pub fn new(size: usize) -> Self {
        Self {
            data: vec![T::zero(); size],
        }
    }

    /// Create buffer with initial data
    #[must_use] pub fn from_data(data: Vec<T>) -> Self {
        Self { data }
    }
}

impl<T: RealField + Copy> ComputeBuffer<T> for CpuBuffer<T> {
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

impl<T: RealField + Copy> std::fmt::Debug for CpuBuffer<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CpuBuffer")
            .field("size", &self.data.len())
            .finish()
    }
}

/// Example CPU kernel for advection
pub struct CpuAdvectionKernel<T: RealField + Copy> {
    _phantom: PhantomData<T>,
}

impl<T: RealField + Copy> Default for CpuAdvectionKernel<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> CpuAdvectionKernel<T> {
    /// Creates a new CPU advection kernel
    #[must_use] pub fn new() -> Self {
        Self {
            _phantom: PhantomData,
        }
    }
}

impl<T: RealField + Copy> ComputeKernel<T> for CpuAdvectionKernel<T> {
    fn name(&self) -> &'static str {
        "CPU Advection (Upwind)"
    }

    fn execute(&self, input: &[T], output: &mut [T], params: KernelParams) -> Result<()> {
        let (nx, ny, _nz) = params.domain_params.grid_dims;
        let (dx, dy, _dz) = params.domain_params.grid_spacing;
        let dt = params.domain_params.dt;

        // Simple upwind advection for demonstration
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                // Extract velocity from params - in real implementation this would come from velocity field
                // For demonstration, using reference velocity from domain params
                let vx = safe_f64_to_t(params.domain_params.reynolds * 0.001, T::zero());
                let vy = safe_f64_to_t(params.domain_params.reynolds * 0.0005, T::zero());

                // Upwind differences
                let gradient_x = if vx > T::zero() {
                    (input[idx] - input[idx - 1]) / safe_f64_to_t(dx, T::one())
                } else {
                    (input[idx + 1] - input[idx]) / safe_f64_to_t(dx, T::one())
                };

                let gradient_y = if vy > T::zero() {
                    (input[idx] - input[idx - nx]) / safe_f64_to_t(dy, T::one())
                } else {
                    (input[idx + nx] - input[idx]) / safe_f64_to_t(dy, T::one())
                };

                // Update
                output[idx] =
                    input[idx] - safe_f64_to_t(dt, T::one()) * (vx * gradient_x + vy * gradient_y);
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

impl<T: RealField + Copy> std::fmt::Debug for CpuAdvectionKernel<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CpuAdvectionKernel").finish()
    }
}
