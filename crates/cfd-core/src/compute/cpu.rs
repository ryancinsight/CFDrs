//! CPU compute backend

use super::traits::{ComputeBackend, ComputeBuffer, ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use std::marker::PhantomData;

/// CPU buffer implementation
pub struct CpuBuffer<T: RealField + Copy> {
    data: Vec<T>,
}

impl<T: RealField + Copy> CpuBuffer<T> {
    /// Create a new CPU buffer
    pub fn new(size: usize) -> Self {
        Self {
            data: vec![T::zero(); size],
        }
    }

    /// Create buffer with initial data
    pub fn new_with_data(data: Vec<T>) -> Self {
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

impl<T: RealField + Copy> CpuAdvectionKernel<T> {
    pub fn new() -> Self {
        Self {
            _phantom: PhantomData,
        }
    }
}

impl<T: RealField + Copy> ComputeKernel<T> for CpuAdvectionKernel<T> {
    fn name(&self) -> &str {
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
                let vx = T::from_f64(params.domain_params.reynolds * 0.001).unwrap_or_else(T::zero);
                let vy =
                    T::from_f64(params.domain_params.reynolds * 0.0005).unwrap_or_else(T::zero);

                // Upwind differences
                let du_dx = if vx > T::zero() {
                    (input[idx] - input[idx - 1]) / T::from_f64(dx).unwrap_or_else(T::one)
                } else {
                    (input[idx + 1] - input[idx]) / T::from_f64(dx).unwrap_or_else(T::one)
                };

                let du_dy = if vy > T::zero() {
                    (input[idx] - input[idx - nx]) / T::from_f64(dy).unwrap_or_else(T::one)
                } else {
                    (input[idx + nx] - input[idx]) / T::from_f64(dy).unwrap_or_else(T::one)
                };

                // Update
                output[idx] =
                    input[idx] - T::from_f64(dt).unwrap_or_else(T::one) * (vx * du_dx + vy * du_dy);
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
