//! GPU-accelerated matrix-free operators using WGSL compute shaders.

#[cfg(feature = "gpu")]
use crate::error::Result;
#[cfg(feature = "gpu")]
use crate::linear_solver::traits::{GpuLinearOperator, LinearOperator};
#[cfg(feature = "gpu")]
pub use cfd_core::compute::gpu::kernels::laplacian::BoundaryType;
#[cfg(feature = "gpu")]
use cfd_core::compute::gpu::{kernels::laplacian::Laplacian2DKernel, GpuBuffer, GpuContext};
#[cfg(feature = "gpu")]
use cfd_core::compute::ComputeBuffer;
#[cfg(feature = "gpu")]
use cfd_core::error::Error;
#[cfg(feature = "gpu")]
use eunomia::{NumericElement, RealField};
#[cfg(feature = "gpu")]
use leto::Array1;
#[cfg(feature = "gpu")]
use std::sync::Arc;

/// Metrics collected for a compute dispatch
#[cfg(feature = "gpu")]
#[derive(Debug, Clone, Copy)]
pub struct DispatchMetrics {
    /// Duration of the operation in milliseconds
    pub duration_ms: f64,
    /// Whether the backend supports timestamp queries
    pub timestamp_supported: bool,
}

/// GPU-accelerated 2D Laplacian operator.
///
/// This operator wraps the mathematically rigorous implementation in `cfd-core`
/// to provide a `LinearOperator` interface for iterative solvers.
#[cfg(feature = "gpu")]
pub struct GpuLaplacianOperator2D<T: RealField + Copy + bytemuck::Pod + bytemuck::Zeroable> {
    gpu_context: Arc<GpuContext>,
    kernel: Laplacian2DKernel,
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    bc: BoundaryType,
}

#[cfg(feature = "gpu")]
impl<T: RealField + Copy + bytemuck::Pod + bytemuck::Zeroable + NumericElement>
    GpuLaplacianOperator2D<T>
{
    /// Create a new GPU-accelerated 2D Laplacian operator.
    pub fn new(
        gpu_context: Arc<GpuContext>,
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        bc_type: BoundaryType,
    ) -> Self {
        let kernel = Laplacian2DKernel::new(gpu_context.clone());

        Self {
            gpu_context,
            kernel,
            nx,
            ny,
            dx,
            dy,
            bc: bc_type,
        }
    }

    /// Internal method to apply the operator using CPU buffers (with GPU acceleration)
    fn apply_gpu(&self, x: &[T], y: &mut [T]) -> Result<()> {
        // Create GPU buffers from CPU data
        let input_buffer = GpuBuffer::from_data(self.gpu_context.clone(), x)?;
        let mut output_buffer = GpuBuffer::new(self.gpu_context.clone(), y.len())?;

        let dx_f32 = <T as NumericElement>::to_f64(self.dx) as f32;
        let dy_f32 = <T as NumericElement>::to_f64(self.dy) as f32;

        // Execute kernel
        self.kernel.execute_on_gpu(
            &input_buffer,
            &mut output_buffer,
            self.nx,
            self.ny,
            dx_f32,
            dy_f32,
            self.bc,
        )?;

        // Read back results
        let result = output_buffer.read()?;
        y.copy_from_slice(&result);
        Ok(())
    }

    /// Execute the operator with performance metrics
    pub fn apply_gpu_with_metrics(&self, x: &[T], y: &mut [T]) -> Result<DispatchMetrics> {
        let start = std::time::Instant::now();
        self.apply_gpu(x, y)?;
        self.gpu_context.synchronize()?;
        let elapsed = start.elapsed();

        Ok(DispatchMetrics {
            duration_ms: elapsed.as_secs_f64() * 1e3,
            timestamp_supported: self.gpu_context.supports_timestamp_queries(),
        })
    }
}

#[cfg(feature = "gpu")]
impl<T: RealField + Copy + bytemuck::Pod + bytemuck::Zeroable + NumericElement> LinearOperator<T>
    for GpuLaplacianOperator2D<T>
{
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        let expected = self.size();
        if x.shape()[0] != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: x.shape()[0],
            });
        }
        if y.shape()[0] != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: y.shape()[0],
            });
        }

        let input: Vec<T> = (0..expected).map(|idx| x[idx]).collect();
        let mut output = vec![<T as NumericElement>::ZERO; expected];
        self.apply_gpu(&input, &mut output)?;
        for idx in 0..expected {
            y[idx] = output[idx];
        }
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }

    fn is_symmetric(&self) -> bool {
        true
    }
}

#[cfg(feature = "gpu")]
impl<T: RealField + Copy + bytemuck::Pod + bytemuck::Zeroable + NumericElement> GpuLinearOperator<T>
    for GpuLaplacianOperator2D<T>
{
    fn apply_gpu(
        &self,
        _gpu_context: &Arc<GpuContext>,
        input_buffer: &GpuBuffer<T>,
        output_buffer: &mut GpuBuffer<T>,
    ) -> Result<()> {
        let dx_f32 = <T as NumericElement>::to_f64(self.dx) as f32;
        let dy_f32 = <T as NumericElement>::to_f64(self.dy) as f32;

        // Execute kernel directly on GPU buffers
        self.kernel.execute_on_gpu(
            input_buffer,
            output_buffer,
            self.nx,
            self.ny,
            dx_f32,
            dy_f32,
            self.bc,
        )
    }

    fn supports_gpu(&self) -> bool {
        true
    }
}
