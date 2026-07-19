//! GPU-accelerated matrix-free operators using Hephaestus kernels.

#[cfg(feature = "gpu")]
use crate::error::Result;
#[cfg(feature = "gpu")]
use crate::linear_solver::traits::{GpuLinearOperator, LinearOperator};
#[cfg(feature = "gpu")]
use aequitas::systems::si::quantities::Length;
#[cfg(feature = "gpu")]
pub use cfd_core::compute::gpu::kernels::laplacian::BoundaryCondition;
#[cfg(feature = "gpu")]
use cfd_core::compute::gpu::{kernels::laplacian::Laplacian2DKernel, GpuBuffer, GpuContext};
#[cfg(feature = "gpu")]
use cfd_core::error::Error;
#[cfg(feature = "gpu")]
use leto::Array1;
#[cfg(feature = "gpu")]
use std::sync::Arc;

/// Metrics collected for a compute dispatch.
#[cfg(feature = "gpu")]
#[derive(Debug, Clone, Copy)]
pub struct DispatchMetrics {
    /// Duration of the operation in milliseconds.
    pub duration_ms: f64,
    /// Whether the backend supports timestamp queries.
    pub timestamp_supported: bool,
}

/// GPU-accelerated two-dimensional Laplacian operator.
///
/// WGSL defines this operation over `f32`; exposing a generic scalar here would
/// claim unsupported precision and reinterpret non-`f32` buffers incorrectly.
#[cfg(feature = "gpu")]
pub struct GpuLaplacianOperator2D {
    gpu_context: Arc<GpuContext>,
    kernel: Laplacian2DKernel,
    nx: usize,
    ny: usize,
    dx: Length<f32>,
    dy: Length<f32>,
    bc: BoundaryCondition,
}

#[cfg(feature = "gpu")]
impl GpuLaplacianOperator2D {
    /// Create a GPU-accelerated two-dimensional Laplacian operator.
    ///
    /// # Errors
    /// Returns a typed provider error when the kernel cannot be compiled.
    pub fn new(
        gpu_context: Arc<GpuContext>,
        nx: usize,
        ny: usize,
        dx: Length<f32>,
        dy: Length<f32>,
        bc: BoundaryCondition,
    ) -> Result<Self> {
        let kernel = Laplacian2DKernel::new(gpu_context.clone())?;
        Ok(Self {
            gpu_context,
            kernel,
            nx,
            ny,
            dx,
            dy,
            bc,
        })
    }

    fn apply_host(&self, input: &[f32], output: &mut [f32]) -> Result<()> {
        let input_buffer = GpuBuffer::from_data(self.gpu_context.clone(), input)?;
        let mut output_buffer = GpuBuffer::new(self.gpu_context.clone(), output.len())?;
        self.kernel.execute_on_gpu(
            &input_buffer,
            &mut output_buffer,
            self.nx,
            self.ny,
            self.dx,
            self.dy,
            self.bc,
        )?;
        use cfd_core::compute::ComputeBuffer;
        output.copy_from_slice(&output_buffer.read()?);
        Ok(())
    }

    /// Execute the operator and return host-observed dispatch metrics.
    ///
    /// # Errors
    /// Returns a typed dimension, configuration, transfer, or provider error.
    pub fn apply_gpu_with_metrics(
        &self,
        input: &[f32],
        output: &mut [f32],
    ) -> Result<DispatchMetrics> {
        let start = std::time::Instant::now();
        self.apply_host(input, output)?;
        self.gpu_context.synchronize()?;
        let elapsed = start.elapsed();

        Ok(DispatchMetrics {
            duration_ms: elapsed.as_secs_f64() * 1e3,
            timestamp_supported: self.gpu_context.supports_timestamp_queries(),
        })
    }
}

#[cfg(feature = "gpu")]
impl LinearOperator<f32> for GpuLaplacianOperator2D {
    fn apply(&self, input: &Array1<f32>, output: &mut Array1<f32>) -> Result<()> {
        let expected = self.size();
        if input.shape()[0] != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: input.shape()[0],
            });
        }
        if output.shape()[0] != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: output.shape()[0],
            });
        }

        let input: Vec<f32> = (0..expected).map(|index| input[index]).collect();
        let mut values = vec![0.0; expected];
        self.apply_host(&input, &mut values)?;
        for (index, value) in values.into_iter().enumerate() {
            output[index] = value;
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
impl GpuLinearOperator<f32> for GpuLaplacianOperator2D {
    fn apply_gpu(
        &self,
        _gpu_context: &Arc<GpuContext>,
        input_buffer: &GpuBuffer<f32>,
        output_buffer: &mut GpuBuffer<f32>,
    ) -> Result<()> {
        self.kernel.execute_on_gpu(
            input_buffer,
            output_buffer,
            self.nx,
            self.ny,
            self.dx,
            self.dy,
            self.bc,
        )
    }
}
