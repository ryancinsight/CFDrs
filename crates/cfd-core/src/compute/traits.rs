//! Core traits for compute operations

use crate::error::Result;
use nalgebra::RealField;
use std::fmt::Debug;

/// Trait for compute kernels that can run on different backends
pub trait ComputeKernel<T: RealField + Copy>: Debug + Send + Sync {
    /// Name of the kernel for identification
    fn name(&self) -> &str;

    /// Execute the kernel with given input and output buffers
    fn execute(&self, input: &[T], output: &mut [T], params: KernelParams) -> Result<()>;

    /// Execute kernel on GPU if implemented
    #[cfg(feature = "gpu")]
    fn execute_gpu(
        &self,
        _gpu_context: &std::sync::Arc<crate::compute::gpu::GpuContext>,
        _input: &[T],
        _output: &mut [T],
        _params: KernelParams,
    ) -> Result<()> {
        Err(crate::error::Error::UnsupportedOperation(
            "GPU execution not implemented for this kernel".to_string(),
        ))
    }

    /// Estimate computational complexity (FLOPs)
    fn complexity(&self, size: usize) -> usize;

    /// Check if this kernel can run on the given backend
    fn supports_backend(&self, backend: &ComputeBackend) -> bool;
}

/// Parameters for kernel execution
#[derive(Debug, Clone, Copy)]
pub struct KernelParams {
    /// Problem size (e.g., number of grid points)
    pub size: usize,
    /// Work group size for GPU/parallel execution
    pub work_group_size: usize,
    /// Additional domain-specific parameters
    pub domain_params: DomainParams,
}

/// Domain-specific parameters for CFD kernels
#[derive(Debug, Clone, Copy)]
pub struct DomainParams {
    /// Grid dimensions (nx, ny, nz)
    pub grid_dims: (usize, usize, usize),
    /// Grid spacing (dx, dy, dz)
    pub grid_spacing: (f64, f64, f64),
    /// Time step
    pub dt: f64,
    /// Reynolds number
    pub reynolds: f64,
}

/// Trait for compute buffers that can be shared between backends
pub trait ComputeBuffer<T: RealField + Copy>: Debug + Send + Sync {
    /// Get buffer size
    fn size(&self) -> usize;

    /// Read data from buffer
    fn read(&self) -> Result<Vec<T>>;

    /// Write data to buffer
    fn write(&mut self, data: &[T]) -> Result<()>;

    /// Map buffer for direct access (if supported)
    fn map(&self) -> Option<&[T]>;

    /// Map buffer for mutable access (if supported)
    fn map_mut(&mut self) -> Option<&mut [T]>;
}

/// Compute backend types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComputeBackend {
    /// CPU backend (scalar operations)
    Cpu,
    /// SIMD backend (vectorized CPU)
    Simd,
    /// GPU backend (wgpu)
    Gpu,
    /// Hybrid backend (CPU + GPU)
    Hybrid,
}

impl ComputeBackend {
    /// Check if backend is available on current system
    #[must_use]
    pub fn is_available(&self) -> bool {
        match self {
            Self::Cpu => true,
            Self::Simd => Self::detect_simd_support(),
            Self::Gpu => Self::detect_gpu_support(),
            Self::Hybrid => Self::detect_simd_support() || Self::detect_gpu_support(),
        }
    }

    /// Detect SIMD support
    fn detect_simd_support() -> bool {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            is_x86_feature_detected!("avx2") || is_x86_feature_detected!("sse4.1")
        }
        #[cfg(target_arch = "aarch64")]
        {
            std::arch::is_aarch64_feature_detected!("neon")
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
        {
            false
        }
    }

    /// Detect GPU support
    fn detect_gpu_support() -> bool {
        #[cfg(feature = "gpu")]
        {
            // Check if wgpu can create an adapter
            pollster::block_on(async {
                let instance = wgpu::Instance::default();
                instance
                    .request_adapter(&wgpu::RequestAdapterOptions::default())
                    .await
                    .is_some()
            })
        }
        #[cfg(not(feature = "gpu"))]
        {
            false
        }
    }
}
