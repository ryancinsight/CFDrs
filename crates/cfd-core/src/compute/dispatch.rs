//! Runtime dispatch for compute operations

use super::backend::{BackendContext, ComputeCapability};
use super::traits::{ComputeBackend, ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use std::sync::Arc;

/// Compute dispatcher that selects appropriate backend at runtime
pub struct ComputeDispatcher {
    /// Backend context
    context: BackendContext,
    /// System capabilities
    capabilities: Arc<ComputeCapability>,
}

impl ComputeDispatcher {
    /// Create a new dispatcher with automatic backend selection
    pub fn new() -> Result<Self> {
        let capabilities = Arc::new(ComputeCapability::detect());
        let backend = capabilities.preferred_backend;
        let context = BackendContext::new(backend)?;

        Ok(Self {
            context,
            capabilities,
        })
    }

    /// Create dispatcher with specific backend
    pub fn with_backend(backend: ComputeBackend) -> Result<Self> {
        let capabilities = Arc::new(ComputeCapability::detect());
        let context = BackendContext::new(backend)?;

        Ok(Self {
            context,
            capabilities,
        })
    }

    /// Execute a kernel with automatic backend selection
    pub fn execute<T: RealField + Copy>(
        &self,
        kernel: &dyn ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        // Select backend based on problem size and kernel support
        let backend = if kernel.supports_backend(&self.context.backend) {
            self.context.backend
        } else {
            // Fallback to CPU if kernel doesn't support current backend
            ComputeBackend::Cpu
        };

        match backend {
            ComputeBackend::Cpu => self.execute_cpu(kernel, input, output, params),
            ComputeBackend::Simd => self.execute_simd(kernel, input, output, params),
            ComputeBackend::Gpu => self.execute_gpu(kernel, input, output, params),
            ComputeBackend::Hybrid => self.execute_hybrid(kernel, input, output, params),
        }
    }

    /// Execute on CPU backend
    fn execute_cpu<T: RealField + Copy>(
        &self,
        kernel: &dyn ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        kernel.execute(input, output, params)
    }

    /// Execute on SIMD backend
    fn execute_simd<T: RealField + Copy>(
        &self,
        kernel: &dyn ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        // Check if we can use SIMD, otherwise fall back to CPU
        if ComputeBackend::Simd.is_available() {
            kernel.execute(input, output, params)
        } else {
            self.execute_cpu(kernel, input, output, params)
        }
    }

    /// Execute on GPU backend
    fn execute_gpu<T: RealField + Copy>(
        &self,
        kernel: &dyn ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        #[cfg(feature = "gpu")]
        {
            if let Some(ref gpu_context) = self.context.gpu_context {
                // GPU execution would go here
                // For now, fall back to CPU
                self.execute_cpu(kernel, input, output, params)
            } else {
                self.execute_cpu(kernel, input, output, params)
            }
        }
        #[cfg(not(feature = "gpu"))]
        {
            self.execute_cpu(kernel, input, output, params)
        }
    }

    /// Execute on hybrid backend (split between CPU and GPU)
    fn execute_hybrid<T: RealField + Copy>(
        &self,
        kernel: &dyn ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        // For large problems, split between GPU and CPU
        // For now, just use the best available backend
        if ComputeBackend::Gpu.is_available() {
            self.execute_gpu(kernel, input, output, params)
        } else if ComputeBackend::Simd.is_available() {
            self.execute_simd(kernel, input, output, params)
        } else {
            self.execute_cpu(kernel, input, output, params)
        }
    }

    /// Get current backend
    pub fn current_backend(&self) -> ComputeBackend {
        self.context.backend
    }

    /// Switch to a different backend
    pub fn switch_backend(&mut self, backend: ComputeBackend) -> Result<()> {
        self.context.switch_backend(backend)
    }

    /// Get system capabilities
    pub fn capabilities(&self) -> &ComputeCapability {
        &self.capabilities
    }
}
