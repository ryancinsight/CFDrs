//! Runtime dispatch for compute operations.

use super::backend::{BackendContext, ComputeCapability};
use super::traits::{ComputeBackend, ComputeKernel, KernelParams};
use crate::error::Result;
use eunomia::RealField;
use std::sync::Arc;

/// Compute dispatcher for an explicit backend context.
pub struct ComputeDispatcher {
    /// Backend context
    context: BackendContext,
    /// System capabilities
    capabilities: Arc<ComputeCapability>,
}

impl ComputeDispatcher {
    /// Create a new dispatcher with capability-selected backend preference.
    ///
    /// # Errors
    /// Returns error if no suitable compute backend is available on this system
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
    ///
    /// # Errors
    /// Returns error if the specified backend is not available on this system
    pub fn with_backend(backend: ComputeBackend) -> Result<Self> {
        let capabilities = Arc::new(ComputeCapability::detect());
        let context = BackendContext::new(backend)?;

        Ok(Self {
            context,
            capabilities,
        })
    }

    /// Execute a kernel on the dispatcher's backend.
    ///
    /// # Errors
    /// Returns error if kernel execution fails or if input/output buffers are invalid
    pub fn execute<T: RealField + Copy>(
        &self,
        kernel: &impl ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        if !kernel.supports_backend(&self.context.backend) {
            return Err(crate::error::Error::UnsupportedOperation(format!(
                "kernel '{}' does not support requested backend {:?}",
                kernel.name(),
                self.context.backend
            )));
        }

        match self.context.backend {
            ComputeBackend::Cpu => Self::execute_cpu(kernel, input, output, params),
            ComputeBackend::Simd => Self::execute_simd(kernel, input, output, params),
            ComputeBackend::Gpu => self.execute_gpu(kernel, input, output, params),
            ComputeBackend::Hybrid => self.execute_hybrid(kernel, input, output, params),
        }
    }

    /// Execute on CPU backend
    fn execute_cpu<T: RealField + Copy>(
        kernel: &impl ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        kernel.execute(input, output, params)
    }

    /// Execute on SIMD backend
    fn execute_simd<T: RealField + Copy>(
        kernel: &impl ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        if !ComputeBackend::Simd.is_available() {
            return Err(crate::error::Error::UnsupportedOperation(
                "SIMD backend is not available on this system".to_string(),
            ));
        }
        kernel.execute(input, output, params)
    }

    /// Execute on GPU backend
    fn execute_gpu<T: RealField + Copy>(
        &self,
        kernel: &impl ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        #[cfg(feature = "gpu")]
        {
            if let Some(ref gpu_context) = self.context.gpu_context {
                if !kernel.supports_backend(&ComputeBackend::Gpu) {
                    return Err(crate::error::Error::UnsupportedOperation(format!(
                        "kernel '{}' does not support requested backend {:?}",
                        kernel.name(),
                        ComputeBackend::Gpu
                    )));
                }
                kernel.execute_gpu(gpu_context, input, output, params)
            } else {
                Err(crate::error::Error::UnsupportedOperation(
                    "GPU backend was requested without an initialized Hephaestus context"
                        .to_string(),
                ))
            }
        }
        #[cfg(not(feature = "gpu"))]
        {
            Err(crate::error::Error::UnsupportedOperation(format!(
                "GPU backend was requested for kernel '{}' with input length {}, output \
                     length {}, and problem size {}, but the gpu feature is disabled",
                kernel.name(),
                input.len(),
                output.len(),
                params.size
            )))
        }
    }

    /// Execute on the strongest supported member of the hybrid backend.
    fn execute_hybrid<T: RealField + Copy>(
        &self,
        kernel: &impl ComputeKernel<T>,
        input: &[T],
        output: &mut [T],
        params: KernelParams,
    ) -> Result<()> {
        if kernel.supports_backend(&ComputeBackend::Gpu) && ComputeBackend::Gpu.is_available() {
            self.execute_gpu(kernel, input, output, params)
        } else if kernel.supports_backend(&ComputeBackend::Simd)
            && ComputeBackend::Simd.is_available()
        {
            Self::execute_simd(kernel, input, output, params)
        } else if kernel.supports_backend(&ComputeBackend::Cpu) {
            Self::execute_cpu(kernel, input, output, params)
        } else {
            Err(crate::error::Error::UnsupportedOperation(format!(
                "kernel '{}' does not support any available hybrid backend member",
                kernel.name()
            )))
        }
    }

    /// Get current backend
    #[must_use]
    pub fn current_backend(&self) -> ComputeBackend {
        self.context.backend
    }

    /// Switch to a different backend
    ///
    /// # Errors
    /// Returns an error if backend switching fails due to hardware incompatibility or initialization errors
    pub fn switch_backend(&mut self, backend: ComputeBackend) -> Result<()> {
        self.context.switch_backend(backend)
    }

    /// Get system capabilities
    #[must_use]
    pub fn capabilities(&self) -> &ComputeCapability {
        &self.capabilities
    }
}
