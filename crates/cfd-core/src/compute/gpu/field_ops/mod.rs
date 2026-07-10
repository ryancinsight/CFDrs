//! GPU field-operation facade.

mod arithmetic;
mod laplacian;

#[cfg(test)]
mod tests;

use crate::compute::gpu::kernels::Laplacian2DKernel;
use crate::compute::gpu::GpuContext;
use crate::error::Result;
use std::sync::Arc;

/// Facade for GPU-accelerated field operations.
pub struct GpuFieldOps {
    pub(super) context: Arc<GpuContext>,
    pub(super) laplacian_kernel: Laplacian2DKernel,
}

impl GpuFieldOps {
    /// Create a field-operation facade over one GPU context.
    ///
    /// # Errors
    /// Returns a typed provider error when the Laplacian kernel cannot be
    /// compiled for the context.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        Ok(Self {
            context: context.clone(),
            laplacian_kernel: Laplacian2DKernel::new(context)?,
        })
    }
}
