//! GPU acceleration using wgpu-rs for cross-platform compatibility
//!
//! Supports both integrated and discrete GPUs

#[cfg(feature = "gpu")]
pub mod field_ops;

#[cfg(feature = "gpu")]
pub mod kernels;

#[cfg(feature = "gpu")]
pub mod shaders;

#[cfg(feature = "gpu")]
pub mod validation_tests;

// Re-export GpuContext from compute module
#[cfg(feature = "gpu")]
pub use crate::compute::gpu::GpuContext;

#[cfg(not(feature = "gpu"))]
/// GPU computation context (no-op implementation when GPU feature is disabled).
///
/// This is a zero-sized type that provides a compile-time guarantee that
/// GPU operations are not available when the GPU feature is disabled.
/// Any attempt to use GPU functionality will result in a compile-time error
/// rather than runtime failure.
///
/// # Design Rationale
///
/// Rather than providing stub methods that fail at runtime, this implementation
/// uses the type system to prevent GPU usage at compile time. Code that requires
/// GPU must be gated with `#[cfg(feature = "gpu")]`.
///
/// # Usage
///
/// ```rust,no_run
/// # use cfd_core::gpu::GpuContext;
/// // This code only compiles when the "gpu" feature is enabled
/// #[cfg(feature = "gpu")]
/// fn use_gpu(ctx: &GpuContext) {
///     // GPU operations here
/// }
/// ```
///
/// # References
///
/// - Rust API Guidelines: C-FEATURE (feature flags for optional functionality)
/// - Zero-cost abstractions: No runtime overhead for disabled features
#[derive(Debug, Clone, Copy, Default)]
pub struct GpuContext;

#[cfg(not(feature = "gpu"))]
impl GpuContext {
    /// Attempt to create a GPU context.
    ///
    /// # Errors
    ///
    /// Always returns an error when the GPU feature is disabled, with a
    /// clear message indicating that the feature must be enabled.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use cfd_core::gpu::GpuContext;
    /// let result = GpuContext::create();
    /// assert!(result.is_err());
    /// ```
    pub fn create() -> crate::error::Result<Self> {
        Err(crate::error::Error::InvalidConfiguration(
            "GPU support is not enabled. Rebuild with --features gpu to enable GPU acceleration."
                .to_string(),
        ))
    }
}

#[cfg(all(test, not(feature = "gpu")))]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_context_disabled() {
        // When GPU feature is disabled, create() should fail with clear error
        let result = GpuContext::create();
        assert!(result.is_err());

        if let Err(crate::error::Error::InvalidConfiguration(msg)) = result {
            assert!(msg.contains("GPU support is not enabled"));
            assert!(msg.contains("--features gpu"));
        } else {
            panic!("Expected InvalidConfiguration error");
        }
    }

    #[test]
    fn test_gpu_context_is_zero_sized() {
        // GpuContext should be zero-sized when GPU is disabled
        assert_eq!(std::mem::size_of::<GpuContext>(), 0);
    }

    #[test]
    fn test_gpu_context_default() {
        // GpuContext should implement Default for convenient initialization
        let _ctx = GpuContext::default();
        // If we get here, Default is implemented correctly
    }

    #[test]
    fn test_gpu_context_copy() {
        // GpuContext should be Copy (zero-sized type)
        let ctx1 = GpuContext::default();
        let ctx2 = ctx1; // Copy, not move
        let _ctx3 = ctx1; // Should still work (proves Copy trait)
        let _ctx4 = ctx2; // Should still work (proves Copy trait)
    }
}
