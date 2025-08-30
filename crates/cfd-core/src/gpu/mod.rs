//! GPU acceleration using wgpu-rs for cross-platform compatibility
//!
//! Supports both integrated and discrete GPUs

#[cfg(feature = "gpu")]
pub mod compute;

#[cfg(feature = "gpu")]
pub mod field_ops;

#[cfg(feature = "gpu")]
pub use compute::GpuContext;

#[cfg(not(feature = "gpu"))]
pub struct GpuContext;

#[cfg(not(feature = "gpu"))]
impl GpuContext {
    pub fn new() -> Option<Self> {
        None
    }
}
