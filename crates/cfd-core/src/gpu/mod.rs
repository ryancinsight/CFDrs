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
pub use crate::compute::gpu::GpuContext;
