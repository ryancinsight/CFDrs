//! Unified compute backend for CPU, GPU, and SIMD operations.
//!
//! This module provides backend capability detection plus explicit dispatch
//! through the selected provider. Unsupported provider/kernel combinations are
//! reported as typed errors instead of being recomputed on another backend.

pub mod backend;

pub mod cpu;
pub mod dispatch;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod simd;
pub mod solver;
pub mod time;
pub mod traits;

#[cfg(test)]
mod tests;

pub use backend::{BackendContext, ComputeCapability};
pub use dispatch::ComputeDispatcher;
pub use traits::{ComputeBackend, ComputeBuffer, ComputeKernel, DomainParams, KernelParams};
