//! Unified compute backend for CPU, GPU, and SIMD operations
//!
//! This module provides a hardware-agnostic compute layer that automatically
//! dispatches operations to the most appropriate backend based on system capabilities.

pub mod backend;

/// Example backend abstraction pattern demonstrating Rust best practices.
///
/// This module shows the recommended pattern for compute backend abstraction:
/// - Zero-cost generics with trait bounds
/// - Backend abstraction over wgpu-rs
/// - Feature-gated GPU compilation
/// - Iterator-based operations for zero-copy
pub mod backend_example;

pub mod dispatch;
pub mod traits;

#[cfg(feature = "gpu")]
pub mod gpu;

pub mod cpu;
pub mod simd;

#[cfg(test)]
mod tests;

pub use backend::{BackendContext, ComputeCapability};
pub use dispatch::ComputeDispatcher;
pub use traits::{ComputeBackend, ComputeBuffer, ComputeKernel, DomainParams, KernelParams};
