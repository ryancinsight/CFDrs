//! Vectorization utilities for SIMD operations in CFD computations.
//!
//! This module provides vectorized operations that can take advantage of SIMD
//! instructions for improved performance in numerical computations.

pub mod operations;
pub mod stencil;

// Re-export main types
pub use operations::VectorizedOps;
pub use stencil::StencilOps;
