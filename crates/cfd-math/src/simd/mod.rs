//! SIMD and SWAR implementations for vectorized operations
//!
//! This module provides architecture-aware SIMD implementations with SWAR fallbacks
//! for portable high-performance computing in CFD simulations.

pub mod arch_detect;
pub mod operations;
pub mod swar;

#[cfg(test)]
mod tests;

use std::arch::is_x86_feature_detected;

/// Architecture detection for SIMD support
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdCapability {
    /// AVX2 support (256-bit vectors)
    Avx2,
    /// SSE4.2 support (128-bit vectors)
    Sse42,
    /// ARM NEON support
    Neon,
    /// No SIMD, use SWAR
    Swar,
}

impl SimdCapability {
    /// Detect the best available SIMD capability at runtime
    #[inline]
    pub fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                return SimdCapability::Avx2;
            }
            if is_x86_feature_detected!("sse4.2") {
                return SimdCapability::Sse42;
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            // ARM NEON is always available on AArch64
            return SimdCapability::Neon;
        }

        // Fallback to SWAR for other architectures
        SimdCapability::Swar
    }
}

// Re-export main types
pub use operations::{SimdOps, VectorOps};
pub use swar::SwarOps;
