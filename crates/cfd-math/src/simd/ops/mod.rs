//! SIMD operations module with architecture-specific implementations
//!
//! Organized following SSOT and modular architecture principles.

mod arm;
mod traits;
mod x86;

pub use traits::{SimdOperation, VectorOps};

use super::{SimdCapability, SwarOps};
use crate::error::Result;

/// Main SIMD operations dispatcher
pub struct SimdOps {
    capability: SimdCapability,
    swar: SwarOps,
}

impl SimdOps {
    /// Create new SIMD operations handler with detected capabilities
    pub fn new() -> Self {
        Self {
            capability: SimdCapability::detect(),
            swar: SwarOps::new(),
        }
    }

    /// Create with specific capability (for testing)
    pub fn with_capability(capability: SimdCapability) -> Self {
        Self {
            capability,
            swar: SwarOps::new(),
        }
    }

    /// Get current SIMD capability
    pub fn capability(&self) -> SimdCapability {
        self.capability
    }
}

impl Default for SimdOps {
    fn default() -> Self {
        Self::new()
    }
}

impl VectorOps for SimdOps {
    #[inline]
    fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::add_avx2_f32(a, b, result) },

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => unsafe { x86::add_sse42_f32(a, b, result) },

            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            SimdCapability::Neon => unsafe { arm::add_neon_f32(a, b, result) },

            _ => self.swar.add(a, b, result),
        }
    }

    #[inline]
    fn sub(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::sub_avx2_f32(a, b, result) },

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => unsafe { x86::sub_sse42_f32(a, b, result) },

            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            SimdCapability::Neon => unsafe { arm::sub_neon_f32(a, b, result) },

            _ => self.swar.sub(a, b, result),
        }
    }

    #[inline]
    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::mul_avx2_f32(a, b, result) },

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => unsafe { x86::mul_sse42_f32(a, b, result) },

            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            SimdCapability::Neon => unsafe { arm::mul_neon_f32(a, b, result) },

            _ => self.swar.mul(a, b, result),
        }
    }

    #[inline]
    fn div(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Division is typically not SIMD-optimized, use SWAR
        self.swar.div(a, b, result)
    }

    #[inline]
    fn fma(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != c.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "fma"))]
            SimdCapability::Avx2 => unsafe { x86::fma_avx2_f32(a, b, c, result) },

            _ => self.swar.fma(a, b, c, result),
        }
    }

    #[inline]
    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::scale_avx2_f32(input, scalar, result) },

            _ => self.swar.scale(input, scalar, result),
        }
    }

    #[inline]
    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        if a.len() != b.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => Ok(unsafe { x86::dot_avx2_f32(a, b) }),

            _ => self.swar.dot(a, b),
        }
    }

    // f64 implementations delegate to appropriate modules
    #[inline]
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::add_avx2_f64(a, b, result) },

            _ => self.swar.add_f64(a, b, result),
        }
    }

    #[inline]
    fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::sub_avx2_f64(a, b, result) },

            _ => self.swar.sub_f64(a, b, result),
        }
    }

    #[inline]
    fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::mul_avx2_f64(a, b, result) },

            _ => self.swar.mul_f64(a, b, result),
        }
    }

    #[inline]
    fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        self.swar.div_f64(a, b, result)
    }

    #[inline]
    fn fma_f64(&self, a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != c.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        self.swar.fma_f64(a, b, c, result)
    }

    #[inline]
    fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { x86::scale_avx2_f64(input, scalar, result) },

            _ => self.swar.scale_f64(input, scalar, result),
        }
    }

    #[inline]
    fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        if a.len() != b.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => Ok(unsafe { x86::dot_avx2_f64(a, b) }),

            _ => self.swar.dot_f64(a, b),
        }
    }

    #[inline]
    fn sum_f32(&self, input: &[f32]) -> Result<f32> {
        // Sum is inherently sequential, use SWAR implementation
        self.swar.sum_f32(input)
    }

    #[inline]
    fn max_f32(&self, input: &[f32]) -> Result<f32> {
        // Max operation: Use SWAR (SIMD Within A Register) as portable fallback.
        // Future enhancement: SIMD intrinsics for supported architectures.
        self.swar.max_f32(input)
    }

    #[inline]
    fn add_u32(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()> {
        // Integer operations use SWAR implementation
        self.swar.add_u32(a, b, result)
    }
}
