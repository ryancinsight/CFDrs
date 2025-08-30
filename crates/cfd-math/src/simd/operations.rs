//! SIMD operations with architecture-specific implementations
//!
//! Single implementation per operation following SSOT principle.
//! No redundant dispatch layers or "improved" versions.

use super::{SimdCapability, SwarOps};
use crate::error::Result;

/// SIMD operation types
#[derive(Debug, Clone, Copy)]
pub enum SimdOperation {
    Add,
    Mul,
    Sub,
    Div,
    FusedMulAdd,
}

/// Trait for vectorized operations
pub trait VectorOps {
    /// Add two vectors element-wise
    fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Subtract two vectors element-wise
    fn sub(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Multiply two vectors element-wise
    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Divide two vectors element-wise
    fn div(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Fused multiply-add: a * b + c
    fn fma(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()>;

    /// Scale a vector by a scalar
    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()>;

    /// Compute dot product of two vectors
    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32>;

    /// Add for f64
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Subtract for f64
    fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Multiply for f64
    fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Divide for f64
    fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Fused multiply-add for f64
    fn fma_f64(&self, a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()>;

    /// Scale for f64
    fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()>;

    /// Dot product for f64
    fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64>;
}

/// SIMD operations with automatic architecture dispatch
pub struct SimdOps {
    capability: SimdCapability,
    swar: SwarOps,
}

impl SimdOps {
    /// Create a new SIMD operations handler
    pub fn new() -> Self {
        Self {
            capability: SimdCapability::detect(),
            swar: SwarOps::new(),
        }
    }

    /// Get current capability
    pub fn capability(&self) -> SimdCapability {
        self.capability
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
            SimdCapability::Avx2 => unsafe { self.add_avx2_f32(a, b, result) },

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => unsafe { self.add_sse42_f32(a, b, result) },

            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            SimdCapability::Neon => unsafe { self.add_neon_f32(a, b, result) },

            _ => self.swar.add_f32(a, b, result),
        }
    }

    #[inline]
    fn sub(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { self.sub_avx2_f32(a, b, result) },

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => unsafe { self.sub_sse42_f32(a, b, result) },

            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            SimdCapability::Neon => unsafe { self.sub_neon_f32(a, b, result) },

            _ => self.swar.sub_f32(a, b, result),
        }
    }

    #[inline]
    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { self.mul_avx2_f32(a, b, result) },

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => unsafe { self.mul_sse42_f32(a, b, result) },

            #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
            SimdCapability::Neon => unsafe { self.mul_neon_f32(a, b, result) },

            _ => self.swar.mul_f32(a, b, result),
        }
    }

    #[inline]
    fn div(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Division is typically not SIMD-optimized, use SWAR
        self.swar.div_f32(a, b, result)
    }

    #[inline]
    fn fma(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != c.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "fma"))]
            SimdCapability::Avx2 => unsafe { self.fma_avx2_f32(a, b, c, result) },

            _ => self.swar.fma_f32(a, b, c, result),
        }
    }

    #[inline]
    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { self.scale_avx2_f32(input, scalar, result) },

            _ => self.swar.scale_f32(input, scalar, result),
        }
    }

    #[inline]
    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        if a.len() != b.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => Ok(unsafe { self.dot_avx2_f32(a, b) }),

            _ => self.swar.dot_f32(a, b),
        }
    }

    // f64 implementations
    #[inline]
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => unsafe { self.add_avx2_f64(a, b, result) },

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
            SimdCapability::Avx2 => unsafe { self.sub_avx2_f64(a, b, result) },

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
            SimdCapability::Avx2 => unsafe { self.mul_avx2_f64(a, b, result) },

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
            SimdCapability::Avx2 => unsafe { self.scale_avx2_f64(input, scalar, result) },

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
            SimdCapability::Avx2 => Ok(unsafe { self.dot_avx2_f64(a, b) }),

            _ => self.swar.dot_f64(a, b),
        }
    }
}

// Architecture-specific implementations
impl SimdOps {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn add_avx2_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !7; // Process 8 elements at a time

        for i in (0..simd_len).step_by(8) {
            let va = _mm256_loadu_ps(a.as_ptr().add(i));
            let vb = _mm256_loadu_ps(b.as_ptr().add(i));
            let vr = _mm256_add_ps(va, vb);
            _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
        }

        // Handle remainder
        for i in simd_len..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn sub_avx2_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !7;

        for i in (0..simd_len).step_by(8) {
            let va = _mm256_loadu_ps(a.as_ptr().add(i));
            let vb = _mm256_loadu_ps(b.as_ptr().add(i));
            let vr = _mm256_sub_ps(va, vb);
            _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = a[i] - b[i];
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn mul_avx2_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !7;

        for i in (0..simd_len).step_by(8) {
            let va = _mm256_loadu_ps(a.as_ptr().add(i));
            let vb = _mm256_loadu_ps(b.as_ptr().add(i));
            let vr = _mm256_mul_ps(va, vb);
            _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn scale_avx2_f32(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = input.len();
        let simd_len = len & !7;
        let vs = _mm256_set1_ps(scalar);

        for i in (0..simd_len).step_by(8) {
            let vi = _mm256_loadu_ps(input.as_ptr().add(i));
            let vr = _mm256_mul_ps(vi, vs);
            _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = input[i] * scalar;
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn dot_avx2_f32(&self, a: &[f32], b: &[f32]) -> f32 {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !7;
        let mut sum = _mm256_setzero_ps();

        for i in (0..simd_len).step_by(8) {
            let va = _mm256_loadu_ps(a.as_ptr().add(i));
            let vb = _mm256_loadu_ps(b.as_ptr().add(i));
            let prod = _mm256_mul_ps(va, vb);
            sum = _mm256_add_ps(sum, prod);
        }

        // Horizontal sum
        let sum_array = std::mem::transmute::<__m256, [f32; 8]>(sum);
        let mut result = sum_array.iter().sum::<f32>();

        // Handle remainder
        for i in simd_len..len {
            result += a[i] * b[i];
        }

        result
    }

    // Similar implementations for f64
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn add_avx2_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !3; // Process 4 elements at a time

        for i in (0..simd_len).step_by(4) {
            let va = _mm256_loadu_pd(a.as_ptr().add(i));
            let vb = _mm256_loadu_pd(b.as_ptr().add(i));
            let vr = _mm256_add_pd(va, vb);
            _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn sub_avx2_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !3;

        for i in (0..simd_len).step_by(4) {
            let va = _mm256_loadu_pd(a.as_ptr().add(i));
            let vb = _mm256_loadu_pd(b.as_ptr().add(i));
            let vr = _mm256_sub_pd(va, vb);
            _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = a[i] - b[i];
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn mul_avx2_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !3;

        for i in (0..simd_len).step_by(4) {
            let va = _mm256_loadu_pd(a.as_ptr().add(i));
            let vb = _mm256_loadu_pd(b.as_ptr().add(i));
            let vr = _mm256_mul_pd(va, vb);
            _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn scale_avx2_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        use std::arch::x86_64::*;

        let len = input.len();
        let simd_len = len & !3;
        let vs = _mm256_set1_pd(scalar);

        for i in (0..simd_len).step_by(4) {
            let vi = _mm256_loadu_pd(input.as_ptr().add(i));
            let vr = _mm256_mul_pd(vi, vs);
            _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
        }

        for i in simd_len..len {
            result[i] = input[i] * scalar;
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    unsafe fn dot_avx2_f64(&self, a: &[f64], b: &[f64]) -> f64 {
        use std::arch::x86_64::*;

        let len = a.len();
        let simd_len = len & !3;
        let mut sum = _mm256_setzero_pd();

        for i in (0..simd_len).step_by(4) {
            let va = _mm256_loadu_pd(a.as_ptr().add(i));
            let vb = _mm256_loadu_pd(b.as_ptr().add(i));
            let prod = _mm256_mul_pd(va, vb);
            sum = _mm256_add_pd(sum, prod);
        }

        // Horizontal sum
        let sum_array = std::mem::transmute::<__m256d, [f64; 4]>(sum);
        let mut result = sum_array.iter().sum::<f64>();

        // Handle remainder
        for i in simd_len..len {
            result += a[i] * b[i];
        }

        result
    }

    // SSE4.2 and NEON implementations would follow similar patterns
    // For brevity, defaulting to SWAR for non-AVX2 architectures
}

impl Default for SimdOps {
    fn default() -> Self {
        Self::new()
    }
}
