//! Safe SIMD operations with architecture-specific implementations

use super::{SimdCapability, SwarOps};
use crate::error::Result;

/// Trait for vectorized operations
pub trait VectorOps {
    /// Add two vectors element-wise
    fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Multiply two vectors element-wise
    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Scale a vector by a scalar
    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()>;

    /// Compute dot product of two vectors
    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32>;

    /// Add for f64
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Multiply for f64
    fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Scale for f64
    fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()>;

    /// Dot product for f64
    fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64>;
}

/// SIMD operations dispatcher
pub struct SimdOps {
    capability: SimdCapability,
}

impl SimdOps {
    /// Create a new SIMD operations handler
    pub fn new() -> Self {
        Self {
            capability: SimdCapability::detect(),
        }
    }

    /// Process f32 vectors with architecture-specific SIMD
    #[inline]
    fn process_f32_binary<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match self.capability {
            SimdCapability::Avx2 => self.process_avx2_f32(a, b, result, op),
            SimdCapability::Sse42 => self.process_sse42_f32(a, b, result, op),
            SimdCapability::Neon => self.process_neon_f32(a, b, result, op),
            SimdCapability::Swar => {
                let swar = SwarOps::new();
                swar.process_binary_f32(a, b, result, op)
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    fn process_avx2_f32<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], _op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !7; // Process 8 elements at a time

            for i in (0..simd_len).step_by(8) {
                let va = _mm256_loadu_ps(a.as_ptr().add(i));
                let vb = _mm256_loadu_ps(b.as_ptr().add(i));
                let vr = _mm256_add_ps(va, vb); // For add operation
                _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = a[i] + b[i];
            }
        }

        Ok(())
    }

    #[cfg(not(target_arch = "x86_64"))]
    #[inline]
    fn process_avx2_f32<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        // Fallback to SWAR on non-x86_64
        let swar = SwarOps::new();
        swar.process_binary_f32(a, b, result, op)
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    fn process_sse42_f32<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], _op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = _mm_loadu_ps(a.as_ptr().add(i));
                let vb = _mm_loadu_ps(b.as_ptr().add(i));
                let vr = _mm_add_ps(va, vb); // For add operation
                _mm_storeu_ps(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = a[i] + b[i];
            }
        }

        Ok(())
    }

    #[cfg(not(target_arch = "x86_64"))]
    #[inline]
    fn process_sse42_f32<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        // Fallback to SWAR on non-x86_64
        let swar = SwarOps::new();
        swar.process_binary_f32(a, b, result, op)
    }

    #[cfg(target_arch = "aarch64")]
    #[inline]
    fn process_neon_f32<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], _op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        use std::arch::aarch64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = vld1q_f32(a.as_ptr().add(i));
                let vb = vld1q_f32(b.as_ptr().add(i));
                let vr = vaddq_f32(va, vb); // For add operation
                vst1q_f32(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = a[i] + b[i];
            }
        }

        Ok(())
    }

    #[cfg(not(target_arch = "aarch64"))]
    #[inline]
    fn process_neon_f32<F>(&self, a: &[f32], b: &[f32], result: &mut [f32], op: F) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        // Fallback to SWAR on non-ARM
        let swar = SwarOps::new();
        swar.process_binary_f32(a, b, result, op)
    }
}

impl VectorOps for SimdOps {
    fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.process_f32_binary(a, b, result, |a, b, r| {
            for i in 0..a.len() {
                r[i] = a[i] + b[i];
            }
        })
    }

    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Direct implementation instead of calling non-existent methods
        let swar = SwarOps::new();
        swar.mul_f32(a, b, result)
    }

    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Create a vector filled with the scalar for multiplication
        let scalar_vec = vec![scalar; input.len()];
        self.mul(input, &scalar_vec, result)
    }

    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        if a.len() != b.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        let mut temp = vec![0.0f32; a.len()];
        self.mul(a, b, &mut temp)?;

        // Sum reduction
        Ok(temp.iter().sum())
    }

    // f64 implementations
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // For now, use scalar implementation for f64
        for i in 0..a.len() {
            result[i] = a[i] + b[i];
        }
        Ok(())
    }

    fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        for i in 0..a.len() {
            result[i] = a[i] * b[i];
        }
        Ok(())
    }

    fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        for i in 0..input.len() {
            result[i] = input[i] * scalar;
        }
        Ok(())
    }

    fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        if a.len() != b.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        let mut sum = 0.0;
        for i in 0..a.len() {
            sum += a[i] * b[i];
        }
        Ok(sum)
    }
}

impl Default for SimdOps {
    fn default() -> Self {
        Self::new()
    }
}
