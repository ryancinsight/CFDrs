//! SWAR (SIMD Within A Register) operations
//!
//! Portable SIMD-like operations using standard arithmetic
//! for platforms without hardware SIMD support.

use crate::error::Result;

const UNROLL_FACTOR: usize = 4;

/// SWAR operations for portable vectorization
pub struct SwarOps {
    unroll_factor: usize,
}

impl SwarOps {
    /// Create new SWAR operations handler
    pub fn new() -> Self {
        Self {
            unroll_factor: UNROLL_FACTOR,
        }
    }

    /// Add f32 arrays with loop unrolling
    #[inline]
    pub fn add_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;
        let remainder = len % self.unroll_factor;

        // Unrolled loop for main portion
        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx] + b[idx];
            result[idx + 1] = a[idx + 1] + b[idx + 1];
            result[idx + 2] = a[idx + 2] + b[idx + 2];
            result[idx + 3] = a[idx + 3] + b[idx + 3];
        }

        // Handle remainder
        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    /// Subtract f32 arrays with loop unrolling
    #[inline]
    pub fn sub_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx] - b[idx];
            result[idx + 1] = a[idx + 1] - b[idx + 1];
            result[idx + 2] = a[idx + 2] - b[idx + 2];
            result[idx + 3] = a[idx + 3] - b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i] - b[i];
        }

        Ok(())
    }

    /// Multiply f32 arrays with loop unrolling
    #[inline]
    pub fn mul_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx] * b[idx];
            result[idx + 1] = a[idx + 1] * b[idx + 1];
            result[idx + 2] = a[idx + 2] * b[idx + 2];
            result[idx + 3] = a[idx + 3] * b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    /// Divide f32 arrays
    #[inline]
    pub fn div_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        for i in 0..a.len() {
            result[i] = a[i] / b[i];
        }
        Ok(())
    }

    /// Fused multiply-add for f32
    #[inline]
    pub fn fma_f32(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx].mul_add(b[idx], c[idx]);
            result[idx + 1] = a[idx + 1].mul_add(b[idx + 1], c[idx + 1]);
            result[idx + 2] = a[idx + 2].mul_add(b[idx + 2], c[idx + 2]);
            result[idx + 3] = a[idx + 3].mul_add(b[idx + 3], c[idx + 3]);
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i].mul_add(b[i], c[i]);
        }

        Ok(())
    }

    /// Scale f32 array by scalar
    #[inline]
    pub fn scale_f32(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        let len = input.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = input[idx] * scalar;
            result[idx + 1] = input[idx + 1] * scalar;
            result[idx + 2] = input[idx + 2] * scalar;
            result[idx + 3] = input[idx + 3] * scalar;
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = input[i] * scalar;
        }

        Ok(())
    }

    /// Dot product for f32
    #[inline]
    pub fn dot_f32(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        let len = a.len();
        let chunks = len / self.unroll_factor;
        let mut sum = 0.0;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            sum += a[idx] * b[idx]
                + a[idx + 1] * b[idx + 1]
                + a[idx + 2] * b[idx + 2]
                + a[idx + 3] * b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            sum += a[i] * b[i];
        }

        Ok(sum)
    }

    /// Add f64 arrays with loop unrolling
    #[inline]
    pub fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx] + b[idx];
            result[idx + 1] = a[idx + 1] + b[idx + 1];
            result[idx + 2] = a[idx + 2] + b[idx + 2];
            result[idx + 3] = a[idx + 3] + b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    /// Subtract f64 arrays with loop unrolling
    #[inline]
    pub fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx] - b[idx];
            result[idx + 1] = a[idx + 1] - b[idx + 1];
            result[idx + 2] = a[idx + 2] - b[idx + 2];
            result[idx + 3] = a[idx + 3] - b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i] - b[i];
        }

        Ok(())
    }

    /// Multiply f64 arrays with loop unrolling
    #[inline]
    pub fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx] * b[idx];
            result[idx + 1] = a[idx + 1] * b[idx + 1];
            result[idx + 2] = a[idx + 2] * b[idx + 2];
            result[idx + 3] = a[idx + 3] * b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    /// Divide f64 arrays
    #[inline]
    pub fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        for i in 0..a.len() {
            result[i] = a[i] / b[i];
        }
        Ok(())
    }

    /// Fused multiply-add for f64
    #[inline]
    pub fn fma_f64(&self, a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx].mul_add(b[idx], c[idx]);
            result[idx + 1] = a[idx + 1].mul_add(b[idx + 1], c[idx + 1]);
            result[idx + 2] = a[idx + 2].mul_add(b[idx + 2], c[idx + 2]);
            result[idx + 3] = a[idx + 3].mul_add(b[idx + 3], c[idx + 3]);
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i].mul_add(b[i], c[i]);
        }

        Ok(())
    }

    /// Scale f64 array by scalar
    #[inline]
    pub fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        let len = input.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = input[idx] * scalar;
            result[idx + 1] = input[idx + 1] * scalar;
            result[idx + 2] = input[idx + 2] * scalar;
            result[idx + 3] = input[idx + 3] * scalar;
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = input[i] * scalar;
        }

        Ok(())
    }

    /// Dot product for f64
    #[inline]
    pub fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        let len = a.len();
        let chunks = len / self.unroll_factor;
        let mut sum = 0.0;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            sum += a[idx] * b[idx]
                + a[idx + 1] * b[idx + 1]
                + a[idx + 2] * b[idx + 2]
                + a[idx + 3] * b[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            sum += a[i] * b[i];
        }

        Ok(sum)
    }

    /// Process binary operation on f32 arrays
    pub fn process_binary_f32<F>(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        op: F,
    ) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        op(a, b, result);
        Ok(())
    }

    /// Sum all elements in f32 array
    #[inline]
    pub fn sum_f32(&self, input: &[f32]) -> f32 {
        let len = input.len();
        let chunks = len / self.unroll_factor;
        let mut sum = 0.0f32;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            sum += input[idx] + input[idx + 1] + input[idx + 2] + input[idx + 3];
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            sum += input[i];
        }

        sum
    }

    /// Find maximum element in f32 array
    #[inline]
    pub fn max_f32(&self, input: &[f32]) -> Result<f32> {
        if input.is_empty() {
            return Err(crate::error::MathError::InvalidInput("Empty array".to_string()).into());
        }

        let len = input.len();
        let chunks = len / self.unroll_factor;
        let mut max = input[0];

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            max = max
                .max(input[idx])
                .max(input[idx + 1])
                .max(input[idx + 2])
                .max(input[idx + 3]);
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            max = max.max(input[i]);
        }

        Ok(max)
    }

    /// Add two u32 arrays element-wise
    #[inline]
    pub fn add_u32(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx].wrapping_add(b[idx]);
            result[idx + 1] = a[idx + 1].wrapping_add(b[idx + 1]);
            result[idx + 2] = a[idx + 2].wrapping_add(b[idx + 2]);
            result[idx + 3] = a[idx + 3].wrapping_add(b[idx + 3]);
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i].wrapping_add(b[i]);
        }

        Ok(())
    }
}

impl Default for SwarOps {
    fn default() -> Self {
        Self::new()
    }
}
