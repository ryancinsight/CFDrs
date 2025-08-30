//! Enhanced SWAR (SIMD Within A Register) operations
//!
//! Provides portable SIMD-like operations using standard arithmetic
//! for platforms without hardware SIMD support.

use crate::error::Result;

/// Enhanced SWAR operations with unrolling and optimization
pub struct SwarOps {
    unroll_factor: usize,
}

impl SwarOps {
    /// Create new SWAR operations handler
    pub fn new() -> Self {
        Self {
            unroll_factor: 4, // Default unrolling factor
        }
    }

    /// Add f32 arrays with loop unrolling
    #[inline]
    pub fn add_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        // Process 4 elements at a time for better instruction-level parallelism
        let chunks = len / 4;
        let remainder = len % 4;

        // Unrolled loop for main portion
        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] + b[idx];
            result[idx + 1] = a[idx + 1] + b[idx + 1];
            result[idx + 2] = a[idx + 2] + b[idx + 2];
            result[idx + 3] = a[idx + 3] + b[idx + 3];
        }

        // Handle remainder
        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    /// Multiply f32 arrays with loop unrolling
    #[inline]
    pub fn mul_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] * b[idx];
            result[idx + 1] = a[idx + 1] * b[idx + 1];
            result[idx + 2] = a[idx + 2] * b[idx + 2];
            result[idx + 3] = a[idx + 3] * b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    /// Subtract f32 arrays with loop unrolling
    #[inline]
    pub fn sub_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] - b[idx];
            result[idx + 1] = a[idx + 1] - b[idx + 1];
            result[idx + 2] = a[idx + 2] - b[idx + 2];
            result[idx + 3] = a[idx + 3] - b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] - b[i];
        }

        Ok(())
    }

    /// Divide f32 arrays with loop unrolling
    #[inline]
    pub fn div_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] / b[idx];
            result[idx + 1] = a[idx + 1] / b[idx + 1];
            result[idx + 2] = a[idx + 2] / b[idx + 2];
            result[idx + 3] = a[idx + 3] / b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] / b[i];
        }

        Ok(())
    }

    /// Fused multiply-add for f32
    #[inline]
    pub fn fma_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx].mul_add(b[idx], 0.0);
            result[idx + 1] = a[idx + 1].mul_add(b[idx + 1], 0.0);
            result[idx + 2] = a[idx + 2].mul_add(b[idx + 2], 0.0);
            result[idx + 3] = a[idx + 3].mul_add(b[idx + 3], 0.0);
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i].mul_add(b[i], 0.0);
        }

        Ok(())
    }

    /// Add f64 arrays with loop unrolling
    #[inline]
    pub fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] + b[idx];
            result[idx + 1] = a[idx + 1] + b[idx + 1];
            result[idx + 2] = a[idx + 2] + b[idx + 2];
            result[idx + 3] = a[idx + 3] + b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    /// Multiply f64 arrays with loop unrolling
    #[inline]
    pub fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] * b[idx];
            result[idx + 1] = a[idx + 1] * b[idx + 1];
            result[idx + 2] = a[idx + 2] * b[idx + 2];
            result[idx + 3] = a[idx + 3] * b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    /// Subtract f64 arrays with loop unrolling
    #[inline]
    pub fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] - b[idx];
            result[idx + 1] = a[idx + 1] - b[idx + 1];
            result[idx + 2] = a[idx + 2] - b[idx + 2];
            result[idx + 3] = a[idx + 3] - b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] - b[i];
        }

        Ok(())
    }

    /// Divide f64 arrays with loop unrolling
    #[inline]
    pub fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx] / b[idx];
            result[idx + 1] = a[idx + 1] / b[idx + 1];
            result[idx + 2] = a[idx + 2] / b[idx + 2];
            result[idx + 3] = a[idx + 3] / b[idx + 3];
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i] / b[i];
        }

        Ok(())
    }

    /// Fused multiply-add for f64
    #[inline]
    pub fn fma_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        let len = a.len();
        if len != b.len() || len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = a[idx].mul_add(b[idx], 0.0);
            result[idx + 1] = a[idx + 1].mul_add(b[idx + 1], 0.0);
            result[idx + 2] = a[idx + 2].mul_add(b[idx + 2], 0.0);
            result[idx + 3] = a[idx + 3].mul_add(b[idx + 3], 0.0);
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = a[i].mul_add(b[i], 0.0);
        }

        Ok(())
    }

    /// Optimized dot product for f32 using unrolling
    #[inline]
    pub fn dot_f32(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        let len = a.len();
        if len != b.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;
        let mut sum = 0.0f32;
        let mut sum1 = 0.0f32;
        let mut sum2 = 0.0f32;
        let mut sum3 = 0.0f32;

        // Unrolled loop with 4 accumulators for better pipelining
        for i in 0..chunks {
            let idx = i * 4;
            sum += a[idx] * b[idx];
            sum1 += a[idx + 1] * b[idx + 1];
            sum2 += a[idx + 2] * b[idx + 2];
            sum3 += a[idx + 3] * b[idx + 3];
        }

        // Combine accumulators
        sum += sum1 + sum2 + sum3;

        // Handle remainder
        let start = chunks * 4;
        for i in start..len {
            sum += a[i] * b[i];
        }

        Ok(sum)
    }

    /// Optimized dot product for f64 using unrolling
    #[inline]
    pub fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        let len = a.len();
        if len != b.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;
        let mut sum = 0.0f64;
        let mut sum1 = 0.0f64;
        let mut sum2 = 0.0f64;
        let mut sum3 = 0.0f64;

        for i in 0..chunks {
            let idx = i * 4;
            sum += a[idx] * b[idx];
            sum1 += a[idx + 1] * b[idx + 1];
            sum2 += a[idx + 2] * b[idx + 2];
            sum3 += a[idx + 3] * b[idx + 3];
        }

        sum += sum1 + sum2 + sum3;

        let start = chunks * 4;
        for i in start..len {
            sum += a[i] * b[i];
        }

        Ok(sum)
    }

    /// Scale a vector by a scalar (f32)
    #[inline]
    pub fn scale_f32(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        let len = input.len();
        if len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = input[idx] * scalar;
            result[idx + 1] = input[idx + 1] * scalar;
            result[idx + 2] = input[idx + 2] * scalar;
            result[idx + 3] = input[idx + 3] * scalar;
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = input[i] * scalar;
        }

        Ok(())
    }

    /// Scale a vector by a scalar (f64)
    #[inline]
    pub fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        let len = input.len();
        if len != result.len() {
            return Err("Vector length mismatch".into());
        }

        let chunks = len / 4;

        for i in 0..chunks {
            let idx = i * 4;
            result[idx] = input[idx] * scalar;
            result[idx + 1] = input[idx + 1] * scalar;
            result[idx + 2] = input[idx + 2] * scalar;
            result[idx + 3] = input[idx + 3] * scalar;
        }

        let start = chunks * 4;
        for i in start..len {
            result[i] = input[i] * scalar;
        }

        Ok(())
    }
}

impl Default for SwarOps {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_swar_add() {
        let swar = SwarOps::new();
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = vec![2.0, 3.0, 4.0, 5.0, 6.0];
        let mut result = vec![0.0; 5];

        swar.add_f32(&a, &b, &mut result).unwrap();

        for i in 0..5 {
            assert_relative_eq!(result[i], a[i] + b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_swar_dot_product() {
        let swar = SwarOps::new();
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![2.0, 3.0, 4.0, 5.0];

        let dot = swar.dot_f32(&a, &b).unwrap();
        let expected = 1.0 * 2.0 + 2.0 * 3.0 + 3.0 * 4.0 + 4.0 * 5.0;

        assert_relative_eq!(dot, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_swar_scale() {
        let swar = SwarOps::new();
        let input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let scalar = 2.5;
        let mut result = vec![0.0; 5];

        swar.scale_f64(&input, scalar, &mut result).unwrap();

        for i in 0..5 {
            assert_relative_eq!(result[i], input[i] * scalar, epsilon = 1e-10);
        }
    }
}
