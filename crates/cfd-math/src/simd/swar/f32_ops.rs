//! f32 SWAR operations

use super::SwarOps;
use crate::error::Result;

impl SwarOps {
    /// Add f32 arrays with loop unrolling
    #[inline]
    pub fn add_f32_arrays(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        let len = a.len();
        let chunks = len / self.unroll_factor;

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
    pub fn sub_f32_arrays(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
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
    pub fn mul_f32_arrays(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
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
    pub fn div_f32_arrays(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        for i in 0..a.len() {
            result[i] = a[i] / b[i];
        }
        Ok(())
    }

    /// Fused multiply-add for f32
    #[inline]
    pub fn fma_f32_arrays(
        &self,
        a: &[f32],
        b: &[f32],
        c: &[f32],
        result: &mut [f32],
    ) -> Result<()> {
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
    pub fn scale_f32_array(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
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
    pub fn dot_f32_arrays(&self, a: &[f32], b: &[f32]) -> Result<f32> {
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
}
