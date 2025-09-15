//! Reduction operations for SWAR

use super::SwarOps;
use crate::error::Result;

impl SwarOps {
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
    pub fn sum_f32_array(&self, input: &[f32]) -> Result<f32> {
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

        Ok(sum)
    }

    /// Find maximum element in f32 array
    #[inline]
    pub fn max_f32_array(&self, input: &[f32]) -> Result<f32> {
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
}
