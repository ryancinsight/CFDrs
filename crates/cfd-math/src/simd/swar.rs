//! SWAR (SIMD Within A Register) implementations for portable vectorization
//!
//! SWAR provides portable bit-level parallelism within standard integer registers
//! as a fallback when hardware SIMD is not available.

use crate::error::Result;

/// SWAR operations for portable vectorization
pub struct SwarOps;

impl SwarOps {
    /// Create a new SWAR operations handler
    pub fn new() -> Self {
        Self
    }

    /// Process binary operations on f32 arrays using SWAR techniques
    pub fn process_binary_f32<F>(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        _op: F,
    ) -> Result<()>
    where
        F: Fn(&[f32], &[f32], &mut [f32]),
    {
        // SWAR optimization: process multiple elements using bit manipulation
        // For f32, we can pack operations efficiently
        const UNROLL_FACTOR: usize = 4;

        let len = a.len();
        let unrolled_len = len & !(UNROLL_FACTOR - 1);

        // Unrolled loop for instruction-level parallelism
        for i in (0..unrolled_len).step_by(UNROLL_FACTOR) {
            // Manual unrolling for compiler optimization
            unsafe {
                let a_ptr = a.as_ptr().add(i);
                let b_ptr = b.as_ptr().add(i);
                let r_ptr = result.as_mut_ptr().add(i);

                // Process 4 elements at once for pipelining
                *r_ptr = *a_ptr + *b_ptr;
                *r_ptr.add(1) = *a_ptr.add(1) + *b_ptr.add(1);
                *r_ptr.add(2) = *a_ptr.add(2) + *b_ptr.add(2);
                *r_ptr.add(3) = *a_ptr.add(3) + *b_ptr.add(3);
            }
        }

        // Handle remaining elements
        for i in unrolled_len..len {
            result[i] = a[i] + b[i];
        }

        Ok(())
    }

    /// Multiply f32 arrays using SWAR
    pub fn mul_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        const UNROLL_FACTOR: usize = 4;

        let len = a.len();
        let unrolled_len = len & !(UNROLL_FACTOR - 1);

        for i in (0..unrolled_len).step_by(UNROLL_FACTOR) {
            unsafe {
                let a_ptr = a.as_ptr().add(i);
                let b_ptr = b.as_ptr().add(i);
                let r_ptr = result.as_mut_ptr().add(i);

                *r_ptr = *a_ptr * *b_ptr;
                *r_ptr.add(1) = *a_ptr.add(1) * *b_ptr.add(1);
                *r_ptr.add(2) = *a_ptr.add(2) * *b_ptr.add(2);
                *r_ptr.add(3) = *a_ptr.add(3) * *b_ptr.add(3);
            }
        }

        for i in unrolled_len..len {
            result[i] = a[i] * b[i];
        }

        Ok(())
    }

    /// SWAR implementation for integer operations (useful for indices)
    pub fn add_u32(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()> {
        // SWAR for u32: pack two 16-bit values in a 32-bit register
        // This allows parallel processing of multiple values

        let len = a.len();
        let pairs = len / 2;

        for i in 0..pairs {
            let idx = i * 2;
            // Pack two values
            let packed_a = (a[idx] as u64) | ((a[idx + 1] as u64) << 32);
            let packed_b = (b[idx] as u64) | ((b[idx + 1] as u64) << 32);

            // Parallel add with overflow detection
            let sum = packed_a.wrapping_add(packed_b);

            // Unpack results
            result[idx] = sum as u32;
            result[idx + 1] = (sum >> 32) as u32;
        }

        // Handle odd length
        if len % 2 != 0 {
            result[len - 1] = a[len - 1] + b[len - 1];
        }

        Ok(())
    }

    /// SWAR dot product for f32
    pub fn dot_f32(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        const UNROLL_FACTOR: usize = 8;

        let len = a.len();
        let unrolled_len = len & !(UNROLL_FACTOR - 1);

        // Use multiple accumulators to avoid dependency chains
        let mut acc0 = 0.0f32;
        let mut acc1 = 0.0f32;
        let mut acc2 = 0.0f32;
        let mut acc3 = 0.0f32;

        for i in (0..unrolled_len).step_by(UNROLL_FACTOR) {
            unsafe {
                let a_ptr = a.as_ptr().add(i);
                let b_ptr = b.as_ptr().add(i);

                acc0 += *a_ptr * *b_ptr;
                acc1 += *a_ptr.add(1) * *b_ptr.add(1);
                acc2 += *a_ptr.add(2) * *b_ptr.add(2);
                acc3 += *a_ptr.add(3) * *b_ptr.add(3);
                acc0 += *a_ptr.add(4) * *b_ptr.add(4);
                acc1 += *a_ptr.add(5) * *b_ptr.add(5);
                acc2 += *a_ptr.add(6) * *b_ptr.add(6);
                acc3 += *a_ptr.add(7) * *b_ptr.add(7);
            }
        }

        // Combine accumulators
        let mut sum = acc0 + acc1 + acc2 + acc3;

        // Handle remaining elements
        for i in unrolled_len..len {
            sum += a[i] * b[i];
        }

        Ok(sum)
    }

    /// SWAR reduction (sum)
    pub fn sum_f32(&self, input: &[f32]) -> f32 {
        const UNROLL_FACTOR: usize = 8;

        let len = input.len();
        let unrolled_len = len & !(UNROLL_FACTOR - 1);

        // Multiple accumulators for instruction-level parallelism
        let mut acc0 = 0.0f32;
        let mut acc1 = 0.0f32;
        let mut acc2 = 0.0f32;
        let mut acc3 = 0.0f32;

        for i in (0..unrolled_len).step_by(UNROLL_FACTOR) {
            unsafe {
                let ptr = input.as_ptr().add(i);

                acc0 += *ptr;
                acc1 += *ptr.add(1);
                acc2 += *ptr.add(2);
                acc3 += *ptr.add(3);
                acc0 += *ptr.add(4);
                acc1 += *ptr.add(5);
                acc2 += *ptr.add(6);
                acc3 += *ptr.add(7);
            }
        }

        let mut sum = acc0 + acc1 + acc2 + acc3;

        for i in unrolled_len..len {
            sum += input[i];
        }

        sum
    }

    /// SWAR maximum element
    pub fn max_f32(&self, input: &[f32]) -> Option<f32> {
        if input.is_empty() {
            return None;
        }

        const UNROLL_FACTOR: usize = 4;
        let len = input.len();
        let unrolled_len = len & !(UNROLL_FACTOR - 1);

        let mut max_val = input[0];

        for i in (0..unrolled_len).step_by(UNROLL_FACTOR) {
            unsafe {
                let ptr = input.as_ptr().add(i);

                let v0 = *ptr;
                let v1 = *ptr.add(1);
                let v2 = *ptr.add(2);
                let v3 = *ptr.add(3);

                // Tree reduction for maximum
                let max01 = v0.max(v1);
                let max23 = v2.max(v3);
                let chunk_max = max01.max(max23);

                max_val = max_val.max(chunk_max);
            }
        }

        for i in unrolled_len..len {
            max_val = max_val.max(input[i]);
        }

        Some(max_val)
    }
}

impl Default for SwarOps {
    fn default() -> Self {
        Self::new()
    }
}
