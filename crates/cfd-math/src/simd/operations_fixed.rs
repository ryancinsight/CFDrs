//! SIMD operations with proper operation dispatch
//!
//! This module provides SIMD-accelerated operations that properly handle
//! different operations instead of hardcoding addition.

use super::SimdCapability;

type Result<T> = std::result::Result<T, &'static str>;

/// SIMD operation type
#[derive(Debug, Clone, Copy)]
pub enum SimdOperation {
    Add,
    Subtract,
    Multiply,
    Divide,
}

/// SIMD operations processor with proper operation dispatch
pub struct SimdProcessor {
    capability: SimdCapability,
}

impl SimdProcessor {
    /// Create a new SIMD processor
    pub fn new(capability: SimdCapability) -> Self {
        Self { capability }
    }

    /// Process binary operation on f32 slices
    pub fn process_binary_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("Slice length mismatch");
        }

        match self.capability {
            SimdCapability::Avx2 => self.process_avx2_f32(a, b, result, operation),
            SimdCapability::Sse42 => self.process_sse42_f32(a, b, result, operation),
            SimdCapability::Neon => self.process_neon_f32(a, b, result, operation),
            SimdCapability::Swar => self.process_scalar_f32(a, b, result, operation),
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    fn process_avx2_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !7; // Process 8 elements at a time

            for i in (0..simd_len).step_by(8) {
                let va = _mm256_loadu_ps(a.as_ptr().add(i));
                let vb = _mm256_loadu_ps(b.as_ptr().add(i));

                let vr = match operation {
                    SimdOperation::Add => _mm256_add_ps(va, vb),
                    SimdOperation::Subtract => _mm256_sub_ps(va, vb),
                    SimdOperation::Multiply => _mm256_mul_ps(va, vb),
                    SimdOperation::Divide => _mm256_div_ps(va, vb),
                };

                _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match operation {
                    SimdOperation::Add => a[i] + b[i],
                    SimdOperation::Subtract => a[i] - b[i],
                    SimdOperation::Multiply => a[i] * b[i],
                    SimdOperation::Divide => a[i] / b[i],
                };
            }
        }

        Ok(())
    }

    #[cfg(not(target_arch = "x86_64"))]
    #[inline]
    fn process_avx2_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        // Fallback to scalar on non-x86_64
        self.process_scalar_f32(a, b, result, operation)
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    fn process_sse42_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = _mm_loadu_ps(a.as_ptr().add(i));
                let vb = _mm_loadu_ps(b.as_ptr().add(i));

                let vr = match operation {
                    SimdOperation::Add => _mm_add_ps(va, vb),
                    SimdOperation::Subtract => _mm_sub_ps(va, vb),
                    SimdOperation::Multiply => _mm_mul_ps(va, vb),
                    SimdOperation::Divide => _mm_div_ps(va, vb),
                };

                _mm_storeu_ps(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match operation {
                    SimdOperation::Add => a[i] + b[i],
                    SimdOperation::Subtract => a[i] - b[i],
                    SimdOperation::Multiply => a[i] * b[i],
                    SimdOperation::Divide => a[i] / b[i],
                };
            }
        }

        Ok(())
    }

    #[cfg(not(target_arch = "x86_64"))]
    #[inline]
    fn process_sse42_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        self.process_scalar_f32(a, b, result, operation)
    }

    #[cfg(target_arch = "aarch64")]
    #[inline]
    fn process_neon_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        use std::arch::aarch64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = vld1q_f32(a.as_ptr().add(i));
                let vb = vld1q_f32(b.as_ptr().add(i));

                let vr = match operation {
                    SimdOperation::Add => vaddq_f32(va, vb),
                    SimdOperation::Subtract => vsubq_f32(va, vb),
                    SimdOperation::Multiply => vmulq_f32(va, vb),
                    SimdOperation::Divide => vdivq_f32(va, vb),
                };

                vst1q_f32(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match operation {
                    SimdOperation::Add => a[i] + b[i],
                    SimdOperation::Subtract => a[i] - b[i],
                    SimdOperation::Multiply => a[i] * b[i],
                    SimdOperation::Divide => a[i] / b[i],
                };
            }
        }

        Ok(())
    }

    #[cfg(not(target_arch = "aarch64"))]
    #[inline]
    fn process_neon_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        self.process_scalar_f32(a, b, result, operation)
    }

    /// Scalar fallback implementation
    #[inline]
    fn process_scalar_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        operation: SimdOperation,
    ) -> Result<()> {
        for i in 0..a.len() {
            result[i] = match operation {
                SimdOperation::Add => a[i] + b[i],
                SimdOperation::Subtract => a[i] - b[i],
                SimdOperation::Multiply => a[i] * b[i],
                SimdOperation::Divide => a[i] / b[i],
            };
        }
        Ok(())
    }
}
