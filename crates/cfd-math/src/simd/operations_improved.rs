//! Improved SIMD operations with proper architecture support and SWAR fallback
//!
//! This module provides safe, portable SIMD operations that:
//! - Use architecture-specific SIMD when available
//! - Fall back to SWAR (SIMD Within A Register) for portability
//! - Maintain zero-cost abstractions
//! - Support both f32 and f64 operations

use super::{SimdCapability, SwarOps};
use crate::error::Result;
use std::ops::{Add, Mul};

/// SIMD operation types
#[derive(Debug, Clone, Copy)]
pub enum SimdOp {
    Add,
    Mul,
    Sub,
    Div,
    FusedMulAdd,
}

/// Improved SIMD operations with proper dispatch
pub struct SimdProcessor {
    capability: SimdCapability,
    swar: SwarOps,
}

impl SimdProcessor {
    /// Create a new SIMD processor
    pub fn new() -> Self {
        Self {
            capability: SimdCapability::detect(),
            swar: SwarOps::new(),
        }
    }

    /// Process f32 arrays with the specified operation
    #[inline]
    pub fn process_f32(&self, a: &[f32], b: &[f32], result: &mut [f32], op: SimdOp) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("Vector length mismatch".into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => self.process_avx2_f32(a, b, result, op),

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => self.process_sse42_f32(a, b, result, op),

            #[cfg(target_arch = "aarch64")]
            SimdCapability::Neon => self.process_neon_f32(a, b, result, op),

            _ => self.process_swar_f32(a, b, result, op),
        }
    }

    /// Process f64 arrays with the specified operation
    #[inline]
    pub fn process_f64(&self, a: &[f64], b: &[f64], result: &mut [f64], op: SimdOp) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("Vector length mismatch".into());
        }

        match self.capability {
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
            SimdCapability::Avx2 => self.process_avx2_f64(a, b, result, op),

            #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
            SimdCapability::Sse42 => self.process_sse42_f64(a, b, result, op),

            #[cfg(target_arch = "aarch64")]
            SimdCapability::Neon => self.process_neon_f64(a, b, result, op),

            _ => self.process_swar_f64(a, b, result, op),
        }
    }

    // AVX2 implementations for x86_64
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    fn process_avx2_f32(&self, a: &[f32], b: &[f32], result: &mut [f32], op: SimdOp) -> Result<()> {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !7; // Process 8 elements at a time

            for i in (0..simd_len).step_by(8) {
                let va = _mm256_loadu_ps(a.as_ptr().add(i));
                let vb = _mm256_loadu_ps(b.as_ptr().add(i));

                let vr = match op {
                    SimdOp::Add => _mm256_add_ps(va, vb),
                    SimdOp::Mul => _mm256_mul_ps(va, vb),
                    SimdOp::Sub => _mm256_sub_ps(va, vb),
                    SimdOp::Div => _mm256_div_ps(va, vb),
                    SimdOp::FusedMulAdd => {
                        // Requires FMA support
                        #[cfg(target_feature = "fma")]
                        {
                            _mm256_fmadd_ps(va, vb, _mm256_setzero_ps())
                        }
                        #[cfg(not(target_feature = "fma"))]
                        {
                            _mm256_mul_ps(va, vb)
                        }
                    }
                };

                _mm256_storeu_ps(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match op {
                    SimdOp::Add => a[i] + b[i],
                    SimdOp::Mul => a[i] * b[i],
                    SimdOp::Sub => a[i] - b[i],
                    SimdOp::Div => a[i] / b[i],
                    SimdOp::FusedMulAdd => a[i] * b[i],
                };
            }
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    #[inline]
    fn process_avx2_f64(&self, a: &[f64], b: &[f64], result: &mut [f64], op: SimdOp) -> Result<()> {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = _mm256_loadu_pd(a.as_ptr().add(i));
                let vb = _mm256_loadu_pd(b.as_ptr().add(i));

                let vr = match op {
                    SimdOp::Add => _mm256_add_pd(va, vb),
                    SimdOp::Mul => _mm256_mul_pd(va, vb),
                    SimdOp::Sub => _mm256_sub_pd(va, vb),
                    SimdOp::Div => _mm256_div_pd(va, vb),
                    SimdOp::FusedMulAdd => {
                        #[cfg(target_feature = "fma")]
                        {
                            _mm256_fmadd_pd(va, vb, _mm256_setzero_pd())
                        }
                        #[cfg(not(target_feature = "fma"))]
                        {
                            _mm256_mul_pd(va, vb)
                        }
                    }
                };

                _mm256_storeu_pd(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match op {
                    SimdOp::Add => a[i] + b[i],
                    SimdOp::Mul => a[i] * b[i],
                    SimdOp::Sub => a[i] - b[i],
                    SimdOp::Div => a[i] / b[i],
                    SimdOp::FusedMulAdd => a[i] * b[i],
                };
            }
        }

        Ok(())
    }

    // SSE4.2 implementations for x86_64
    #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
    #[inline]
    fn process_sse42_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        op: SimdOp,
    ) -> Result<()> {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = _mm_loadu_ps(a.as_ptr().add(i));
                let vb = _mm_loadu_ps(b.as_ptr().add(i));

                let vr = match op {
                    SimdOp::Add => _mm_add_ps(va, vb),
                    SimdOp::Mul => _mm_mul_ps(va, vb),
                    SimdOp::Sub => _mm_sub_ps(va, vb),
                    SimdOp::Div => _mm_div_ps(va, vb),
                    SimdOp::FusedMulAdd => _mm_mul_ps(va, vb),
                };

                _mm_storeu_ps(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match op {
                    SimdOp::Add => a[i] + b[i],
                    SimdOp::Mul => a[i] * b[i],
                    SimdOp::Sub => a[i] - b[i],
                    SimdOp::Div => a[i] / b[i],
                    SimdOp::FusedMulAdd => a[i] * b[i],
                };
            }
        }

        Ok(())
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "sse4.2"))]
    #[inline]
    fn process_sse42_f64(
        &self,
        a: &[f64],
        b: &[f64],
        result: &mut [f64],
        op: SimdOp,
    ) -> Result<()> {
        use std::arch::x86_64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !1; // Process 2 elements at a time

            for i in (0..simd_len).step_by(2) {
                let va = _mm_loadu_pd(a.as_ptr().add(i));
                let vb = _mm_loadu_pd(b.as_ptr().add(i));

                let vr = match op {
                    SimdOp::Add => _mm_add_pd(va, vb),
                    SimdOp::Mul => _mm_mul_pd(va, vb),
                    SimdOp::Sub => _mm_sub_pd(va, vb),
                    SimdOp::Div => _mm_div_pd(va, vb),
                    SimdOp::FusedMulAdd => _mm_mul_pd(va, vb),
                };

                _mm_storeu_pd(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match op {
                    SimdOp::Add => a[i] + b[i],
                    SimdOp::Mul => a[i] * b[i],
                    SimdOp::Sub => a[i] - b[i],
                    SimdOp::Div => a[i] / b[i],
                    SimdOp::FusedMulAdd => a[i] * b[i],
                };
            }
        }

        Ok(())
    }

    // NEON implementations for ARM
    #[cfg(target_arch = "aarch64")]
    #[inline]
    fn process_neon_f32(&self, a: &[f32], b: &[f32], result: &mut [f32], op: SimdOp) -> Result<()> {
        use std::arch::aarch64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !3; // Process 4 elements at a time

            for i in (0..simd_len).step_by(4) {
                let va = vld1q_f32(a.as_ptr().add(i));
                let vb = vld1q_f32(b.as_ptr().add(i));

                let vr = match op {
                    SimdOp::Add => vaddq_f32(va, vb),
                    SimdOp::Mul => vmulq_f32(va, vb),
                    SimdOp::Sub => vsubq_f32(va, vb),
                    SimdOp::Div => vdivq_f32(va, vb),
                    SimdOp::FusedMulAdd => vmlaq_f32(vdupq_n_f32(0.0), va, vb),
                };

                vst1q_f32(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match op {
                    SimdOp::Add => a[i] + b[i],
                    SimdOp::Mul => a[i] * b[i],
                    SimdOp::Sub => a[i] - b[i],
                    SimdOp::Div => a[i] / b[i],
                    SimdOp::FusedMulAdd => a[i] * b[i],
                };
            }
        }

        Ok(())
    }

    #[cfg(target_arch = "aarch64")]
    #[inline]
    fn process_neon_f64(&self, a: &[f64], b: &[f64], result: &mut [f64], op: SimdOp) -> Result<()> {
        use std::arch::aarch64::*;

        unsafe {
            let len = a.len();
            let simd_len = len & !1; // Process 2 elements at a time

            for i in (0..simd_len).step_by(2) {
                let va = vld1q_f64(a.as_ptr().add(i));
                let vb = vld1q_f64(b.as_ptr().add(i));

                let vr = match op {
                    SimdOp::Add => vaddq_f64(va, vb),
                    SimdOp::Mul => vmulq_f64(va, vb),
                    SimdOp::Sub => vsubq_f64(va, vb),
                    SimdOp::Div => vdivq_f64(va, vb),
                    SimdOp::FusedMulAdd => vmlaq_f64(vdupq_n_f64(0.0), va, vb),
                };

                vst1q_f64(result.as_mut_ptr().add(i), vr);
            }

            // Handle remaining elements
            for i in simd_len..len {
                result[i] = match op {
                    SimdOp::Add => a[i] + b[i],
                    SimdOp::Mul => a[i] * b[i],
                    SimdOp::Sub => a[i] - b[i],
                    SimdOp::Div => a[i] / b[i],
                    SimdOp::FusedMulAdd => a[i] * b[i],
                };
            }
        }

        Ok(())
    }

    // SWAR fallback for all architectures
    #[inline]
    fn process_swar_f32(&self, a: &[f32], b: &[f32], result: &mut [f32], op: SimdOp) -> Result<()> {
        match op {
            SimdOp::Add => self.swar.add_f32(a, b, result),
            SimdOp::Mul => self.swar.mul_f32(a, b, result),
            SimdOp::Sub => self.swar.sub_f32(a, b, result),
            SimdOp::Div => self.swar.div_f32(a, b, result),
            SimdOp::FusedMulAdd => self.swar.fma_f32(a, b, result),
        }
    }

    #[inline]
    fn process_swar_f64(&self, a: &[f64], b: &[f64], result: &mut [f64], op: SimdOp) -> Result<()> {
        match op {
            SimdOp::Add => self.swar.add_f64(a, b, result),
            SimdOp::Mul => self.swar.mul_f64(a, b, result),
            SimdOp::Sub => self.swar.sub_f64(a, b, result),
            SimdOp::Div => self.swar.div_f64(a, b, result),
            SimdOp::FusedMulAdd => self.swar.fma_f64(a, b, result),
        }
    }

    // Fallback implementations for non-supported architectures
    #[cfg(not(any(
        all(
            target_arch = "x86_64",
            any(target_feature = "avx2", target_feature = "sse4.2")
        ),
        target_arch = "aarch64"
    )))]
    #[inline]
    fn process_avx2_f32(&self, a: &[f32], b: &[f32], result: &mut [f32], op: SimdOp) -> Result<()> {
        self.process_swar_f32(a, b, result, op)
    }

    #[cfg(not(any(
        all(
            target_arch = "x86_64",
            any(target_feature = "avx2", target_feature = "sse4.2")
        ),
        target_arch = "aarch64"
    )))]
    #[inline]
    fn process_avx2_f64(&self, a: &[f64], b: &[f64], result: &mut [f64], op: SimdOp) -> Result<()> {
        self.process_swar_f64(a, b, result, op)
    }

    #[cfg(not(any(
        all(target_arch = "x86_64", target_feature = "sse4.2"),
        target_arch = "aarch64"
    )))]
    #[inline]
    fn process_sse42_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        op: SimdOp,
    ) -> Result<()> {
        self.process_swar_f32(a, b, result, op)
    }

    #[cfg(not(any(
        all(target_arch = "x86_64", target_feature = "sse4.2"),
        target_arch = "aarch64"
    )))]
    #[inline]
    fn process_sse42_f64(
        &self,
        a: &[f64],
        b: &[f64],
        result: &mut [f64],
        op: SimdOp,
    ) -> Result<()> {
        self.process_swar_f64(a, b, result, op)
    }

    #[cfg(not(target_arch = "aarch64"))]
    #[inline]
    fn process_neon_f32(&self, a: &[f32], b: &[f32], result: &mut [f32], op: SimdOp) -> Result<()> {
        self.process_swar_f32(a, b, result, op)
    }

    #[cfg(not(target_arch = "aarch64"))]
    #[inline]
    fn process_neon_f64(&self, a: &[f64], b: &[f64], result: &mut [f64], op: SimdOp) -> Result<()> {
        self.process_swar_f64(a, b, result, op)
    }
}

/// High-level SIMD operations interface
impl SimdProcessor {
    /// Vectorized matrix-vector multiplication
    pub fn matvec_f32(
        &self,
        matrix: &[f32],
        vector: &[f32],
        result: &mut [f32],
        m: usize,
        n: usize,
    ) -> Result<()> {
        if matrix.len() != m * n || vector.len() != n || result.len() != m {
            return Err("Dimension mismatch".into());
        }

        for i in 0..m {
            let row_start = i * n;
            let row = &matrix[row_start..row_start + n];

            // Compute dot product using SIMD
            result[i] = self.dot_product_f32(row, vector)?;
        }

        Ok(())
    }

    /// Compute dot product with SIMD
    pub fn dot_product_f32(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        if a.len() != b.len() {
            return Err("Vector length mismatch".into());
        }

        let mut temp = vec![0.0f32; a.len()];
        self.process_f32(a, b, &mut temp, SimdOp::Mul)?;

        // Sum reduction (could be optimized with horizontal SIMD ops)
        Ok(temp.iter().sum())
    }

    /// Compute dot product with SIMD for f64
    pub fn dot_product_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        if a.len() != b.len() {
            return Err("Vector length mismatch".into());
        }

        let mut temp = vec![0.0f64; a.len()];
        self.process_f64(a, b, &mut temp, SimdOp::Mul)?;

        Ok(temp.iter().sum())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_simd_add_f32() {
        let processor = SimdProcessor::new();
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let b = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let mut result = vec![0.0; 8];

        processor
            .process_f32(&a, &b, &mut result, SimdOp::Add)
            .unwrap();

        for i in 0..8 {
            assert_relative_eq!(result[i], a[i] + b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_simd_mul_f64() {
        let processor = SimdProcessor::new();
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![2.0, 3.0, 4.0, 5.0];
        let mut result = vec![0.0; 4];

        processor
            .process_f64(&a, &b, &mut result, SimdOp::Mul)
            .unwrap();

        for i in 0..4 {
            assert_relative_eq!(result[i], a[i] * b[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_dot_product() {
        let processor = SimdProcessor::new();
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![2.0, 3.0, 4.0, 5.0];

        let dot = processor.dot_product_f32(&a, &b).unwrap();
        let expected = 1.0 * 2.0 + 2.0 * 3.0 + 3.0 * 4.0 + 4.0 * 5.0;

        assert_relative_eq!(dot, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_matvec_multiplication() {
        let processor = SimdProcessor::new();
        let matrix = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let vector = vec![1.0, 2.0, 3.0];
        let mut result = vec![0.0; 2];

        processor
            .matvec_f32(&matrix, &vector, &mut result, 2, 3)
            .unwrap();

        assert_relative_eq!(result[0], 14.0, epsilon = 1e-6); // 1*1 + 2*2 + 3*3
        assert_relative_eq!(result[1], 32.0, epsilon = 1e-6); // 4*1 + 5*2 + 6*3
    }
}
