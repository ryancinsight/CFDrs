//! SIMD optimizations for numerical operations
//!
//! Architecture-conditional SIMD using safe abstractions

mod arch_detect;
mod operations;
mod operations_dispatch;
mod operations_improved;
mod swar;
mod swar_enhanced;

#[cfg(test)]
mod tests;

/// SIMD capability levels
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdCapability {
    /// AVX2 (256-bit vectors)
    Avx2,
    /// SSE4.2 (128-bit vectors)
    Sse42,
    /// ARM NEON (128-bit vectors)
    Neon,
    /// Software SIMD within a register (fallback)
    Swar,
}

impl SimdCapability {
    /// Detect the best available SIMD capability
    pub fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx2") {
                return Self::Avx2;
            }
            if std::arch::is_x86_feature_detected!("sse4.2") {
                return Self::Sse42;
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            if std::arch::is_aarch64_feature_detected!("neon") {
                return Self::Neon;
            }
        }

        Self::Swar
    }
}

pub use arch_detect::ArchDetect;
pub use operations::{SimdOps, VectorOps};
pub use operations_improved::{SimdOp, SimdProcessor};
pub use swar::SwarOps;
pub use swar_enhanced::SwarOps as EnhancedSwarOps;

use std::arch::is_x86_feature_detected;

/// Dot product with architecture-specific SIMD
#[inline]
pub fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len(), "Vectors must have same length");

    // Use SIMD on x86_64 with AVX2
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { dot_product_avx2(a, b) };
        }
    }

    // Fallback to SWAR (SIMD Within A Register) for portability
    dot_product_swar(a, b)
}

/// SWAR implementation for portable SIMD-like optimization
#[inline]
fn dot_product_swar(a: &[f64], b: &[f64]) -> f64 {
    let chunks = a.len() / 4;
    let remainder = a.len() % 4;

    let mut sum = 0.0;

    // Process 4 elements at a time
    for i in 0..chunks {
        let idx = i * 4;
        let a0 = a[idx];
        let a1 = a[idx + 1];
        let a2 = a[idx + 2];
        let a3 = a[idx + 3];

        let b0 = b[idx];
        let b1 = b[idx + 1];
        let b2 = b[idx + 2];
        let b3 = b[idx + 3];

        // Compute 4 products in parallel (compiler can vectorize)
        sum += a0 * b0 + a1 * b1 + a2 * b2 + a3 * b3;
    }

    // Handle remainder
    let start = chunks * 4;
    for i in 0..remainder {
        sum += a[start + i] * b[start + i];
    }

    sum
}

/// AVX2 implementation for x86_64
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dot_product_avx2(a: &[f64], b: &[f64]) -> f64 {
    use std::arch::x86_64::*;

    let mut sum = _mm256_setzero_pd();
    let chunks = a.len() / 4;

    for i in 0..chunks {
        let idx = i * 4;
        let a_vec = _mm256_loadu_pd(a.as_ptr().add(idx));
        let b_vec = _mm256_loadu_pd(b.as_ptr().add(idx));
        let prod = _mm256_mul_pd(a_vec, b_vec);
        sum = _mm256_add_pd(sum, prod);
    }

    // Horizontal sum
    let sum_array = std::mem::transmute::<__m256d, [f64; 4]>(sum);
    let mut result = sum_array[0] + sum_array[1] + sum_array[2] + sum_array[3];

    // Handle remainder
    let start = chunks * 4;
    for i in start..a.len() {
        result += a[i] * b[i];
    }

    result
}

/// Element-wise vector addition with SIMD
#[inline]
pub fn vector_add(a: &[f64], b: &[f64], result: &mut [f64]) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), result.len());

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { vector_add_avx2(a, b, result) };
            return;
        }
    }

    // SWAR fallback
    vector_add_swar(a, b, result);
}

/// SWAR vector addition
#[inline]
fn vector_add_swar(a: &[f64], b: &[f64], result: &mut [f64]) {
    let chunks = a.len() / 4;

    // Process 4 elements at a time
    for i in 0..chunks {
        let idx = i * 4;
        result[idx] = a[idx] + b[idx];
        result[idx + 1] = a[idx + 1] + b[idx + 1];
        result[idx + 2] = a[idx + 2] + b[idx + 2];
        result[idx + 3] = a[idx + 3] + b[idx + 3];
    }

    // Handle remainder
    let start = chunks * 4;
    for i in start..a.len() {
        result[i] = a[i] + b[i];
    }
}

/// AVX2 vector addition
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn vector_add_avx2(a: &[f64], b: &[f64], result: &mut [f64]) {
    use std::arch::x86_64::*;

    let chunks = a.len() / 4;

    for i in 0..chunks {
        let idx = i * 4;
        let a_vec = _mm256_loadu_pd(a.as_ptr().add(idx));
        let b_vec = _mm256_loadu_pd(b.as_ptr().add(idx));
        let sum = _mm256_add_pd(a_vec, b_vec);
        _mm256_storeu_pd(result.as_mut_ptr().add(idx), sum);
    }

    // Handle remainder
    let start = chunks * 4;
    for i in start..a.len() {
        result[i] = a[i] + b[i];
    }
}

// Tests are in the separate tests.rs file
