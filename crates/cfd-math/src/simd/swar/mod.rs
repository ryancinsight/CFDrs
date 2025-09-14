//! SWAR (SIMD Within A Register) operations
//!
//! Portable SIMD-like operations using standard arithmetic
//! for platforms without hardware SIMD support.

mod f32_ops;
mod f64_ops;
mod integer_ops;
mod reductions;

use crate::error::Result;
use crate::simd::VectorOps;

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

    /// Get the unroll factor used for optimization
    pub fn unroll_factor(&self) -> usize {
        self.unroll_factor
    }
}

impl Default for SwarOps {
    fn default() -> Self {
        Self::new()
    }
}

impl VectorOps for SwarOps {
    #[inline]
    fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.add_f32_arrays(a, b, result)
    }

    #[inline]
    fn sub(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.sub_f32_arrays(a, b, result)
    }

    #[inline]
    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.mul_f32_arrays(a, b, result)
    }

    #[inline]
    fn div(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.div_f32_arrays(a, b, result)
    }

    #[inline]
    fn fma(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
        self.fma_f32_arrays(a, b, c, result)
    }

    #[inline]
    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        self.scale_f32_array(input, scalar, result)
    }

    #[inline]
    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        self.dot_f32_arrays(a, b)
    }

    #[inline]
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        self.add_f64_arrays(a, b, result)
    }

    #[inline]
    fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        self.sub_f64_arrays(a, b, result)
    }

    #[inline]
    fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        self.mul_f64_arrays(a, b, result)
    }

    #[inline]
    fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        self.div_f64_arrays(a, b, result)
    }

    #[inline]
    fn fma_f64(&self, a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()> {
        self.fma_f64_arrays(a, b, c, result)
    }

    #[inline]
    fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        self.scale_f64_array(input, scalar, result)
    }

    #[inline]
    fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        self.dot_f64_arrays(a, b)
    }

    #[inline]
    fn sum_f32(&self, input: &[f32]) -> Result<f32> {
        self.sum_f32_array(input)
    }

    #[inline]
    fn max_f32(&self, input: &[f32]) -> Result<f32> {
        self.max_f32_array(input)
    }

    #[inline]
    fn add_u32(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()> {
        self.add_u32_arrays(a, b, result)
    }
}