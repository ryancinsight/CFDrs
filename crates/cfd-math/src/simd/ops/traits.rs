//! SIMD operation traits and types

use crate::error::Result;

/// SIMD operation types
#[derive(Debug, Clone, Copy)]
pub enum SimdOperation {
    /// Vector addition
    Add,
    /// Vector multiplication
    Mul,
    /// Vector subtraction
    Sub,
    /// Vector division
    Div,
    /// Fused multiply-add
    FusedMulAdd,
}

/// Trait for vectorized operations
pub trait VectorOps {
    /// Add two vectors element-wise
    fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Subtract two vectors element-wise
    fn sub(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Multiply two vectors element-wise
    fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Divide two vectors element-wise
    fn div(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()>;

    /// Fused multiply-add: a * b + c
    fn fma(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()>;

    /// Scale a vector by a scalar
    fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()>;

    /// Compute dot product of two vectors
    fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32>;

    /// Add for f64
    fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Subtract for f64
    fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Multiply for f64
    fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Divide for f64
    fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()>;

    /// Fused multiply-add for f64
    fn fma_f64(&self, a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()>;

    /// Scale for f64
    fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()>;

    /// Dot product for f64
    fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64>;
    
    /// Sum all elements in f32 array
    fn sum_f32(&self, input: &[f32]) -> Result<f32>;
    
    /// Find maximum element in f32 array
    fn max_f32(&self, input: &[f32]) -> Result<f32>;
    
    /// Add two u32 arrays element-wise
    fn add_u32(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()>;
}
