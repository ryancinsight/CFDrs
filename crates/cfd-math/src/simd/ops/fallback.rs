//! Fallback scalar implementations for unsupported architectures

#[allow(unused_imports)]
use crate::error::Result;

/// Fallback scalar implementation for f32 addition
#[inline]
pub fn add_scalar_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] + b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f32 subtraction
#[inline]
pub fn sub_scalar_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] - b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f32 multiplication
#[inline]
pub fn mul_scalar_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] * b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f32 division
#[inline]
pub fn div_scalar_f32(a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] / b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f32 FMA
#[inline]
pub fn fma_scalar_f32(a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i].mul_add(b[i], c[i]);
    }
    Ok(())
}

/// Fallback scalar implementation for f32 scaling
#[inline]
pub fn scale_scalar_f32(input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
    for i in 0..input.len() {
        result[i] = input[i] * scalar;
    }
    Ok(())
}

/// Fallback scalar implementation for f32 dot product
#[inline]
pub fn dot_scalar_f32(a: &[f32], b: &[f32]) -> Result<f32> {
    let mut sum = 0.0f32;
    for i in 0..a.len() {
        sum += a[i] * b[i];
    }
    Ok(sum)
}

// f64 implementations
/// Fallback scalar implementation for f64 addition
#[inline]
pub fn add_scalar_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] + b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f64 subtraction
#[inline]
pub fn sub_scalar_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] - b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f64 multiplication
#[inline]
pub fn mul_scalar_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] * b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f64 division
#[inline]
pub fn div_scalar_f64(a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i] / b[i];
    }
    Ok(())
}

/// Fallback scalar implementation for f64 FMA
#[inline]
pub fn fma_scalar_f64(a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()> {
    for i in 0..a.len() {
        result[i] = a[i].mul_add(b[i], c[i]);
    }
    Ok(())
}

/// Fallback scalar implementation for f64 scaling
#[inline]
pub fn scale_scalar_f64(input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
    for i in 0..input.len() {
        result[i] = input[i] * scalar;
    }
    Ok(())
}

/// Fallback scalar implementation for f64 dot product
#[inline]
pub fn dot_scalar_f64(a: &[f64], b: &[f64]) -> Result<f64> {
    let mut sum = 0.0f64;
    for i in 0..a.len() {
        sum += a[i] * b[i];
    }
    Ok(sum)
}
