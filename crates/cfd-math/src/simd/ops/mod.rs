//! SIMD operations backed by hermes-simd.
//!
//! Delegates all architecture detection and dispatch to hermes, which handles
//! AVX2/SSE4.2/NEON/Scalar fallback internally via its `SimdOps` sealed trait.
//! This eliminates the prior x86/arm/SWAR duplication.

use crate::error::Result;
use hermes_simd::dispatch::SimdOps as HermesOps;

/// SIMD operations dispatcher backed by hermes-simd.
///
/// Zero-sized type — all methods delegate to hermes' generic runtime-dispatched
/// SIMD kernels. No manual architecture detection or SWAR fallback needed.
pub struct SimdOps;

impl SimdOps {
    /// Create a new SIMD operations handler.
    #[inline]
    pub fn new() -> Self {
        Self
    }

    // ── f32 operations ──────────────────────────────────────────────────

    /// Element-wise addition: `result[i] = a[i] + b[i]`.
    #[inline]
    pub fn add(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f32::elementwise_add(a, b, result).map_err(simd_err)
    }

    /// Element-wise subtraction: `result[i] = a[i] - b[i]`.
    #[inline]
    pub fn sub(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f32::elementwise_sub(a, b, result).map_err(simd_err)
    }

    /// Element-wise multiplication: `result[i] = a[i] * b[i]`.
    #[inline]
    pub fn mul(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f32::elementwise_mul(a, b, result).map_err(simd_err)
    }

    /// Element-wise division: `result[i] = a[i] / b[i]`.
    #[inline]
    pub fn div(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f32::elementwise_div(a, b, result).map_err(simd_err)
    }

    /// Fused multiply-add: `result[i] = a[i] * b[i] + c[i]`.
    ///
    /// Uses hardware FMA via `mul_add` when available.
    #[inline]
    pub fn fma(&self, a: &[f32], b: &[f32], c: &[f32], result: &mut [f32]) -> Result<()> {
        validate_same_len(a.len(), b.len(), c.len())?;
        validate_same_len(a.len(), result.len(), result.len())?;
        for ((r, &ai), (&bi, &ci)) in result.iter_mut().zip(a.iter()).zip(b.iter().zip(c.iter())) {
            *r = ai.mul_add(bi, ci);
        }
        Ok(())
    }

    /// Scalar multiplication: `result[i] = input[i] * scalar`.
    #[inline]
    pub fn scale(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        validate_same_len(input.len(), result.len(), result.len())?;
        for (r, &v) in result.iter_mut().zip(input.iter()) {
            *r = v * scalar;
        }
        Ok(())
    }

    /// Dot product: `sum(a[i] * b[i])`.
    #[inline]
    pub fn dot(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        validate_len_match(a.len(), b.len())?;
        f32::dot(a, b).map_err(simd_err)
    }

    // ── f64 operations ──────────────────────────────────────────────────

    /// Element-wise addition (f64): `result[i] = a[i] + b[i]`.
    #[inline]
    pub fn add_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f64::elementwise_add(a, b, result).map_err(simd_err)
    }

    /// Element-wise subtraction (f64): `result[i] = a[i] - b[i]`.
    #[inline]
    pub fn sub_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f64::elementwise_sub(a, b, result).map_err(simd_err)
    }

    /// Element-wise multiplication (f64): `result[i] = a[i] * b[i]`.
    #[inline]
    pub fn mul_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f64::elementwise_mul(a, b, result).map_err(simd_err)
    }

    /// Element-wise division (f64): `result[i] = a[i] / b[i]`.
    #[inline]
    pub fn div_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        f64::elementwise_div(a, b, result).map_err(simd_err)
    }

    /// Fused multiply-add (f64): `result[i] = a[i] * b[i] + c[i]`.
    #[inline]
    pub fn fma_f64(&self, a: &[f64], b: &[f64], c: &[f64], result: &mut [f64]) -> Result<()> {
        validate_same_len(a.len(), b.len(), c.len())?;
        validate_same_len(a.len(), result.len(), result.len())?;
        for ((r, &ai), (&bi, &ci)) in result.iter_mut().zip(a.iter()).zip(b.iter().zip(c.iter())) {
            *r = ai.mul_add(bi, ci);
        }
        Ok(())
    }

    /// Scalar multiplication (f64): `result[i] = input[i] * scalar`.
    #[inline]
    pub fn scale_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        validate_same_len(input.len(), result.len(), result.len())?;
        for (r, &v) in result.iter_mut().zip(input.iter()) {
            *r = v * scalar;
        }
        Ok(())
    }

    /// Dot product (f64): `sum(a[i] * b[i])`.
    #[inline]
    pub fn dot_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        validate_len_match(a.len(), b.len())?;
        f64::dot(a, b).map_err(simd_err)
    }

    // ── Reductions ──────────────────────────────────────────────────────

    /// Sum of all elements (f32).
    #[inline]
    pub fn sum_f32(&self, input: &[f32]) -> Result<f32> {
        Ok(<f32 as HermesOps>::sum(input))
    }

    /// Maximum element (f32).
    #[inline]
    pub fn max_f32(&self, input: &[f32]) -> Result<f32> {
        Ok(<f32 as HermesOps>::max(input))
    }

    // ── Integer operations (scalar fallback) ────────────────────────────

    /// Element-wise u32 addition (scalar — hermes does not seal u32).
    #[inline]
    pub fn add_u32(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()> {
        validate_same_len(a.len(), b.len(), result.len())?;
        for ((r, &ai), &bi) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
            *r = ai.wrapping_add(bi);
        }
        Ok(())
    }
}

impl Default for SimdOps {
    fn default() -> Self {
        Self
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────

fn validate_same_len(a: usize, b: usize, c: usize) -> Result<()> {
    if a == b && a == c {
        Ok(())
    } else {
        Err(cfd_core::error::Error::InvalidInput(
            "Dimension mismatch".to_string(),
        ))
    }
}

fn validate_len_match(a: usize, b: usize) -> Result<()> {
    if a == b {
        Ok(())
    } else {
        Err(cfd_core::error::Error::InvalidInput(
            "Dimension mismatch".to_string(),
        ))
    }
}

fn simd_err(e: hermes_simd::SimdError) -> cfd_core::error::Error {
    cfd_core::error::Error::InvalidInput(format!("SIMD error: {e:?}"))
}
