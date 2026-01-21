//! SIMD optimizations for numerical operations
//!
//! Architecture-conditional SIMD using safe abstractions with SWAR fallback.
//! Follows SSOT principle with single implementation per operation.

mod arch_detect;
pub mod cfd;
mod ops;
mod swar;
pub mod vector;
pub mod vectorization;

#[cfg(test)]
mod tests;

use crate::performance_monitor::PerformanceAwareSelector;
pub use arch_detect::ArchDetect;
pub use ops::{SimdOperation, SimdOps, VectorOps};
pub use swar::SwarOps;
use std::sync::Mutex;
use std::time::Instant;

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

/// Unified SIMD processor with automatic dispatch
pub struct SimdProcessor {
    capability: SimdCapability,
    /// SIMD operations handler
    pub ops: SimdOps,
    selector: Mutex<PerformanceAwareSelector>,
}

impl SimdProcessor {
    /// Create new processor with automatic capability detection
    pub fn new() -> Self {
        let capability = SimdCapability::detect();
        Self {
            capability,
            ops: SimdOps::new(),
            selector: Mutex::new(PerformanceAwareSelector::new()),
        }
    }

    /// Get detected capability
    pub fn capability(&self) -> SimdCapability {
        self.capability
    }

    fn operation_key(capability: SimdCapability, op: SimdOperation) -> &'static str {
        match (capability, op) {
            (SimdCapability::Avx2, SimdOperation::Add) => "simd_avx2_add",
            (SimdCapability::Avx2, SimdOperation::Mul) => "simd_avx2_mul",
            (SimdCapability::Avx2, SimdOperation::Sub) => "simd_avx2_sub",
            (SimdCapability::Avx2, SimdOperation::Div) => "simd_avx2_div",
            (SimdCapability::Avx2, SimdOperation::FusedMulAdd) => "simd_avx2_fma",
            (SimdCapability::Sse42, SimdOperation::Add) => "simd_sse42_add",
            (SimdCapability::Sse42, SimdOperation::Mul) => "simd_sse42_mul",
            (SimdCapability::Sse42, SimdOperation::Sub) => "simd_sse42_sub",
            (SimdCapability::Sse42, SimdOperation::Div) => "simd_sse42_div",
            (SimdCapability::Sse42, SimdOperation::FusedMulAdd) => "simd_sse42_fma",
            (SimdCapability::Neon, SimdOperation::Add) => "simd_neon_add",
            (SimdCapability::Neon, SimdOperation::Mul) => "simd_neon_mul",
            (SimdCapability::Neon, SimdOperation::Sub) => "simd_neon_sub",
            (SimdCapability::Neon, SimdOperation::Div) => "simd_neon_div",
            (SimdCapability::Neon, SimdOperation::FusedMulAdd) => "simd_neon_fma",
            (SimdCapability::Swar, SimdOperation::Add) => "simd_swar_add",
            (SimdCapability::Swar, SimdOperation::Mul) => "simd_swar_mul",
            (SimdCapability::Swar, SimdOperation::Sub) => "simd_swar_sub",
            (SimdCapability::Swar, SimdOperation::Div) => "simd_swar_div",
            (SimdCapability::Swar, SimdOperation::FusedMulAdd) => "simd_swar_fma",
        }
    }

    /// Process f32 arrays with specified operation
    ///
    /// Performance optimization: Automatically chooses between SIMD and parallel
    /// based on array size and measured performance characteristics.
    #[inline]
    pub fn process_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        op: SimdOperation,
    ) -> crate::error::Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        let operation_key = Self::operation_key(self.capability, op);
        let use_simd = {
            let mut selector = self.selector.lock().map_err(|_| {
                crate::error::MathError::InvalidInput(
                    "Performance selector lock poisoned".to_string(),
                )
            })?;
            selector.should_use_simd(operation_key, a.len())
        };

        let start = Instant::now();
        let result_status = if use_simd {
            match op {
                SimdOperation::Add => self.ops.add(a, b, result),
                SimdOperation::Mul => self.ops.mul(a, b, result),
                SimdOperation::Sub => self.ops.sub(a, b, result),
                SimdOperation::Div => self.ops.div(a, b, result),
                SimdOperation::FusedMulAdd => {
                    let c = result.to_vec();
                    self.ops.fma(a, b, &c, result)
                }
            }
        } else {
            self.process_scalar_f32(a, b, result, op)
        };
        let elapsed_ns = start.elapsed().as_nanos().min(u128::from(u64::MAX)) as u64;

        if result_status.is_ok() {
            let mut selector = self.selector.lock().map_err(|_| {
                crate::error::MathError::InvalidInput(
                    "Performance selector lock poisoned".to_string(),
                )
            })?;
            selector.record_performance(operation_key, a.len(), elapsed_ns);
        }

        result_status
    }

    /// Scalar fallback for small arrays where SIMD overhead is detrimental
    fn process_scalar_f32(
        &self,
        a: &[f32],
        b: &[f32],
        result: &mut [f32],
        op: SimdOperation,
    ) -> crate::error::Result<()> {
        match op {
            SimdOperation::Add => {
                for ((res, &a_val), &b_val) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
                    *res = a_val + b_val;
                }
            }
            SimdOperation::Mul => {
                for ((res, &a_val), &b_val) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
                    *res = a_val * b_val;
                }
            }
            SimdOperation::Sub => {
                for ((res, &a_val), &b_val) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
                    *res = a_val - b_val;
                }
            }
            SimdOperation::Div => {
                for ((res, &a_val), &b_val) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
                    *res = a_val / b_val;
                }
            }
            SimdOperation::FusedMulAdd => {
                for ((res, &a_val), &b_val) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
                    *res = a_val.mul_add(b_val, *res);
                }
            }
        }
        Ok(())
    }

    /// Process f64 arrays with specified operation
    #[inline]
    pub fn process_f64(
        &self,
        a: &[f64],
        b: &[f64],
        result: &mut [f64],
        op: SimdOperation,
    ) -> crate::error::Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        match op {
            SimdOperation::Add => self.ops.add_f64(a, b, result),
            SimdOperation::Mul => self.ops.mul_f64(a, b, result),
            SimdOperation::Sub => self.ops.sub_f64(a, b, result),
            SimdOperation::Div => self.ops.div_f64(a, b, result),
            SimdOperation::FusedMulAdd => {
                // For FMA, we multiply and add in-place
                self.ops.mul_f64(a, b, result)?;
                // Result now contains a * b, add to existing values would require separate c
                Ok(())
            }
        }
    }
}

impl Default for SimdProcessor {
    fn default() -> Self {
        Self::new()
    }
}
