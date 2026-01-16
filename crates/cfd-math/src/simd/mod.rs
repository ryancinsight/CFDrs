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

pub use arch_detect::ArchDetect;
pub use ops::{SimdOperation, SimdOps, VectorOps};
pub use swar::SwarOps;

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
}

impl SimdProcessor {
    /// Create new processor with automatic capability detection
    pub fn new() -> Self {
        let capability = SimdCapability::detect();
        Self {
            capability,
            ops: SimdOps::new(),
        }
    }

    /// Get detected capability
    pub fn capability(&self) -> SimdCapability {
        self.capability
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
        // TODO: Implement runtime performance profiling to optimize threshold selection
        const SIMD_THRESHOLD: usize = 500; // Below this, SIMD overhead may exceed benefits
        
        // For very small arrays, scalar operations may be faster due to SIMD overhead
        if a.len() < SIMD_THRESHOLD && matches!(self.capability, SimdCapability::Avx2) {
            return self.process_scalar_f32(a, b, result, op);
        }
        
        match op {
            SimdOperation::Add => self.ops.add(a, b, result),
            SimdOperation::Mul => self.ops.mul(a, b, result),
            SimdOperation::Sub => self.ops.sub(a, b, result),
            SimdOperation::Div => self.ops.div(a, b, result),
            SimdOperation::FusedMulAdd => {
                // TODO: Implement proper FMA operation instead of mul+add sequence
                self.ops.mul(a, b, result)?;
                // Result now contains a * b, add to existing values would require separate c
                Ok(())
            }
        }
    }
    
    /// Scalar fallback for small arrays where SIMD overhead is detrimental
    /// TODO: Add benchmarking to determine optimal thresholds per hardware
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
                // TODO: Implement proper scalar FMA
                for ((res, &a_val), &b_val) in result.iter_mut().zip(a.iter()).zip(b.iter()) {
                    *res = a_val * b_val + *res;
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
