//! Integer operations for SWAR

use super::SwarOps;
use crate::error::Result;

impl SwarOps {
    /// Add two u32 arrays element-wise
    #[inline]
    pub fn add_u32_arrays(&self, a: &[u32], b: &[u32], result: &mut [u32]) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        let len = a.len();
        let chunks = len / self.unroll_factor;

        for i in 0..chunks {
            let idx = i * self.unroll_factor;
            result[idx] = a[idx].wrapping_add(b[idx]);
            result[idx + 1] = a[idx + 1].wrapping_add(b[idx + 1]);
            result[idx + 2] = a[idx + 2].wrapping_add(b[idx + 2]);
            result[idx + 3] = a[idx + 3].wrapping_add(b[idx + 3]);
        }

        let start = chunks * self.unroll_factor;
        for i in start..len {
            result[i] = a[i].wrapping_add(b[i]);
        }

        Ok(())
    }
}
