//! Architecture detection utilities for SIMD support

use super::SimdCapability;

/// Runtime architecture detection
pub struct ArchDetect {
    capability: SimdCapability,
}

impl ArchDetect {
    /// Create a new architecture detector
    pub fn new() -> Self {
        Self {
            capability: SimdCapability::detect(),
        }
    }

    /// Get the detected SIMD capability
    pub fn capability(&self) -> SimdCapability {
        self.capability
    }

    /// Check if AVX2 is available
    #[inline]
    pub fn has_avx2(&self) -> bool {
        self.capability == SimdCapability::Avx2
    }

    /// Check if SSE4.2 is available
    #[inline]
    pub fn has_sse42(&self) -> bool {
        matches!(
            self.capability,
            SimdCapability::Sse42 | SimdCapability::Avx2
        )
    }

    /// Check if NEON is available
    #[inline]
    pub fn has_neon(&self) -> bool {
        self.capability == SimdCapability::Neon
    }

    /// Get vector width in elements for f32
    pub fn vector_width_f32(&self) -> usize {
        match self.capability {
            SimdCapability::Avx2 => 8,  // 256-bit / 32-bit
            SimdCapability::Sse42 | SimdCapability::Neon => 4, // 128-bit / 32-bit
            SimdCapability::Swar => 1,  // Scalar fallback
        }
    }

    /// Get vector width in elements for f64
    pub fn vector_width_f64(&self) -> usize {
        match self.capability {
            SimdCapability::Avx2 => 4,  // 256-bit / 64-bit
            SimdCapability::Sse42 | SimdCapability::Neon => 2, // 128-bit / 64-bit
            SimdCapability::Swar => 1,  // Scalar fallback
        }
    }
}

impl Default for ArchDetect {
    fn default() -> Self {
        Self::new()
    }
}
