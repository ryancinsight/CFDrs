//! Compression strategies for checkpoints

use serde::{Deserialize, Serialize};

/// Compression strategy for checkpoints
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum CompressionStrategy {
    /// No compression
    None,
    /// Zstandard compression with level (1-22, default 3)
    Zstd(u8),
}

impl Default for CompressionStrategy {
    fn default() -> Self {
        Self::None
    }
}

impl CompressionStrategy {
    /// Get recommended compression for checkpoint size
    pub fn recommended_for_size(size_bytes: usize) -> Self {
        if size_bytes > 10_000_000 {
            // > 10 MB: use compression
            Self::Zstd(3)
        } else {
            // Small files: no compression
            Self::None
        }
    }

    /// Get compression ratio estimate
    pub fn estimated_ratio(&self) -> f64 {
        match self {
            Self::None => 1.0,
            Self::Zstd(level) => {
                // Rough estimates based on typical CFD data
                match level {
                    1..=3 => 0.7,   // Fast compression
                    4..=9 => 0.5,   // Balanced
                    10..=15 => 0.4, // High compression
                    _ => 0.3,       // Maximum compression
                }
            }
        }
    }
}
