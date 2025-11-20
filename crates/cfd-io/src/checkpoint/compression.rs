//! Compression strategies for checkpoints

use serde::{Deserialize, Serialize};

/// Compression ratio constants for different compression levels
mod compression_ratios {
    /// No compression ratio
    pub const NONE: f64 = 1.0;
    /// Low compression level ratio (levels 1-3)
    pub const LOW_LEVELS: f64 = 0.7;
    /// Balanced compression ratio (levels 4-9)
    pub const BALANCED: f64 = 0.5;
    /// High compression ratio (levels 10-15)
    pub const HIGH: f64 = 0.4;
    /// Maximum compression ratio (levels 16+)
    pub const MAXIMUM: f64 = 0.3;
}

/// Compression strategy for checkpoints
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
pub enum CompressionStrategy {
    /// No compression
    #[default]
    None,
    /// Zstandard compression with level (1-22, default 3)
    Zstd(u8),
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
            Self::None => compression_ratios::NONE,
            Self::Zstd(level) => {
                // Rough estimates based on typical CFD data
                match level {
                    1..=3 => compression_ratios::LOW_LEVELS, // Low compression levels
                    4..=9 => compression_ratios::BALANCED,   // Balanced
                    10..=15 => compression_ratios::HIGH,     // High compression
                    _ => compression_ratios::MAXIMUM,        // Maximum compression
                }
            }
        }
    }
}
