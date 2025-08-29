//! Checkpoint metadata structures

use serde::{Deserialize, Serialize};

/// Version for checkpoint format compatibility
pub const CHECKPOINT_VERSION: u32 = 1;

/// Checkpoint metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckpointMetadata {
    /// Version of checkpoint format
    pub version: u32,
    /// Simulation time
    pub time: f64,
    /// Iteration number
    pub iteration: usize,
    /// Grid dimensions
    pub dimensions: (usize, usize),
    /// Physical domain size
    pub domain_size: (f64, f64),
    /// Timestamp of checkpoint creation
    pub timestamp: u64,
}

impl CheckpointMetadata {
    /// Create new metadata
    pub fn new(
        time: f64,
        iteration: usize,
        dimensions: (usize, usize),
        domain_size: (f64, f64),
    ) -> Self {
        Self {
            version: CHECKPOINT_VERSION,
            time,
            iteration,
            dimensions,
            domain_size,
            timestamp: std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0),
        }
    }

    /// Check if version is compatible
    pub fn is_version_compatible(&self) -> bool {
        self.version == CHECKPOINT_VERSION
    }

    /// Validate metadata consistency
    pub fn validate(&self) -> Result<(), String> {
        if !self.is_version_compatible() {
            return Err(format!(
                "Incompatible checkpoint version: {} (expected {})",
                self.version, CHECKPOINT_VERSION
            ));
        }

        if self.dimensions.0 == 0 || self.dimensions.1 == 0 {
            return Err("Invalid grid dimensions".to_string());
        }

        if self.domain_size.0 <= 0.0 || self.domain_size.1 <= 0.0 {
            return Err("Invalid domain size".to_string());
        }

        Ok(())
    }
}
