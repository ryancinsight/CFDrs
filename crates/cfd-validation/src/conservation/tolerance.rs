//! Tolerance settings for conservation checks

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Tolerance settings for conservation checks
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationTolerance<T: RealField + Copy> {
    /// Absolute tolerance
    pub absolute: T,
    /// Relative tolerance
    pub relative: T,
}

impl<T: RealField + Copy> ConservationTolerance<T> {
    /// Create new tolerance settings
    pub fn new(absolute: T, relative: T) -> Self {
        Self { absolute, relative }
    }

    /// Check if error is within tolerance
    pub fn is_satisfied(&self, error: T, reference: T) -> bool {
        error <= self.absolute || error <= self.relative * reference.abs()
    }
}
