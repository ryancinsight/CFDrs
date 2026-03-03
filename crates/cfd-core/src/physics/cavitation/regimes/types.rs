//! Cavitation regime types.

use serde::{Deserialize, Serialize};

/// Cavitation regime types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CavitationRegime {
    /// No cavitation occurring
    None,
    /// Stable oscillating cavitation
    Stable,
    /// Transient inertial cavitation
    Inertial,
    /// Mixed regime (both stable and inertial present)
    Mixed,
}
