//! Strategy selection thresholds and parameters

use super::primitives;

/// Strategy selection thresholds and parameters
///
/// Plain compile-time constants: nothing in the crate ever mutates, validates
/// or adapts them, so the parameter wrapper only ever contributed allocation
/// and indirection. The values themselves live in [`primitives`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrategyThresholds {
    /// Minimum length threshold for frustum channel selection in smart mode
    pub frustum_min_length_threshold: f64,

    /// Maximum length threshold for frustum channel selection in smart mode
    pub frustum_max_length_threshold: f64,

    /// Maximum angle threshold for frustum channel selection (horizontal preference)
    pub frustum_angle_threshold: f64,
}

impl StrategyThresholds {
    /// Canonical default thresholds.
    pub const DEFAULT: Self = Self {
        frustum_min_length_threshold: primitives::FRUSTUM_MIN_LENGTH_THRESHOLD,
        frustum_max_length_threshold: primitives::FRUSTUM_MAX_LENGTH_THRESHOLD,
        frustum_angle_threshold: primitives::FRUSTUM_ANGLE_THRESHOLD,
    };
}

impl Default for StrategyThresholds {
    fn default() -> Self {
        Self::DEFAULT
    }
}
