//! Strategy selection thresholds and parameters

/// Strategy selection thresholds and parameters
///
/// These are plain compile-time constants, not wrapped parameters: nothing in
/// the crate ever mutates, validates or adapts them, so the parameter wrapper
/// only ever contributed allocation and indirection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrategyThresholds {
    /// Minimum curvature factor to use arc strategy instead of straight
    pub arc_curvature_threshold: f64,

    /// Maximum fill factor to use serpentine strategy
    pub serpentine_fill_threshold: f64,

    /// Minimum channel length for complex strategies
    pub min_complex_strategy_length: f64,

    /// Branch count threshold for adaptive behavior
    pub adaptive_branch_threshold: usize,

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
        arc_curvature_threshold: 0.1,
        serpentine_fill_threshold: 0.95,
        min_complex_strategy_length: 10.0,
        adaptive_branch_threshold: 4,
        frustum_min_length_threshold: 0.3,
        frustum_max_length_threshold: 0.7,
        frustum_angle_threshold: 0.5,
    };
}

impl Default for StrategyThresholds {
    fn default() -> Self {
        Self::DEFAULT
    }
}
