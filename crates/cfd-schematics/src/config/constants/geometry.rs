//! Geometry generation constants for point and channel generation

use super::primitives;

/// Geometry generation constants previously hardcoded
///
/// The values themselves live in [`primitives`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryGenerationConstants {
    /// Default wall clearance
    pub default_wall_clearance: f64,

    /// Channel width multiplier for short channel detection
    pub short_channel_width_multiplier: f64,

    /// Horizontal angle threshold for strategy selection
    pub horizontal_angle_threshold: f64,

    /// Long horizontal threshold for strategy selection
    pub long_horizontal_threshold: f64,

    /// Minimum arc length threshold for strategy selection
    pub min_arc_length_threshold: f64,

    /// Maximum curvature reduction factor for adaptive arcs
    pub max_curvature_reduction_factor: f64,

    /// Minimum curvature factor for adaptive arcs
    pub min_curvature_factor: f64,
}

impl GeometryGenerationConstants {
    /// Canonical default geometry generation constants.
    pub const DEFAULT: Self = Self {
        default_wall_clearance: primitives::DEFAULT_WALL_CLEARANCE,
        short_channel_width_multiplier: primitives::SHORT_CHANNEL_WIDTH_MULTIPLIER,
        horizontal_angle_threshold: primitives::strategy_thresholds::HORIZONTAL_ANGLE_THRESHOLD,
        long_horizontal_threshold: primitives::strategy_thresholds::LONG_HORIZONTAL_THRESHOLD,
        min_arc_length_threshold: primitives::strategy_thresholds::MIN_ARC_LENGTH_THRESHOLD,
        max_curvature_reduction_factor: primitives::DEFAULT_MAX_CURVATURE_REDUCTION,
        min_curvature_factor: primitives::MIN_ADAPTIVE_CURVATURE_FACTOR,
    };
}

impl Default for GeometryGenerationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
