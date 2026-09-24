//! Geometry generation constants for point and channel generation

/// Geometry generation constants previously hardcoded
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryGenerationConstants {
    /// Default number of points for serpentine path generation
    pub default_serpentine_points: usize,

    /// Minimum number of points for serpentine path generation
    pub min_serpentine_points: usize,

    /// Maximum number of points for serpentine path generation
    pub max_serpentine_points: usize,

    /// Default wall clearance
    pub default_wall_clearance: f64,

    /// Default channel width
    pub default_channel_width: f64,

    /// Default channel height
    pub default_channel_height: f64,

    /// Channel width multiplier for short channel detection
    pub short_channel_width_multiplier: f64,

    /// Default middle points for smooth straight channels
    pub smooth_straight_middle_points: usize,

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
        default_serpentine_points: 200,
        min_serpentine_points: 10,
        max_serpentine_points: 1000,
        default_wall_clearance: 0.5,
        default_channel_width: 1.0,
        default_channel_height: 1.0,
        short_channel_width_multiplier: 2.0,
        smooth_straight_middle_points: 10,
        horizontal_angle_threshold: 0.5,
        long_horizontal_threshold: 0.6,
        min_arc_length_threshold: 0.3,
        max_curvature_reduction_factor: 0.5,
        min_curvature_factor: 0.1,
    };
}

impl Default for GeometryGenerationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
