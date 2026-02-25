//! Geometry generation constants for point and channel generation

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

/// Geometry generation constants previously hardcoded
pub struct GeometryGenerationConstants {
    /// Default number of points for serpentine path generation
    pub default_serpentine_points: ConfigurableParameter<usize>,

    /// Minimum number of points for serpentine path generation
    pub min_serpentine_points: ConfigurableParameter<usize>,

    /// Maximum number of points for serpentine path generation
    pub max_serpentine_points: ConfigurableParameter<usize>,

    /// Default wall clearance
    pub default_wall_clearance: ConfigurableParameter<f64>,

    /// Default channel width
    pub default_channel_width: ConfigurableParameter<f64>,

    /// Default channel height
    pub default_channel_height: ConfigurableParameter<f64>,

    /// Channel width multiplier for short channel detection
    pub short_channel_width_multiplier: ConfigurableParameter<f64>,

    /// Default middle points for smooth straight channels
    pub smooth_straight_middle_points: ConfigurableParameter<usize>,

    /// Horizontal angle threshold for strategy selection
    pub horizontal_angle_threshold: ConfigurableParameter<f64>,

    /// Long horizontal threshold for strategy selection
    pub long_horizontal_threshold: ConfigurableParameter<f64>,

    /// Minimum arc length threshold for strategy selection
    pub min_arc_length_threshold: ConfigurableParameter<f64>,

    /// Maximum curvature reduction factor for adaptive arcs
    pub max_curvature_reduction_factor: ConfigurableParameter<f64>,

    /// Minimum curvature factor for adaptive arcs
    pub min_curvature_factor: ConfigurableParameter<f64>,
}

impl GeometryGenerationConstants {
    /// Create default geometry generation constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            default_serpentine_points: ConfigurableParameter::new(
                200usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(10, 1000),
                ]),
                ParameterMetadata::new(
                    "default_serpentine_points",
                    "Default number of points for serpentine path generation",
                    "geometry_generation",
                )
                .with_units("points"),
            ),

            min_serpentine_points: ConfigurableParameter::new(
                10usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(3, 100),
                ]),
                ParameterMetadata::new(
                    "min_serpentine_points",
                    "Minimum number of points for serpentine path generation",
                    "geometry_generation",
                )
                .with_units("points"),
            ),

            max_serpentine_points: ConfigurableParameter::new(
                1000usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(100, 5000),
                ]),
                ParameterMetadata::new(
                    "max_serpentine_points",
                    "Maximum number of points for serpentine path generation",
                    "geometry_generation",
                )
                .with_units("points"),
            ),

            default_wall_clearance: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 100.0),
                ]),
                ParameterMetadata::new(
                    "default_wall_clearance",
                    "Default wall clearance for geometry generation",
                    "geometry_generation",
                )
                .with_units("mm"),
            ),

            default_channel_width: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 50.0),
                ]),
                ParameterMetadata::new(
                    "default_channel_width",
                    "Default channel width for geometry generation",
                    "geometry_generation",
                )
                .with_units("mm"),
            ),

            default_channel_height: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 50.0),
                ]),
                ParameterMetadata::new(
                    "default_channel_height",
                    "Default channel height for geometry generation",
                    "geometry_generation",
                )
                .with_units("mm"),
            ),

            short_channel_width_multiplier: ConfigurableParameter::new(
                2.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 10.0),
                ]),
                ParameterMetadata::new(
                    "short_channel_width_multiplier",
                    "Multiplier for channel width to detect short channels",
                    "geometry_generation",
                )
                .with_units("factor"),
            ),

            smooth_straight_middle_points: ConfigurableParameter::new(
                10,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(5, 100),
                ]),
                ParameterMetadata::new(
                    "smooth_straight_middle_points",
                    "Number of middle points for smooth straight channels",
                    "geometry_generation",
                )
                .with_units("points"),
            ),

            horizontal_angle_threshold: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 2.0),
                ]),
                ParameterMetadata::new(
                    "horizontal_angle_threshold",
                    "Threshold for detecting horizontal channels in strategy selection",
                    "geometry_generation",
                )
                .with_units("ratio"),
            ),

            long_horizontal_threshold: ConfigurableParameter::new(
                0.6,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 1.0),
                ]),
                ParameterMetadata::new(
                    "long_horizontal_threshold",
                    "Threshold for detecting long horizontal channels",
                    "geometry_generation",
                )
                .with_units("ratio"),
            ),

            min_arc_length_threshold: ConfigurableParameter::new(
                0.3,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 1.0),
                ]),
                ParameterMetadata::new(
                    "min_arc_length_threshold",
                    "Minimum length threshold for arc channel selection",
                    "geometry_generation",
                )
                .with_units("ratio"),
            ),

            max_curvature_reduction_factor: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "max_curvature_reduction_factor",
                    "Maximum curvature reduction factor for adaptive arcs",
                    "geometry_generation",
                )
                .with_units("factor"),
            ),

            min_curvature_factor: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 1.0),
                ]),
                ParameterMetadata::new(
                    "min_curvature_factor",
                    "Minimum curvature factor for adaptive arcs",
                    "geometry_generation",
                )
                .with_units("factor"),
            ),
        }
    }
}

impl Default for GeometryGenerationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}
