//! Centralized configuration constants
//!
//! This module extracts all hardcoded values and magic numbers from throughout
//! the codebase into configurable parameters with proper validation and documentation.
//! This follows the SSOT (Single Source of Truth) principle and eliminates magic numbers.

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

/// Configuration constants for geometry validation and defaults (Primitives)
pub mod primitives {
    /// Minimum allowed wall clearance (mm)
    pub const MIN_WALL_CLEARANCE: f64 = 0.1;
    /// Maximum allowed wall clearance (mm)
    pub const MAX_WALL_CLEARANCE: f64 = 100.0;
    /// Default wall clearance (mm)
    pub const DEFAULT_WALL_CLEARANCE: f64 = 0.5;

    // Rendering and geometry generation constants
    /// Default number of points for serpentine path generation
    pub const DEFAULT_SERPENTINE_POINTS: usize = 200;
    /// Minimum number of points for serpentine path generation
    pub const MIN_SERPENTINE_POINTS: usize = 10;
    /// Maximum number of points for serpentine path generation
    pub const MAX_SERPENTINE_POINTS: usize = 1000;

    /// Default number of points for optimization path generation
    pub const DEFAULT_OPTIMIZATION_POINTS: usize = 50;
    /// Minimum number of points for optimization path generation
    pub const MIN_OPTIMIZATION_POINTS: usize = 10;
    /// Maximum number of points for optimization path generation
    pub const MAX_OPTIMIZATION_POINTS: usize = 200;

    /// Default number of middle points for smooth straight channels
    pub const DEFAULT_SMOOTH_STRAIGHT_MIDDLE_POINTS: usize = 10;
    /// Minimum number of middle points for smooth straight channels
    pub const MIN_SMOOTH_STRAIGHT_MIDDLE_POINTS: usize = 2;
    /// Maximum number of middle points for smooth straight channels
    pub const MAX_SMOOTH_STRAIGHT_MIDDLE_POINTS: usize = 50;

    /// Default wave multiplier for smooth transitions (2π for one complete wave)
    pub const DEFAULT_TRANSITION_WAVE_MULTIPLIER: f64 = 2.0;
    /// Minimum wave multiplier for smooth transitions
    pub const MIN_TRANSITION_WAVE_MULTIPLIER: f64 = 0.5;
    /// Maximum wave multiplier for smooth transitions
    pub const MAX_TRANSITION_WAVE_MULTIPLIER: f64 = 10.0;

    // Adaptive serpentine control constants
    /// Default distance normalization factor for node proximity effects
    pub const DEFAULT_NODE_DISTANCE_NORMALIZATION: f64 = 10.0;
    /// Minimum distance normalization factor
    pub const MIN_NODE_DISTANCE_NORMALIZATION: f64 = 1.0;
    /// Maximum distance normalization factor
    pub const MAX_NODE_DISTANCE_NORMALIZATION: f64 = 50.0;

    /// Default plateau width factor for horizontal channels (fraction of channel length)
    pub const DEFAULT_PLATEAU_WIDTH_FACTOR: f64 = 0.4;
    /// Minimum plateau width factor
    pub const MIN_PLATEAU_WIDTH_FACTOR: f64 = 0.1;
    /// Maximum plateau width factor
    pub const MAX_PLATEAU_WIDTH_FACTOR: f64 = 0.8;

    /// Default horizontal ratio threshold for middle section detection
    pub const DEFAULT_HORIZONTAL_RATIO_THRESHOLD: f64 = 1.01;
    /// Minimum horizontal ratio threshold
    pub const MIN_HORIZONTAL_RATIO_THRESHOLD: f64 = 0.5;
    /// Maximum horizontal ratio threshold
    pub const MAX_HORIZONTAL_RATIO_THRESHOLD: f64 = 0.95;

    /// Default middle section amplitude factor
    pub const DEFAULT_MIDDLE_SECTION_AMPLITUDE_FACTOR: f64 = 0.7;
    /// Minimum middle section amplitude factor
    pub const MIN_MIDDLE_SECTION_AMPLITUDE_FACTOR: f64 = 0.1;
    /// Maximum middle section amplitude factor
    pub const MAX_MIDDLE_SECTION_AMPLITUDE_FACTOR: f64 = 1.0;

    /// Default plateau amplitude factor
    pub const DEFAULT_PLATEAU_AMPLITUDE_FACTOR: f64 = 0.8;
    /// Minimum plateau amplitude factor
    pub const MIN_PLATEAU_AMPLITUDE_FACTOR: f64 = 0.5;
    /// Maximum plateau amplitude factor
    pub const MAX_PLATEAU_AMPLITUDE_FACTOR: f64 = 1.0;

    /// Minimum allowed channel width (mm)
    pub const MIN_CHANNEL_WIDTH: f64 = 0.01;
    /// Maximum allowed channel width (mm)
    pub const MAX_CHANNEL_WIDTH: f64 = 1000.0;
    /// Default channel width (mm)
    pub const DEFAULT_CHANNEL_WIDTH: f64 = 1.0;

    /// Minimum allowed channel height (mm)
    pub const MIN_CHANNEL_HEIGHT: f64 = 0.01;
    /// Maximum allowed channel height (mm)
    pub const MAX_CHANNEL_HEIGHT: f64 = 1000.0;
    /// Default channel height (mm)
    pub const DEFAULT_CHANNEL_HEIGHT: f64 = 0.5;

    /// Minimum fill factor for serpentine channels
    pub const MIN_FILL_FACTOR: f64 = 0.1;
    /// Maximum fill factor for serpentine channels
    pub const MAX_FILL_FACTOR: f64 = 0.95;
    /// Default fill factor for serpentine channels
    pub const DEFAULT_FILL_FACTOR: f64 = 0.8;

    /// Minimum wavelength factor for serpentine channels
    pub const MIN_WAVELENGTH_FACTOR: f64 = 1.0;
    /// Maximum wavelength factor for serpentine channels
    pub const MAX_WAVELENGTH_FACTOR: f64 = 10.0;
    /// Default wavelength factor for serpentine channels
    pub const DEFAULT_WAVELENGTH_FACTOR: f64 = 6.0;

    /// Minimum Gaussian width factor for serpentine channels
    pub const MIN_GAUSSIAN_WIDTH_FACTOR: f64 = 2.0;
    /// Maximum Gaussian width factor for serpentine channels
    pub const MAX_GAUSSIAN_WIDTH_FACTOR: f64 = 20.0;
    /// Default Gaussian width factor for serpentine channels
    pub const DEFAULT_GAUSSIAN_WIDTH_FACTOR: f64 = 6.0;

    /// Minimum wave density factor for serpentine channels
    pub const MIN_WAVE_DENSITY_FACTOR: f64 = 0.5;
    /// Maximum wave density factor for serpentine channels
    pub const MAX_WAVE_DENSITY_FACTOR: f64 = 10.0;
    /// Default wave density factor for serpentine channels
    pub const DEFAULT_WAVE_DENSITY_FACTOR: f64 = 1.5;

    /// Minimum curvature factor for arc channels
    pub const MIN_CURVATURE_FACTOR: f64 = 0.0;
    /// Maximum curvature factor for arc channels
    pub const MAX_CURVATURE_FACTOR: f64 = 2.0;
    /// Default curvature factor for arc channels
    pub const DEFAULT_CURVATURE_FACTOR: f64 = 0.3;

    /// Minimum smoothness for arc channels
    pub const MIN_SMOOTHNESS: usize = 3;
    /// Maximum smoothness for arc channels
    pub const MAX_SMOOTHNESS: usize = 1000;
    /// Default smoothness for arc channels
    pub const DEFAULT_SMOOTHNESS: usize = 20;

    /// Minimum separation distance between arc channels (mm)
    pub const MIN_SEPARATION_DISTANCE: f64 = 0.1;
    /// Maximum separation distance between arc channels (mm)
    pub const MAX_SEPARATION_DISTANCE: f64 = 10.0;
    /// Default minimum separation distance between arc channels (mm)
    pub const DEFAULT_MIN_SEPARATION_DISTANCE: f64 = 1.0;

    /// Minimum curvature reduction factor for collision prevention
    pub const MIN_CURVATURE_REDUCTION: f64 = 0.1;
    /// Maximum curvature reduction factor for collision prevention
    pub const MAX_CURVATURE_REDUCTION_LIMIT: f64 = 1.0;
    /// Default maximum curvature reduction factor
    pub const DEFAULT_MAX_CURVATURE_REDUCTION: f64 = 0.5;

    /// Strategy thresholds for smart channel type selection
    pub mod strategy_thresholds {
        /// Threshold for long horizontal channels (fraction of box width)
        pub const LONG_HORIZONTAL_THRESHOLD: f64 = 0.3;
        /// Threshold for minimum arc length (fraction of box width)
        pub const MIN_ARC_LENGTH_THRESHOLD: f64 = 0.1;
        /// Threshold for horizontal vs angled channel detection
        pub const HORIZONTAL_ANGLE_THRESHOLD: f64 = 0.3;
        /// Threshold for angled channel detection (slope)
        pub const ANGLED_CHANNEL_SLOPE_THRESHOLD: f64 = 0.1;
        /// Default middle zone fraction for mixed by position
        pub const DEFAULT_MIDDLE_ZONE_FRACTION: f64 = 0.4;
    }

    // Adaptive collision constants (Added during refactor recovery)
    pub const DEFAULT_MAX_ADJUSTMENT_FACTOR: f64 = 0.5;
    pub const DEFAULT_MIN_CHANNEL_DISTANCE: f64 = 0.5;
    pub const DEFAULT_SAFETY_MARGIN_FACTOR: f64 = 1.1;
    pub const DEFAULT_MAX_REDUCTION_FACTOR: f64 = 0.8;
    pub const DEFAULT_DETECTION_SENSITIVITY: f64 = 0.1;
    pub const DEFAULT_PROXIMITY_DIVISOR: f64 = 10.0;
    pub const DEFAULT_MIN_PROXIMITY_FACTOR: f64 = 0.2;
    pub const DEFAULT_MAX_PROXIMITY_FACTOR: f64 = 2.0;
    pub const DEFAULT_BRANCH_ADJUSTMENT_DIVISOR: f64 = 5.0;
    pub const DEFAULT_MAX_SENSITIVITY_MULTIPLIER: f64 = 3.0;
    pub const DEFAULT_LONG_CHANNEL_REDUCTION_MULTIPLIER: f64 = 0.9;
    pub const DEFAULT_MAX_REDUCTION_LIMIT: f64 = 0.5;
}

/// Strategy selection thresholds and parameters
pub struct StrategyThresholds {
    /// Minimum curvature factor to use arc strategy instead of straight
    pub arc_curvature_threshold: ConfigurableParameter<f64>,

    /// Maximum fill factor to use serpentine strategy
    pub serpentine_fill_threshold: ConfigurableParameter<f64>,

    /// Minimum channel length for complex strategies
    pub min_complex_strategy_length: ConfigurableParameter<f64>,

    /// Branch count threshold for adaptive behavior
    pub adaptive_branch_threshold: ConfigurableParameter<usize>,

    /// Minimum length threshold for frustum channel selection in smart mode
    pub frustum_min_length_threshold: ConfigurableParameter<f64>,

    /// Maximum length threshold for frustum channel selection in smart mode
    pub frustum_max_length_threshold: ConfigurableParameter<f64>,

    /// Maximum angle threshold for frustum channel selection (horizontal preference)
    pub frustum_angle_threshold: ConfigurableParameter<f64>,
}

impl StrategyThresholds {
    /// Create default strategy thresholds
    #[must_use]
    fn make_default() -> Self {
        Self {
            arc_curvature_threshold: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::non_negative(),
                    ParameterConstraints::range(0.0, 2.0),
                ]),
                ParameterMetadata::new(
                    "arc_curvature_threshold",
                    "Minimum curvature factor to trigger arc strategy selection",
                    "strategy_selection"
                ).with_units("factor").affects_others()
            ),
            
            serpentine_fill_threshold: ConfigurableParameter::new(
                0.95,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "serpentine_fill_threshold",
                    "Maximum fill factor for serpentine strategy selection",
                    "strategy_selection"
                ).affects_others()
            ),
            
            min_complex_strategy_length: ConfigurableParameter::new(
                10.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 1000.0),
                ]),
                ParameterMetadata::new(
                    "min_complex_strategy_length",
                    "Minimum channel length to use complex strategies",
                    "strategy_selection"
                ).with_units("mm")
            ),
            
            adaptive_branch_threshold: ConfigurableParameter::new(
                4usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(1, 100),
                ]),
                ParameterMetadata::new(
                    "adaptive_branch_threshold",
                    "Branch count threshold to enable adaptive parameter behavior",
                    "strategy_selection"
                )
            ),

            frustum_min_length_threshold: ConfigurableParameter::new(
                0.3,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 1.0),
                ]),
                ParameterMetadata::new(
                    "frustum_min_length_threshold",
                    "Minimum length threshold (as fraction of box width) for frustum channel selection",
                    "strategy_selection"
                ).with_units("fraction")
            ),

            frustum_max_length_threshold: ConfigurableParameter::new(
                0.7,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.2, 1.0),
                ]),
                ParameterMetadata::new(
                    "frustum_max_length_threshold",
                    "Maximum length threshold (as fraction of box width) for frustum channel selection",
                    "strategy_selection"
                ).with_units("fraction")
            ),

            frustum_angle_threshold: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 2.0),
                ]),
                ParameterMetadata::new(
                    "frustum_angle_threshold",
                    "Maximum angle threshold for frustum channel selection (dy/dx ratio)",
                    "strategy_selection"
                ).with_units("ratio")
            ),
        }
    }
}

/// Wave generation constants previously hardcoded in strategies
pub struct WaveGenerationConstants {
    /// Sharpness factor for square wave generation
    pub square_wave_sharpness: ConfigurableParameter<f64>,

    /// Transition zone factor for smooth endpoints
    pub transition_zone_factor: ConfigurableParameter<f64>,

    /// Gaussian envelope scaling factor
    pub gaussian_envelope_scale: ConfigurableParameter<f64>,

    /// Phase direction calculation threshold
    pub phase_direction_threshold: ConfigurableParameter<f64>,

    /// Wave amplitude safety margin
    pub amplitude_safety_margin: ConfigurableParameter<f64>,

    /// Smooth endpoint transition start threshold
    pub smooth_endpoint_start_threshold: ConfigurableParameter<f64>,

    /// Smooth endpoint transition end threshold
    pub smooth_endpoint_end_threshold: ConfigurableParameter<f64>,

    /// Default transition length factor for smooth transitions
    pub default_transition_length_factor: ConfigurableParameter<f64>,

    /// Default transition amplitude factor
    pub default_transition_amplitude_factor: ConfigurableParameter<f64>,

    /// Default transition smoothness points
    pub default_transition_smoothness: ConfigurableParameter<usize>,

    /// Default wave multiplier for transitions
    pub default_wave_multiplier: ConfigurableParameter<f64>,

    /// Wall proximity scaling factor
    pub wall_proximity_scaling_factor: ConfigurableParameter<f64>,

    /// Neighbor avoidance scaling factor
    pub neighbor_avoidance_scaling_factor: ConfigurableParameter<f64>,

    /// Geometric tolerance for distance comparisons
    pub geometric_tolerance: ConfigurableParameter<f64>,
}

impl Default for StrategyThresholds {
    fn default() -> Self {
        Self::make_default()
    }
}

impl WaveGenerationConstants {
    /// Create default wave generation constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            square_wave_sharpness: ConfigurableParameter::new(
                5.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 20.0),
                ]),
                ParameterMetadata::new(
                    "square_wave_sharpness",
                    "Sharpness factor for square wave generation using tanh",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            transition_zone_factor: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 0.5),
                ]),
                ParameterMetadata::new(
                    "transition_zone_factor",
                    "Factor for smooth transition zones at wave endpoints",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            gaussian_envelope_scale: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 5.0),
                ]),
                ParameterMetadata::new(
                    "gaussian_envelope_scale",
                    "Scaling factor for Gaussian envelope calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            phase_direction_threshold: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "phase_direction_threshold",
                    "Threshold for phase direction calculation in bilateral symmetry",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            amplitude_safety_margin: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "amplitude_safety_margin",
                    "Safety margin factor for amplitude calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            smooth_endpoint_start_threshold: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 0.5),
                ]),
                ParameterMetadata::new(
                    "smooth_endpoint_start_threshold",
                    "Threshold for smooth endpoint transition start",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            smooth_endpoint_end_threshold: ConfigurableParameter::new(
                0.9,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.5, 0.99),
                ]),
                ParameterMetadata::new(
                    "smooth_endpoint_end_threshold",
                    "Threshold for smooth endpoint transition end",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            default_transition_length_factor: ConfigurableParameter::new(
                0.15,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.05, 0.5),
                ]),
                ParameterMetadata::new(
                    "default_transition_length_factor",
                    "Default length factor for smooth transitions",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            default_transition_amplitude_factor: ConfigurableParameter::new(
                0.3,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "default_transition_amplitude_factor",
                    "Default amplitude factor for smooth transitions",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            default_transition_smoothness: ConfigurableParameter::new(
                20,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(5, 100),
                ]),
                ParameterMetadata::new(
                    "default_transition_smoothness",
                    "Default number of points for transition smoothing",
                    "wave_generation",
                )
                .with_units("points"),
            ),

            default_wave_multiplier: ConfigurableParameter::new(
                2.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.5, 10.0),
                ]),
                ParameterMetadata::new(
                    "default_wave_multiplier",
                    "Default wave multiplier for transitions",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            wall_proximity_scaling_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "wall_proximity_scaling_factor",
                    "Scaling factor for wall proximity calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            neighbor_avoidance_scaling_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "neighbor_avoidance_scaling_factor",
                    "Scaling factor for neighbor avoidance calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            geometric_tolerance: ConfigurableParameter::new(
                1e-6,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1e-12, 1e-3),
                ]),
                ParameterMetadata::new(
                    "geometric_tolerance",
                    "Tolerance for geometric distance comparisons",
                    "wave_generation",
                )
                .with_units("units"),
            ),
        }
    }
}

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

impl Default for WaveGenerationConstants {
    fn default() -> Self {
        Self::make_default()
    }
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

/// Optimization algorithm constants previously hardcoded
pub struct OptimizationConstants {
    /// Branch factor scaling exponent (was hardcoded as 0.75)
    pub branch_factor_exponent: ConfigurableParameter<f64>,

    /// Fill factor enhancement multiplier (was hardcoded as 1.5)
    pub fill_factor_enhancement: ConfigurableParameter<f64>,

    /// Maximum optimization iterations
    pub max_optimization_iterations: ConfigurableParameter<usize>,

    /// Convergence tolerance for optimization
    pub convergence_tolerance: ConfigurableParameter<f64>,

    /// Fast optimization wavelength factors
    pub fast_wavelength_factors: ConfigurableParameter<Vec<f64>>,

    /// Fast optimization wave density factors
    pub fast_wave_density_factors: ConfigurableParameter<Vec<f64>>,

    /// Fast optimization fill factors
    pub fast_fill_factors: ConfigurableParameter<Vec<f64>>,
}

impl Default for GeometryGenerationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}

impl OptimizationConstants {
    /// Create default optimization constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            branch_factor_exponent: ConfigurableParameter::new(
                0.75,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 2.0),
                ]),
                ParameterMetadata::new(
                    "branch_factor_exponent",
                    "Exponent for branch factor scaling (was hardcoded as 0.75)",
                    "optimization",
                )
                .with_units("exponent"),
            ),

            fill_factor_enhancement: ConfigurableParameter::new(
                1.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 3.0),
                ]),
                ParameterMetadata::new(
                    "fill_factor_enhancement",
                    "Enhancement multiplier for fill factor (was hardcoded as 1.5)",
                    "optimization",
                )
                .with_units("multiplier"),
            ),

            max_optimization_iterations: ConfigurableParameter::new(
                100usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(10, 1000),
                ]),
                ParameterMetadata::new(
                    "max_optimization_iterations",
                    "Maximum number of iterations for optimization algorithms",
                    "optimization",
                )
                .with_units("iterations"),
            ),

            convergence_tolerance: ConfigurableParameter::new(
                1e-6,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1e-10, 1e-2),
                ]),
                ParameterMetadata::new(
                    "convergence_tolerance",
                    "Tolerance for optimization convergence detection",
                    "optimization",
                )
                .with_units("tolerance"),
            ),

            fast_wavelength_factors: ConfigurableParameter::new(
                vec![1.0, 2.0, 3.0, 4.0],
                ParameterConstraints::all(vec![]),
                ParameterMetadata::new(
                    "fast_wavelength_factors",
                    "Wavelength factors for fast optimization (was hardcoded array)",
                    "optimization",
                )
                .with_units("factors"),
            ),

            fast_wave_density_factors: ConfigurableParameter::new(
                vec![1.0, 2.0, 3.0],
                ParameterConstraints::all(vec![]),
                ParameterMetadata::new(
                    "fast_wave_density_factors",
                    "Wave density factors for fast optimization (was hardcoded array)",
                    "optimization",
                )
                .with_units("factors"),
            ),

            fast_fill_factors: ConfigurableParameter::new(
                vec![0.7, 0.8, 0.9],
                ParameterConstraints::all(vec![]),
                ParameterMetadata::new(
                    "fast_fill_factors",
                    "Fill factors for fast optimization (was hardcoded array)",
                    "optimization",
                )
                .with_units("factors"),
            ),
        }
    }
}

/// Visualization constants previously hardcoded
pub struct VisualizationConstants {
    /// Default margin for chart rendering
    pub default_chart_margin: ConfigurableParameter<u32>,

    /// Default right margin for chart rendering
    pub default_chart_right_margin: ConfigurableParameter<u32>,

    /// Default label area size for x-axis
    pub default_x_label_area_size: ConfigurableParameter<u32>,

    /// Default label area size for y-axis
    pub default_y_label_area_size: ConfigurableParameter<u32>,

    /// Default buffer factor for chart boundaries
    pub default_boundary_buffer_factor: ConfigurableParameter<f64>,
}

impl Default for OptimizationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}

impl VisualizationConstants {
    /// Create default visualization constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            default_chart_margin: ConfigurableParameter::new(
                20u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(5, 100),
                ]),
                ParameterMetadata::new(
                    "default_chart_margin",
                    "Default margin for chart rendering",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_chart_right_margin: ConfigurableParameter::new(
                150u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(50, 500),
                ]),
                ParameterMetadata::new(
                    "default_chart_right_margin",
                    "Default right margin for chart rendering",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_x_label_area_size: ConfigurableParameter::new(
                30u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(10, 100),
                ]),
                ParameterMetadata::new(
                    "default_x_label_area_size",
                    "Default label area size for x-axis",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_y_label_area_size: ConfigurableParameter::new(
                30u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(10, 100),
                ]),
                ParameterMetadata::new(
                    "default_y_label_area_size",
                    "Default label area size for y-axis",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_boundary_buffer_factor: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 0.5),
                ]),
                ParameterMetadata::new(
                    "default_boundary_buffer_factor",
                    "Default buffer factor for chart boundaries",
                    "visualization",
                )
                .with_units("ratio"),
            ),
        }
    }
}

impl Default for VisualizationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}


/// Registry for all configuration constants
pub struct ConstantsRegistry {
    pub strategies: StrategyThresholds,
    pub waves: WaveGenerationConstants,
    pub geometry: GeometryGenerationConstants,
    pub optimization: OptimizationConstants,
    pub visualization: VisualizationConstants,
}

impl ConstantsRegistry {
    pub fn new() -> Self {
        Self {
            strategies: StrategyThresholds::default(),
            waves: WaveGenerationConstants::default(),
            geometry: GeometryGenerationConstants::default(),
            optimization: OptimizationConstants::default(),
            visualization: VisualizationConstants::default(),
        }
    }

    // --- Adaptive Collision ---
    pub fn get_max_adjustment_factor(&self) -> f64 { primitives::DEFAULT_MAX_ADJUSTMENT_FACTOR }
    pub fn get_min_channel_distance(&self) -> f64 { primitives::DEFAULT_MIN_CHANNEL_DISTANCE }
    pub fn get_min_wall_distance(&self) -> f64 { self.geometry.default_wall_clearance.get_raw_value().clone() }
    pub fn get_safety_margin_factor(&self) -> f64 { primitives::DEFAULT_SAFETY_MARGIN_FACTOR }
    pub fn get_max_reduction_factor(&self) -> f64 { primitives::DEFAULT_MAX_REDUCTION_FACTOR }
    pub fn get_detection_sensitivity(&self) -> f64 { primitives::DEFAULT_DETECTION_SENSITIVITY }
    pub fn get_proximity_divisor(&self) -> f64 { primitives::DEFAULT_PROXIMITY_DIVISOR }
    pub fn get_min_proximity_factor(&self) -> f64 { primitives::DEFAULT_MIN_PROXIMITY_FACTOR }
    pub fn get_max_proximity_factor(&self) -> f64 { primitives::DEFAULT_MAX_PROXIMITY_FACTOR }
    pub fn get_branch_factor_exponent(&self) -> f64 { *self.optimization.branch_factor_exponent.get_raw_value() }
    pub fn get_branch_adjustment_divisor(&self) -> f64 { primitives::DEFAULT_BRANCH_ADJUSTMENT_DIVISOR }
    pub fn get_max_sensitivity_multiplier(&self) -> f64 { primitives::DEFAULT_MAX_SENSITIVITY_MULTIPLIER }
    pub fn get_long_channel_threshold(&self) -> f64 { *self.geometry.long_horizontal_threshold.get_raw_value() }
    pub fn get_long_channel_reduction_multiplier(&self) -> f64 { primitives::DEFAULT_LONG_CHANNEL_REDUCTION_MULTIPLIER }
    pub fn get_max_reduction_limit(&self) -> f64 { primitives::DEFAULT_MAX_REDUCTION_LIMIT }

    // --- Optimization ---
    pub fn get_max_optimization_iterations(&self) -> usize { *self.optimization.max_optimization_iterations.get_raw_value() }
    pub fn get_optimization_tolerance(&self) -> f64 { *self.optimization.convergence_tolerance.get_raw_value() }
    pub fn get_fast_wavelength_factors(&self) -> Vec<f64> { self.optimization.fast_wavelength_factors.get_raw_value().clone() }
    pub fn get_fast_wave_density_factors(&self) -> Vec<f64> { self.optimization.fast_wave_density_factors.get_raw_value().clone() }
    pub fn get_fast_fill_factors(&self) -> Vec<f64> { self.optimization.fast_fill_factors.get_raw_value().clone() }

    // --- Wave Generation ---
    pub fn get_smooth_endpoint_start_threshold(&self) -> f64 { *self.waves.smooth_endpoint_start_threshold.get_raw_value() }
    pub fn get_smooth_endpoint_end_threshold(&self) -> f64 { *self.waves.smooth_endpoint_end_threshold.get_raw_value() }
    pub fn get_default_transition_length_factor(&self) -> f64 { *self.waves.default_transition_length_factor.get_raw_value() }
    pub fn get_default_transition_amplitude_factor(&self) -> f64 { *self.waves.default_transition_amplitude_factor.get_raw_value() }
    pub fn get_default_transition_smoothness(&self) -> usize { *self.waves.default_transition_smoothness.get_raw_value() }
    pub fn get_default_wave_multiplier(&self) -> f64 { *self.waves.default_wave_multiplier.get_raw_value() }
    pub fn get_square_wave_sharpness(&self) -> f64 { *self.waves.square_wave_sharpness.get_raw_value() }
    pub fn get_neighbor_scale_factor(&self) -> f64 { *self.waves.neighbor_avoidance_scaling_factor.get_raw_value() }
    pub fn get_transition_zone_factor(&self) -> f64 { *self.waves.transition_zone_factor.get_raw_value() }

    // --- Geometry ---
    pub fn get_short_channel_width_multiplier(&self) -> f64 { *self.geometry.short_channel_width_multiplier.get_raw_value() }
    pub fn get_geometric_tolerance(&self) -> f64 { *self.waves.geometric_tolerance.get_raw_value() } // Note: was in waves constants
    pub fn get_max_curvature_reduction_factor(&self) -> f64 { *self.geometry.max_curvature_reduction_factor.get_raw_value() }
    pub fn get_min_curvature_factor(&self) -> f64 { *self.geometry.min_curvature_factor.get_raw_value() }
    pub fn get_min_distance_threshold(&self) -> f64 { *self.waves.geometric_tolerance.get_raw_value() } 

    // --- Strategies ---
    pub fn get_long_horizontal_threshold(&self) -> f64 { *self.geometry.long_horizontal_threshold.get_raw_value() }
    pub fn get_horizontal_angle_threshold(&self) -> f64 { *self.geometry.horizontal_angle_threshold.get_raw_value() }
    pub fn get_frustum_min_length_threshold(&self) -> f64 { *self.strategies.frustum_min_length_threshold.get_raw_value() }
    pub fn get_frustum_max_length_threshold(&self) -> f64 { *self.strategies.frustum_max_length_threshold.get_raw_value() }
    pub fn get_frustum_angle_threshold(&self) -> f64 { *self.strategies.frustum_angle_threshold.get_raw_value() }
    pub fn get_min_arc_length_threshold(&self) -> f64 { *self.geometry.min_arc_length_threshold.get_raw_value() }
    
    // --- Visualization ---
    pub fn get_default_chart_margin(&self) -> u32 { *self.visualization.default_chart_margin.get_raw_value() }
    pub fn get_default_chart_right_margin(&self) -> u32 { *self.visualization.default_chart_right_margin.get_raw_value() }
    pub fn get_default_x_label_area_size(&self) -> u32 { *self.visualization.default_x_label_area_size.get_raw_value() }
    pub fn get_default_y_label_area_size(&self) -> u32 { *self.visualization.default_y_label_area_size.get_raw_value() }
}

pub use primitives::*;
