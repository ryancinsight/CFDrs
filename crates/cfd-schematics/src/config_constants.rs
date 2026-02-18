//! Centralized configuration constants
//!
//! This module extracts all hardcoded values and magic numbers from throughout
//! the codebase into configurable parameters with proper validation and documentation.
//! This follows the SSOT (Single Source of Truth) principle and eliminates magic numbers.

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

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
                40u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(20, 200),
                ]),
                ParameterMetadata::new(
                    "default_x_label_area_size",
                    "Default label area size for x-axis",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_y_label_area_size: ConfigurableParameter::new(
                40u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(20, 200),
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
                    ParameterConstraints::range(0.01, 1.0),
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

/// Collision detection constants previously hardcoded
pub struct CollisionConstants {
    /// Minimum channel distance for collision detection
    pub min_channel_distance: ConfigurableParameter<f64>,

    /// Minimum wall distance for collision detection
    pub min_wall_distance: ConfigurableParameter<f64>,

    /// Minimum wall thickness for manufacturing constraints
    pub min_wall_thickness: ConfigurableParameter<f64>,

    /// Safety margin factor for collision detection
    pub safety_margin_factor: ConfigurableParameter<f64>,

    /// Maximum reduction factor for collision avoidance
    pub max_reduction_factor: ConfigurableParameter<f64>,

    /// Detection sensitivity for collision detection
    pub detection_sensitivity: ConfigurableParameter<f64>,

    /// Neighbor scale factor for adaptive collision detection
    pub neighbor_scale_factor: ConfigurableParameter<f64>,

    /// Minimum distance threshold for adaptive collision detection
    pub min_distance_threshold: ConfigurableParameter<f64>,

    /// Maximum adjustment factor for adaptive collision detection
    pub max_adjustment_factor: ConfigurableParameter<f64>,

    /// Proximity divisor for collision detection
    pub proximity_divisor: ConfigurableParameter<f64>,

    /// Minimum proximity factor for collision detection
    pub min_proximity_factor: ConfigurableParameter<f64>,

    /// Maximum proximity factor for collision detection
    pub max_proximity_factor: ConfigurableParameter<f64>,

    /// Branch adjustment divisor for collision detection
    pub branch_adjustment_divisor: ConfigurableParameter<f64>,

    /// Maximum sensitivity multiplier for collision detection
    pub max_sensitivity_multiplier: ConfigurableParameter<f64>,

    /// Long channel threshold for collision detection
    pub long_channel_threshold: ConfigurableParameter<f64>,

    /// Long channel reduction multiplier for collision detection
    pub long_channel_reduction_multiplier: ConfigurableParameter<f64>,

    /// Maximum reduction limit for collision detection
    pub max_reduction_limit: ConfigurableParameter<f64>,
}

impl Default for VisualizationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}

impl CollisionConstants {
    /// Create default collision constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            min_channel_distance: ConfigurableParameter::new(
                2.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 10.0),
                ]),
                ParameterMetadata::new(
                    "min_channel_distance",
                    "Minimum distance between channels for collision detection",
                    "collision",
                )
                .with_units("mm"),
            ),

            min_wall_thickness: ConfigurableParameter::new(
                0.45,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 5.0),
                ]),
                ParameterMetadata::new(
                    "min_wall_thickness",
                    "Minimum wall thickness between channels for manufacturing constraints",
                    "collision",
                )
                .with_units("mm"),
            ),

            min_wall_distance: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 5.0),
                ]),
                ParameterMetadata::new(
                    "min_wall_distance",
                    "Minimum distance from walls for collision detection",
                    "collision",
                )
                .with_units("mm"),
            ),

            safety_margin_factor: ConfigurableParameter::new(
                1.2,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 3.0),
                ]),
                ParameterMetadata::new(
                    "safety_margin_factor",
                    "Safety margin multiplier for collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            max_reduction_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 1.0),
                ]),
                ParameterMetadata::new(
                    "max_reduction_factor",
                    "Maximum reduction factor for collision avoidance",
                    "collision",
                )
                .with_units("ratio"),
            ),

            detection_sensitivity: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 5.0),
                ]),
                ParameterMetadata::new(
                    "detection_sensitivity",
                    "Sensitivity factor for collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            neighbor_scale_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 2.0),
                ]),
                ParameterMetadata::new(
                    "neighbor_scale_factor",
                    "Scale factor for neighbor-based collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            min_distance_threshold: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 5.0),
                ]),
                ParameterMetadata::new(
                    "min_distance_threshold",
                    "Minimum distance threshold for adaptive collision detection",
                    "collision",
                )
                .with_units("mm"),
            ),

            max_adjustment_factor: ConfigurableParameter::new(
                2.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 5.0),
                ]),
                ParameterMetadata::new(
                    "max_adjustment_factor",
                    "Maximum adjustment factor for adaptive collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            proximity_divisor: ConfigurableParameter::new(
                10.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 50.0),
                ]),
                ParameterMetadata::new(
                    "proximity_divisor",
                    "Divisor for proximity calculations in collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            min_proximity_factor: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 1.0),
                ]),
                ParameterMetadata::new(
                    "min_proximity_factor",
                    "Minimum proximity factor for collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            max_proximity_factor: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.5, 2.0),
                ]),
                ParameterMetadata::new(
                    "max_proximity_factor",
                    "Maximum proximity factor for collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            branch_adjustment_divisor: ConfigurableParameter::new(
                4.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 20.0),
                ]),
                ParameterMetadata::new(
                    "branch_adjustment_divisor",
                    "Divisor for branch-based adjustments in collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            max_sensitivity_multiplier: ConfigurableParameter::new(
                2.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 10.0),
                ]),
                ParameterMetadata::new(
                    "max_sensitivity_multiplier",
                    "Maximum sensitivity multiplier for collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            long_channel_threshold: ConfigurableParameter::new(
                50.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(10.0, 200.0),
                ]),
                ParameterMetadata::new(
                    "long_channel_threshold",
                    "Threshold for considering a channel as long in collision detection",
                    "collision",
                )
                .with_units("mm"),
            ),

            long_channel_reduction_multiplier: ConfigurableParameter::new(
                1.2,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 3.0),
                ]),
                ParameterMetadata::new(
                    "long_channel_reduction_multiplier",
                    "Reduction multiplier for long channels in collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),

            max_reduction_limit: ConfigurableParameter::new(
                0.95,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.5, 1.0),
                ]),
                ParameterMetadata::new(
                    "max_reduction_limit",
                    "Maximum reduction limit for collision detection",
                    "collision",
                )
                .with_units("ratio"),
            ),
        }
    }
}

impl Default for CollisionConstants {
    fn default() -> Self {
        Self::make_default()
    }
}

/// Central constants registry containing all configuration constants
pub struct ConstantsRegistry {
    /// Strategy selection thresholds
    pub strategy_thresholds: StrategyThresholds,

    /// Wave generation constants
    pub wave_generation: WaveGenerationConstants,

    /// Geometry generation constants
    pub geometry_generation: GeometryGenerationConstants,

    /// Optimization constants
    pub optimization: OptimizationConstants,

    /// Visualization constants
    pub visualization: VisualizationConstants,

    /// Collision detection constants
    pub collision: CollisionConstants,
}

impl ConstantsRegistry {
    /// Create a new constants registry with default values
    #[must_use]
    pub fn new() -> Self {
        Self {
            strategy_thresholds: StrategyThresholds::default(),
            wave_generation: WaveGenerationConstants::default(),
            geometry_generation: GeometryGenerationConstants::default(),
            optimization: OptimizationConstants::default(),
            visualization: VisualizationConstants::default(),
            collision: CollisionConstants::default(),
        }
    }

    /// Validate all constants
    pub const fn validate_all(&self) -> crate::error::SchemeResult<()> {
        // Validate all parameter groups
        // This would be expanded to validate all parameters in each group
        Ok(())
    }
}

impl Default for ConstantsRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper functions to get commonly used constant values
impl ConstantsRegistry {
    /// Get the square wave sharpness factor
    #[must_use]
    pub const fn get_square_wave_sharpness(&self) -> f64 {
        *self.wave_generation.square_wave_sharpness.get_raw_value()
    }

    /// Get the transition zone factor
    #[must_use]
    pub const fn get_transition_zone_factor(&self) -> f64 {
        *self.wave_generation.transition_zone_factor.get_raw_value()
    }

    /// Get the arc curvature threshold
    #[must_use]
    pub const fn get_arc_curvature_threshold(&self) -> f64 {
        *self
            .strategy_thresholds
            .arc_curvature_threshold
            .get_raw_value()
    }

    /// Get the serpentine fill threshold
    #[must_use]
    pub const fn get_serpentine_fill_threshold(&self) -> f64 {
        *self
            .strategy_thresholds
            .serpentine_fill_threshold
            .get_raw_value()
    }

    /// Get the adaptive branch threshold
    #[must_use]
    pub const fn get_adaptive_branch_threshold(&self) -> usize {
        *self
            .strategy_thresholds
            .adaptive_branch_threshold
            .get_raw_value()
    }

    /// Get the branch factor exponent
    #[must_use]
    pub const fn get_branch_factor_exponent(&self) -> f64 {
        *self.optimization.branch_factor_exponent.get_raw_value()
    }

    /// Get the minimum channel distance for collision detection
    #[must_use]
    pub const fn get_min_channel_distance(&self) -> f64 {
        *self.collision.min_channel_distance.get_raw_value()
    }

    /// Get the minimum wall distance for collision detection
    #[must_use]
    pub const fn get_min_wall_distance(&self) -> f64 {
        *self.collision.min_wall_distance.get_raw_value()
    }

    /// Get the minimum wall thickness for manufacturing constraints
    #[must_use]
    pub const fn get_min_wall_thickness(&self) -> f64 {
        *self.collision.min_wall_thickness.get_raw_value()
    }

    /// Get the safety margin factor for collision detection
    #[must_use]
    pub const fn get_safety_margin_factor(&self) -> f64 {
        *self.collision.safety_margin_factor.get_raw_value()
    }

    /// Get the maximum reduction factor for collision avoidance
    #[must_use]
    pub const fn get_max_reduction_factor(&self) -> f64 {
        *self.collision.max_reduction_factor.get_raw_value()
    }

    /// Get the detection sensitivity for collision detection
    #[must_use]
    pub const fn get_detection_sensitivity(&self) -> f64 {
        *self.collision.detection_sensitivity.get_raw_value()
    }

    /// Get the neighbor scale factor for adaptive collision detection
    #[must_use]
    pub const fn get_neighbor_scale_factor(&self) -> f64 {
        *self.collision.neighbor_scale_factor.get_raw_value()
    }

    /// Get the minimum distance threshold for adaptive collision detection
    #[must_use]
    pub const fn get_min_distance_threshold(&self) -> f64 {
        *self.collision.min_distance_threshold.get_raw_value()
    }

    /// Get the maximum adjustment factor for adaptive collision detection
    #[must_use]
    pub const fn get_max_adjustment_factor(&self) -> f64 {
        *self.collision.max_adjustment_factor.get_raw_value()
    }

    /// Get the maximum number of optimization iterations
    #[must_use]
    pub const fn get_max_optimization_iterations(&self) -> usize {
        *self
            .optimization
            .max_optimization_iterations
            .get_raw_value()
    }

    /// Get the optimization tolerance
    #[must_use]
    pub const fn get_optimization_tolerance(&self) -> f64 {
        *self.optimization.convergence_tolerance.get_raw_value()
    }

    /// Get the proximity divisor for collision detection
    #[must_use]
    pub const fn get_proximity_divisor(&self) -> f64 {
        *self.collision.proximity_divisor.get_raw_value()
    }

    /// Get the minimum proximity factor for collision detection
    #[must_use]
    pub const fn get_min_proximity_factor(&self) -> f64 {
        *self.collision.min_proximity_factor.get_raw_value()
    }

    /// Get the maximum proximity factor for collision detection
    #[must_use]
    pub const fn get_max_proximity_factor(&self) -> f64 {
        *self.collision.max_proximity_factor.get_raw_value()
    }

    /// Get the branch adjustment divisor for collision detection
    #[must_use]
    pub const fn get_branch_adjustment_divisor(&self) -> f64 {
        *self.collision.branch_adjustment_divisor.get_raw_value()
    }

    /// Get the maximum sensitivity multiplier for collision detection
    #[must_use]
    pub const fn get_max_sensitivity_multiplier(&self) -> f64 {
        *self.collision.max_sensitivity_multiplier.get_raw_value()
    }

    /// Get the long channel threshold for collision detection
    #[must_use]
    pub const fn get_long_channel_threshold(&self) -> f64 {
        *self.collision.long_channel_threshold.get_raw_value()
    }

    /// Get the long channel reduction multiplier for collision detection
    #[must_use]
    pub const fn get_long_channel_reduction_multiplier(&self) -> f64 {
        *self
            .collision
            .long_channel_reduction_multiplier
            .get_raw_value()
    }

    /// Get the maximum reduction limit for collision detection
    #[must_use]
    pub const fn get_max_reduction_limit(&self) -> f64 {
        *self.collision.max_reduction_limit.get_raw_value()
    }

    /// Get the fill factor enhancement multiplier
    #[must_use]
    pub const fn get_fill_factor_enhancement(&self) -> f64 {
        *self.optimization.fill_factor_enhancement.get_raw_value()
    }

    /// Get the default serpentine points
    #[must_use]
    pub const fn get_default_serpentine_points(&self) -> usize {
        *self
            .geometry_generation
            .default_serpentine_points
            .get_raw_value()
    }

    /// Get the default wall clearance
    #[must_use]
    pub const fn get_default_wall_clearance(&self) -> f64 {
        *self
            .geometry_generation
            .default_wall_clearance
            .get_raw_value()
    }

    /// Get the default channel width
    #[must_use]
    pub const fn get_default_channel_width(&self) -> f64 {
        *self
            .geometry_generation
            .default_channel_width
            .get_raw_value()
    }

    /// Get the fast optimization wavelength factors
    #[must_use]
    pub const fn get_fast_wavelength_factors(&self) -> &Vec<f64> {
        self.optimization.fast_wavelength_factors.get_raw_value()
    }

    /// Get the fast optimization wave density factors
    #[must_use]
    pub const fn get_fast_wave_density_factors(&self) -> &Vec<f64> {
        self.optimization.fast_wave_density_factors.get_raw_value()
    }

    /// Get the fast optimization fill factors
    #[must_use]
    pub const fn get_fast_fill_factors(&self) -> &Vec<f64> {
        self.optimization.fast_fill_factors.get_raw_value()
    }

    /// Get the smooth endpoint start threshold
    #[must_use]
    pub const fn get_smooth_endpoint_start_threshold(&self) -> f64 {
        *self
            .wave_generation
            .smooth_endpoint_start_threshold
            .get_raw_value()
    }

    /// Get the smooth endpoint end threshold
    #[must_use]
    pub const fn get_smooth_endpoint_end_threshold(&self) -> f64 {
        *self
            .wave_generation
            .smooth_endpoint_end_threshold
            .get_raw_value()
    }

    /// Get the default transition length factor
    #[must_use]
    pub const fn get_default_transition_length_factor(&self) -> f64 {
        *self
            .wave_generation
            .default_transition_length_factor
            .get_raw_value()
    }

    /// Get the default transition amplitude factor
    #[must_use]
    pub const fn get_default_transition_amplitude_factor(&self) -> f64 {
        *self
            .wave_generation
            .default_transition_amplitude_factor
            .get_raw_value()
    }

    /// Get the default transition smoothness
    #[must_use]
    pub const fn get_default_transition_smoothness(&self) -> usize {
        *self
            .wave_generation
            .default_transition_smoothness
            .get_raw_value()
    }

    /// Get the default wave multiplier
    #[must_use]
    pub const fn get_default_wave_multiplier(&self) -> f64 {
        *self.wave_generation.default_wave_multiplier.get_raw_value()
    }

    /// Get the wall proximity scaling factor
    #[must_use]
    pub const fn get_wall_proximity_scaling_factor(&self) -> f64 {
        *self
            .wave_generation
            .wall_proximity_scaling_factor
            .get_raw_value()
    }

    /// Get the neighbor avoidance scaling factor
    #[must_use]
    pub const fn get_neighbor_avoidance_scaling_factor(&self) -> f64 {
        *self
            .wave_generation
            .neighbor_avoidance_scaling_factor
            .get_raw_value()
    }

    /// Get the geometric tolerance
    #[must_use]
    pub const fn get_geometric_tolerance(&self) -> f64 {
        *self.wave_generation.geometric_tolerance.get_raw_value()
    }

    /// Get the short channel width multiplier
    #[must_use]
    pub const fn get_short_channel_width_multiplier(&self) -> f64 {
        *self
            .geometry_generation
            .short_channel_width_multiplier
            .get_raw_value()
    }

    /// Get the smooth straight middle points
    #[must_use]
    pub const fn get_smooth_straight_middle_points(&self) -> usize {
        *self
            .geometry_generation
            .smooth_straight_middle_points
            .get_raw_value()
    }

    /// Get the horizontal angle threshold
    #[must_use]
    pub const fn get_horizontal_angle_threshold(&self) -> f64 {
        *self
            .geometry_generation
            .horizontal_angle_threshold
            .get_raw_value()
    }

    /// Get the long horizontal threshold
    #[must_use]
    pub const fn get_long_horizontal_threshold(&self) -> f64 {
        *self
            .geometry_generation
            .long_horizontal_threshold
            .get_raw_value()
    }

    /// Get the minimum arc length threshold
    #[must_use]
    pub const fn get_min_arc_length_threshold(&self) -> f64 {
        *self
            .geometry_generation
            .min_arc_length_threshold
            .get_raw_value()
    }

    /// Get the maximum curvature reduction factor
    #[must_use]
    pub const fn get_max_curvature_reduction_factor(&self) -> f64 {
        *self
            .geometry_generation
            .max_curvature_reduction_factor
            .get_raw_value()
    }

    /// Get the minimum curvature factor
    #[must_use]
    pub const fn get_min_curvature_factor(&self) -> f64 {
        *self
            .geometry_generation
            .min_curvature_factor
            .get_raw_value()
    }

    /// Get the default chart margin
    #[must_use]
    pub const fn get_default_chart_margin(&self) -> u32 {
        *self.visualization.default_chart_margin.get_raw_value()
    }

    /// Get the default chart right margin
    #[must_use]
    pub const fn get_default_chart_right_margin(&self) -> u32 {
        *self
            .visualization
            .default_chart_right_margin
            .get_raw_value()
    }

    /// Get the default x-label area size
    #[must_use]
    pub const fn get_default_x_label_area_size(&self) -> u32 {
        *self.visualization.default_x_label_area_size.get_raw_value()
    }

    /// Get the default y-label area size
    #[must_use]
    pub const fn get_default_y_label_area_size(&self) -> u32 {
        *self.visualization.default_y_label_area_size.get_raw_value()
    }

    /// Get the default boundary buffer factor
    #[must_use]
    pub const fn get_default_boundary_buffer_factor(&self) -> f64 {
        *self
            .visualization
            .default_boundary_buffer_factor
            .get_raw_value()
    }

    /// Get the frustum minimum length threshold
    #[must_use]
    pub const fn get_frustum_min_length_threshold(&self) -> f64 {
        *self
            .strategy_thresholds
            .frustum_min_length_threshold
            .get_raw_value()
    }

    /// Get the frustum maximum length threshold
    #[must_use]
    pub const fn get_frustum_max_length_threshold(&self) -> f64 {
        *self
            .strategy_thresholds
            .frustum_max_length_threshold
            .get_raw_value()
    }

    /// Get the frustum angle threshold
    #[must_use]
    pub const fn get_frustum_angle_threshold(&self) -> f64 {
        *self
            .strategy_thresholds
            .frustum_angle_threshold
            .get_raw_value()
    }
}
