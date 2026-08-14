//! Centralized configuration constants
//!
//! This module extracts all hardcoded values and magic numbers from throughout
//! the codebase into configurable parameters with proper validation and documentation.
//! This follows the SSOT (Single Source of Truth) principle and eliminates magic numbers.
//!
//! # Submodules
//!
//! - [`primitives`]: Raw `const` values (min/max/default bounds)
//! - `strategy`: [`StrategyThresholds`] for channel-type selection
//! - `wave`: [`WaveGenerationConstants`] for serpentine wave shaping
//! - `geometry`: [`GeometryGenerationConstants`] for point generation
//! - `optimization`: [`OptimizationConstants`] for solver tuning
//! - `visualization`: [`VisualizationConstants`] for chart rendering

mod geometry;
mod optimization;
mod strategy;
mod visualization;
mod wave;

pub use geometry::GeometryGenerationConstants;
pub use optimization::OptimizationConstants;
pub use strategy::StrategyThresholds;
pub use visualization::VisualizationConstants;
pub use wave::WaveGenerationConstants;

/// Configuration constants for geometry validation and defaults (Primitives)
pub mod primitives {
    /// Minimum allowed wall clearance (mm)
    pub const MIN_WALL_CLEARANCE: f64 = 0.01;
    /// Maximum allowed wall clearance (mm)
    pub const MAX_WALL_CLEARANCE: f64 = 10.0;
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

    /// Default wave multiplier for smooth transitions (2pi for one complete wave)
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
    /// This is a cosine-based ratio (`|dx| / distance`), bounded to [0, 1].
    /// A value of 0.75 corresponds to channels within ±41.4° of horizontal.
    pub const DEFAULT_HORIZONTAL_RATIO_THRESHOLD: f64 = 0.75;
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
    pub const MAX_CHANNEL_WIDTH: f64 = 50.0;
    /// Default channel width (mm)
    pub const DEFAULT_CHANNEL_WIDTH: f64 = 1.0;

    /// Minimum allowed channel height (mm)
    pub const MIN_CHANNEL_HEIGHT: f64 = 0.01;
    /// Maximum allowed channel height (mm)
    pub const MAX_CHANNEL_HEIGHT: f64 = 25.0;
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

    // Adaptive collision constants.
    /// Default maximum collision adjustment factor.
    pub const DEFAULT_MAX_ADJUSTMENT_FACTOR: f64 = 0.5;
    /// Default minimum distance between channels.
    pub const DEFAULT_MIN_CHANNEL_DISTANCE: f64 = 0.5;
    /// Default safety margin applied to collision checks.
    pub const DEFAULT_SAFETY_MARGIN_FACTOR: f64 = 1.1;
    /// Default maximum reduction factor for collision recovery.
    pub const DEFAULT_MAX_REDUCTION_FACTOR: f64 = 0.8;
    /// Default sensitivity used by collision detection.
    pub const DEFAULT_DETECTION_SENSITIVITY: f64 = 0.1;
    /// Default divisor for proximity normalization.
    pub const DEFAULT_PROXIMITY_DIVISOR: f64 = 10.0;
    /// Default minimum proximity factor.
    pub const DEFAULT_MIN_PROXIMITY_FACTOR: f64 = 0.2;
    /// Default maximum proximity factor.
    pub const DEFAULT_MAX_PROXIMITY_FACTOR: f64 = 2.0;
    /// Default divisor for branch adjustment scaling.
    pub const DEFAULT_BRANCH_ADJUSTMENT_DIVISOR: f64 = 5.0;
    /// Default maximum sensitivity multiplier.
    pub const DEFAULT_MAX_SENSITIVITY_MULTIPLIER: f64 = 3.0;
    /// Default multiplier for long-channel reduction.
    pub const DEFAULT_LONG_CHANNEL_REDUCTION_MULTIPLIER: f64 = 0.9;
    /// Default upper bound for collision reduction.
    pub const DEFAULT_MAX_REDUCTION_LIMIT: f64 = 0.5;
}

/// Registry for all validated configuration constant groups.
pub struct ConstantsRegistry {
    /// Thresholds used to select channel strategies.
    pub strategies: StrategyThresholds,
    /// Constants used by serpentine wave generation.
    pub waves: WaveGenerationConstants,
    /// Constants used by geometry generation.
    pub geometry: GeometryGenerationConstants,
    /// Constants used by geometry optimization.
    pub optimization: OptimizationConstants,
    /// Constants used by visualization rendering.
    pub visualization: VisualizationConstants,
}

impl ConstantsRegistry {
    /// Constructs a registry populated with the canonical defaults.
    #[must_use]
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
    /// Returns the default maximum collision-adjustment factor.
    #[must_use]
    pub fn get_max_adjustment_factor(&self) -> f64 {
        primitives::DEFAULT_MAX_ADJUSTMENT_FACTOR
    }
    /// Returns the default minimum channel separation distance.
    #[must_use]
    pub fn get_min_channel_distance(&self) -> f64 {
        primitives::DEFAULT_MIN_CHANNEL_DISTANCE
    }
    /// Returns the configured minimum wall clearance.
    #[must_use]
    pub fn get_min_wall_distance(&self) -> f64 {
        *self.geometry.default_wall_clearance.get_raw_value()
    }
    /// Returns the default collision safety-margin factor.
    #[must_use]
    pub fn get_safety_margin_factor(&self) -> f64 {
        primitives::DEFAULT_SAFETY_MARGIN_FACTOR
    }
    /// Returns the default maximum curvature-reduction factor.
    #[must_use]
    pub fn get_max_reduction_factor(&self) -> f64 {
        primitives::DEFAULT_MAX_REDUCTION_FACTOR
    }
    /// Returns the default collision-detection sensitivity.
    #[must_use]
    pub fn get_detection_sensitivity(&self) -> f64 {
        primitives::DEFAULT_DETECTION_SENSITIVITY
    }
    /// Returns the distance normalization divisor used by proximity checks.
    #[must_use]
    pub fn get_proximity_divisor(&self) -> f64 {
        primitives::DEFAULT_PROXIMITY_DIVISOR
    }
    /// Returns the minimum proximity adjustment factor.
    #[must_use]
    pub fn get_min_proximity_factor(&self) -> f64 {
        primitives::DEFAULT_MIN_PROXIMITY_FACTOR
    }
    /// Returns the maximum proximity adjustment factor.
    #[must_use]
    pub fn get_max_proximity_factor(&self) -> f64 {
        primitives::DEFAULT_MAX_PROXIMITY_FACTOR
    }
    /// Returns the configured branch-factor exponent.
    #[must_use]
    pub fn get_branch_factor_exponent(&self) -> f64 {
        *self.optimization.branch_factor_exponent.get_raw_value()
    }
    /// Returns the divisor used to scale branch adjustments.
    #[must_use]
    pub fn get_branch_adjustment_divisor(&self) -> f64 {
        primitives::DEFAULT_BRANCH_ADJUSTMENT_DIVISOR
    }
    /// Returns the maximum sensitivity multiplier.
    #[must_use]
    pub fn get_max_sensitivity_multiplier(&self) -> f64 {
        primitives::DEFAULT_MAX_SENSITIVITY_MULTIPLIER
    }
    /// Returns the threshold for classifying a route as long.
    #[must_use]
    pub fn get_long_channel_threshold(&self) -> f64 {
        *self.geometry.long_horizontal_threshold.get_raw_value()
    }
    /// Returns the default multiplier for long-channel reduction.
    #[must_use]
    pub fn get_long_channel_reduction_multiplier(&self) -> f64 {
        primitives::DEFAULT_LONG_CHANNEL_REDUCTION_MULTIPLIER
    }
    /// Returns the configured maximum reduction limit.
    #[must_use]
    pub fn get_max_reduction_limit(&self) -> f64 {
        primitives::DEFAULT_MAX_REDUCTION_LIMIT
    }

    // --- Optimization ---
    /// Returns the maximum number of optimization iterations.
    #[must_use]
    pub fn get_max_optimization_iterations(&self) -> usize {
        *self
            .optimization
            .max_optimization_iterations
            .get_raw_value()
    }
    /// Returns the convergence tolerance used by optimization.
    #[must_use]
    pub fn get_optimization_tolerance(&self) -> f64 {
        *self.optimization.convergence_tolerance.get_raw_value()
    }
    /// Returns the configured fast-path wavelength factors.
    #[must_use]
    pub fn get_fast_wavelength_factors(&self) -> Vec<f64> {
        self.optimization
            .fast_wavelength_factors
            .get_raw_value()
            .clone()
    }
    /// Returns the configured fast-path wave-density factors.
    #[must_use]
    pub fn get_fast_wave_density_factors(&self) -> Vec<f64> {
        self.optimization
            .fast_wave_density_factors
            .get_raw_value()
            .clone()
    }
    /// Returns the configured fast-path fill factors.
    #[must_use]
    pub fn get_fast_fill_factors(&self) -> Vec<f64> {
        self.optimization.fast_fill_factors.get_raw_value().clone()
    }

    // --- Wave Generation ---
    /// Returns the lower threshold for smoothing a wave endpoint.
    #[must_use]
    pub fn get_smooth_endpoint_start_threshold(&self) -> f64 {
        *self.waves.smooth_endpoint_start_threshold.get_raw_value()
    }
    /// Returns the upper threshold for smoothing a wave endpoint.
    #[must_use]
    pub fn get_smooth_endpoint_end_threshold(&self) -> f64 {
        *self.waves.smooth_endpoint_end_threshold.get_raw_value()
    }
    /// Returns the default transition-length factor.
    #[must_use]
    pub fn get_default_transition_length_factor(&self) -> f64 {
        *self.waves.default_transition_length_factor.get_raw_value()
    }
    /// Returns the default transition-amplitude factor.
    #[must_use]
    pub fn get_default_transition_amplitude_factor(&self) -> f64 {
        *self
            .waves
            .default_transition_amplitude_factor
            .get_raw_value()
    }
    /// Returns the default transition smoothness.
    #[must_use]
    pub fn get_default_transition_smoothness(&self) -> usize {
        *self.waves.default_transition_smoothness.get_raw_value()
    }
    /// Returns the default wave multiplier.
    #[must_use]
    pub fn get_default_wave_multiplier(&self) -> f64 {
        *self.waves.default_wave_multiplier.get_raw_value()
    }
    /// Returns the sharpness used for square-wave generation.
    #[must_use]
    pub fn get_square_wave_sharpness(&self) -> f64 {
        *self.waves.square_wave_sharpness.get_raw_value()
    }
    /// Returns the neighboring-channel avoidance scale factor.
    #[must_use]
    pub fn get_neighbor_scale_factor(&self) -> f64 {
        *self.waves.neighbor_avoidance_scaling_factor.get_raw_value()
    }
    /// Returns the transition-zone scale factor.
    #[must_use]
    pub fn get_transition_zone_factor(&self) -> f64 {
        *self.waves.transition_zone_factor.get_raw_value()
    }

    // --- Geometry ---
    /// Returns the multiplier for short-channel widths.
    #[must_use]
    pub fn get_short_channel_width_multiplier(&self) -> f64 {
        *self.geometry.short_channel_width_multiplier.get_raw_value()
    }
    /// Returns the configured geometric tolerance.
    #[must_use]
    pub fn get_geometric_tolerance(&self) -> f64 {
        *self.waves.geometric_tolerance.get_raw_value()
    }
    /// Returns the maximum curvature-reduction factor.
    #[must_use]
    pub fn get_max_curvature_reduction_factor(&self) -> f64 {
        *self.geometry.max_curvature_reduction_factor.get_raw_value()
    }
    /// Returns the minimum curvature factor.
    #[must_use]
    pub fn get_min_curvature_factor(&self) -> f64 {
        *self.geometry.min_curvature_factor.get_raw_value()
    }
    /// Returns the minimum-distance threshold used by geometry checks.
    #[must_use]
    pub fn get_min_distance_threshold(&self) -> f64 {
        *self.waves.geometric_tolerance.get_raw_value()
    }

    // --- Strategies ---
    /// Returns the horizontal-route length threshold.
    #[must_use]
    pub fn get_long_horizontal_threshold(&self) -> f64 {
        *self.geometry.long_horizontal_threshold.get_raw_value()
    }
    /// Returns the horizontal-angle threshold.
    #[must_use]
    pub fn get_horizontal_angle_threshold(&self) -> f64 {
        *self.geometry.horizontal_angle_threshold.get_raw_value()
    }
    /// Returns the minimum frustum length threshold.
    #[must_use]
    pub fn get_frustum_min_length_threshold(&self) -> f64 {
        *self.strategies.frustum_min_length_threshold.get_raw_value()
    }
    /// Returns the maximum frustum length threshold.
    #[must_use]
    pub fn get_frustum_max_length_threshold(&self) -> f64 {
        *self.strategies.frustum_max_length_threshold.get_raw_value()
    }
    /// Returns the frustum angle threshold.
    #[must_use]
    pub fn get_frustum_angle_threshold(&self) -> f64 {
        *self.strategies.frustum_angle_threshold.get_raw_value()
    }
    /// Returns the minimum arc-length threshold.
    #[must_use]
    pub fn get_min_arc_length_threshold(&self) -> f64 {
        *self.geometry.min_arc_length_threshold.get_raw_value()
    }

    // --- Visualization ---
    /// Returns the default chart margin.
    #[must_use]
    pub fn get_default_chart_margin(&self) -> u32 {
        *self.visualization.default_chart_margin.get_raw_value()
    }
    /// Returns the default chart right margin.
    #[must_use]
    pub fn get_default_chart_right_margin(&self) -> u32 {
        *self
            .visualization
            .default_chart_right_margin
            .get_raw_value()
    }
    /// Returns the default x-axis label area size.
    #[must_use]
    pub fn get_default_x_label_area_size(&self) -> u32 {
        *self.visualization.default_x_label_area_size.get_raw_value()
    }
    /// Returns the default y-axis label area size.
    #[must_use]
    pub fn get_default_y_label_area_size(&self) -> u32 {
        *self.visualization.default_y_label_area_size.get_raw_value()
    }
}

impl Default for ConstantsRegistry {
    fn default() -> Self {
        Self::new()
    }
}

pub use primitives::*;
