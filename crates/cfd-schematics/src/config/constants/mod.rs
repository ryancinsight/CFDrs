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

use std::sync::OnceLock;

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
pub mod primitives;

/// Registry for all configuration constant groups.
///
/// A plain aggregate of compile-time constants: every group is `Copy`, so the
/// whole registry is `Copy` and can be constructed for free.
#[derive(Debug, Clone, Copy, PartialEq)]
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
    ///
    /// The registry is a plain aggregate of compile-time constants, so this is
    /// a `const fn` and costs nothing at run time.
    pub const fn new() -> Self {
        Self {
            strategies: StrategyThresholds::DEFAULT,
            waves: WaveGenerationConstants::DEFAULT,
            geometry: GeometryGenerationConstants::DEFAULT,
            optimization: OptimizationConstants::DEFAULT,
            visualization: VisualizationConstants::DEFAULT,
        }
    }

    /// Returns the process-wide canonical registry.
    ///
    /// The registry is immutable once constructed and every accessor takes
    /// `&self`, so a single shared instance is behaviourally identical to a
    /// per-call [`Self::new`]. Since [`Self::new`] is now a `const fn`, this
    /// accessor exists only to hand out a stable `&'static` handle without
    /// copying the aggregate.
    pub fn shared() -> &'static Self {
        static SHARED: OnceLock<ConstantsRegistry> = OnceLock::new();
        SHARED.get_or_init(Self::new)
    }

    // --- Adaptive Collision ---
    /// Returns the default maximum collision-adjustment factor.
    pub fn get_max_adjustment_factor(&self) -> f64 {
        primitives::DEFAULT_MAX_ADJUSTMENT_FACTOR
    }
    /// Returns the default minimum channel separation distance.
    pub fn get_min_channel_distance(&self) -> f64 {
        primitives::DEFAULT_MIN_CHANNEL_DISTANCE
    }
    /// Returns the configured minimum wall clearance.
    pub fn get_min_wall_distance(&self) -> f64 {
        self.geometry.default_wall_clearance
    }
    /// Returns the default collision safety-margin factor.
    pub fn get_safety_margin_factor(&self) -> f64 {
        primitives::DEFAULT_SAFETY_MARGIN_FACTOR
    }
    /// Returns the default maximum curvature-reduction factor.
    pub fn get_max_reduction_factor(&self) -> f64 {
        primitives::DEFAULT_MAX_REDUCTION_FACTOR
    }
    /// Returns the default collision-detection sensitivity.
    pub fn get_detection_sensitivity(&self) -> f64 {
        primitives::DEFAULT_DETECTION_SENSITIVITY
    }
    /// Returns the distance normalization divisor used by proximity checks.
    pub fn get_proximity_divisor(&self) -> f64 {
        primitives::DEFAULT_PROXIMITY_DIVISOR
    }
    /// Returns the minimum proximity adjustment factor.
    pub fn get_min_proximity_factor(&self) -> f64 {
        primitives::DEFAULT_MIN_PROXIMITY_FACTOR
    }
    /// Returns the maximum proximity adjustment factor.
    pub fn get_max_proximity_factor(&self) -> f64 {
        primitives::DEFAULT_MAX_PROXIMITY_FACTOR
    }
    /// Returns the configured branch-factor exponent.
    pub fn get_branch_factor_exponent(&self) -> f64 {
        self.optimization.branch_factor_exponent
    }
    /// Returns the divisor used to scale branch adjustments.
    pub fn get_branch_adjustment_divisor(&self) -> f64 {
        primitives::DEFAULT_BRANCH_ADJUSTMENT_DIVISOR
    }
    /// Returns the maximum sensitivity multiplier.
    pub fn get_max_sensitivity_multiplier(&self) -> f64 {
        primitives::DEFAULT_MAX_SENSITIVITY_MULTIPLIER
    }
    /// Returns the threshold for classifying a route as long.
    pub fn get_long_channel_threshold(&self) -> f64 {
        self.geometry.long_horizontal_threshold
    }
    /// Returns the default multiplier for long-channel reduction.
    pub fn get_long_channel_reduction_multiplier(&self) -> f64 {
        primitives::DEFAULT_LONG_CHANNEL_REDUCTION_MULTIPLIER
    }
    /// Returns the configured maximum reduction limit.
    pub fn get_max_reduction_limit(&self) -> f64 {
        primitives::DEFAULT_MAX_REDUCTION_LIMIT
    }

    // --- Optimization ---
    /// Returns the maximum number of optimization iterations.
    pub fn get_max_optimization_iterations(&self) -> usize {
        self.optimization.max_optimization_iterations
    }
    /// Returns the convergence tolerance used by optimization.
    pub fn get_optimization_tolerance(&self) -> f64 {
        self.optimization.convergence_tolerance
    }
    /// Returns the configured fast-path wavelength factors.
    ///
    /// Borrowed rather than cloned: the canonical factor tables are fixed for
    /// the lifetime of the registry, so callers can iterate the slice in place.
    pub fn get_fast_wavelength_factors(&self) -> &[f64] {
        self.optimization.fast_wavelength_factors
    }
    /// Returns the configured fast-path wave-density factors.
    pub fn get_fast_wave_density_factors(&self) -> &[f64] {
        self.optimization.fast_wave_density_factors
    }
    /// Returns the configured fast-path fill factors.
    pub fn get_fast_fill_factors(&self) -> &[f64] {
        self.optimization.fast_fill_factors
    }

    // --- Wave Generation ---
    /// Returns the lower threshold for smoothing a wave endpoint.
    pub fn get_smooth_endpoint_start_threshold(&self) -> f64 {
        self.waves.smooth_endpoint_start_threshold
    }
    /// Returns the upper threshold for smoothing a wave endpoint.
    pub fn get_smooth_endpoint_end_threshold(&self) -> f64 {
        self.waves.smooth_endpoint_end_threshold
    }
    /// Returns the default transition-length factor.
    pub fn get_default_transition_length_factor(&self) -> f64 {
        self.waves.default_transition_length_factor
    }
    /// Returns the default transition-amplitude factor.
    pub fn get_default_transition_amplitude_factor(&self) -> f64 {
        self.waves.default_transition_amplitude_factor
    }
    /// Returns the default transition smoothness.
    pub fn get_default_transition_smoothness(&self) -> usize {
        self.waves.default_transition_smoothness
    }
    /// Returns the default wave multiplier.
    pub fn get_default_wave_multiplier(&self) -> f64 {
        self.waves.default_wave_multiplier
    }
    /// Returns the sharpness used for square-wave generation.
    pub fn get_square_wave_sharpness(&self) -> f64 {
        self.waves.square_wave_sharpness
    }
    /// Returns the neighboring-channel avoidance scale factor.
    pub fn get_neighbor_scale_factor(&self) -> f64 {
        self.waves.neighbor_avoidance_scaling_factor
    }
    /// Returns the transition-zone scale factor.
    pub fn get_transition_zone_factor(&self) -> f64 {
        self.waves.transition_zone_factor
    }

    // --- Geometry ---
    /// Returns the multiplier for short-channel widths.
    pub fn get_short_channel_width_multiplier(&self) -> f64 {
        self.geometry.short_channel_width_multiplier
    }
    /// Returns the configured geometric tolerance.
    pub fn get_geometric_tolerance(&self) -> f64 {
        self.waves.geometric_tolerance
    }
    /// Returns the maximum curvature-reduction factor.
    pub fn get_max_curvature_reduction_factor(&self) -> f64 {
        self.geometry.max_curvature_reduction_factor
    }
    /// Returns the minimum curvature factor.
    pub fn get_min_curvature_factor(&self) -> f64 {
        self.geometry.min_curvature_factor
    }
    /// Returns the minimum-distance threshold used by geometry checks.
    pub fn get_min_distance_threshold(&self) -> f64 {
        self.waves.geometric_tolerance
    }

    // --- Strategies ---
    /// Returns the horizontal-route length threshold.
    pub fn get_long_horizontal_threshold(&self) -> f64 {
        self.geometry.long_horizontal_threshold
    }
    /// Returns the horizontal-angle threshold.
    pub fn get_horizontal_angle_threshold(&self) -> f64 {
        self.geometry.horizontal_angle_threshold
    }
    /// Returns the minimum frustum length threshold.
    pub fn get_frustum_min_length_threshold(&self) -> f64 {
        self.strategies.frustum_min_length_threshold
    }
    /// Returns the maximum frustum length threshold.
    pub fn get_frustum_max_length_threshold(&self) -> f64 {
        self.strategies.frustum_max_length_threshold
    }
    /// Returns the frustum angle threshold.
    pub fn get_frustum_angle_threshold(&self) -> f64 {
        self.strategies.frustum_angle_threshold
    }
    /// Returns the minimum arc-length threshold.
    pub fn get_min_arc_length_threshold(&self) -> f64 {
        self.geometry.min_arc_length_threshold
    }

    // --- Visualization ---
    /// Returns the default chart margin.
    pub fn get_default_chart_margin(&self) -> u32 {
        self.visualization.default_chart_margin
    }
    /// Returns the default chart right margin.
    pub fn get_default_chart_right_margin(&self) -> u32 {
        self.visualization.default_chart_right_margin
    }
    /// Returns the default x-axis label area size.
    pub fn get_default_x_label_area_size(&self) -> u32 {
        self.visualization.default_x_label_area_size
    }
    /// Returns the default y-axis label area size.
    pub fn get_default_y_label_area_size(&self) -> u32 {
        self.visualization.default_y_label_area_size
    }
}

impl Default for ConstantsRegistry {
    fn default() -> Self {
        Self::new()
    }
}

pub use primitives::*;
