//! config.rs - 2D Schematic Configuration
//!
//! This module provides comprehensive configuration management for 2D microfluidic
//! schematic generation. It centralizes all configuration logic and eliminates
//! Single Source of Truth (SSOT) violations by providing validated configuration
//! types with clear constraints and relationships.
//!
//! # Design Principles
//!
//! - **Single Source of Truth**: All configuration values are defined once
//! - **Validation**: All configurations are validated at creation time
//! - **Immutability**: Configurations are immutable after creation
//! - **Composability**: Complex configurations are built from simpler ones
//! - **Discoverability**: Presets and builders make common configurations easy

use crate::error::{ConfigurationError, ConfigurationResult};
use crate::geometry::strategies::SmoothTransitionConfig;
use crate::geometry::ChannelType;
use serde::{Deserialize, Serialize};

// Re-export constants module
pub use crate::config_constants::*;

/// Configuration constants for geometry validation and defaults
pub mod constants {
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
}

/// Configuration for adaptive serpentine channel behavior
///
/// This configuration controls how serpentine channels adapt their wave properties
/// based on geometric constraints such as distance from nodes, walls, and other channels.
///
/// # Examples
///
/// ```rust
/// use scheme::config::AdaptiveSerpentineConfig;
///
/// // Create with default adaptive behavior
/// let config = AdaptiveSerpentineConfig::default();
///
/// // Create with custom adaptive parameters
/// let custom_config = AdaptiveSerpentineConfig {
///     node_distance_normalization: 15.0,
///     plateau_width_factor: 0.5,
///     horizontal_ratio_threshold: 0.75,
///     middle_section_amplitude_factor: 0.4,
///     plateau_amplitude_factor: 0.9,
///     enable_distance_based_scaling: true,
///     enable_wall_proximity_scaling: true,
///     enable_neighbor_avoidance: true,
/// };
/// ```
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveSerpentineConfig {
    /// Distance normalization factor for node proximity effects (1.0-50.0)
    pub node_distance_normalization: f64,
    /// Plateau width factor for horizontal channels (0.1-0.8, fraction of channel length)
    pub plateau_width_factor: f64,
    /// Horizontal ratio threshold for middle section detection (0.5-0.95)
    pub horizontal_ratio_threshold: f64,
    /// Amplitude factor for middle sections of horizontal channels (0.1-1.0)
    pub middle_section_amplitude_factor: f64,
    /// Amplitude factor for plateau regions (0.5-1.0)
    pub plateau_amplitude_factor: f64,
    /// Enable distance-based amplitude scaling near nodes
    pub enable_distance_based_scaling: bool,
    /// Enable amplitude scaling based on wall proximity
    pub enable_wall_proximity_scaling: bool,
    /// Enable amplitude reduction for neighbor channel avoidance
    pub enable_neighbor_avoidance: bool,
}

impl Default for AdaptiveSerpentineConfig {
    fn default() -> Self {
        Self {
            node_distance_normalization: constants::DEFAULT_NODE_DISTANCE_NORMALIZATION,
            plateau_width_factor: constants::DEFAULT_PLATEAU_WIDTH_FACTOR,
            horizontal_ratio_threshold: constants::DEFAULT_HORIZONTAL_RATIO_THRESHOLD,
            middle_section_amplitude_factor: constants::DEFAULT_MIDDLE_SECTION_AMPLITUDE_FACTOR,
            plateau_amplitude_factor: constants::DEFAULT_PLATEAU_AMPLITUDE_FACTOR,
            enable_distance_based_scaling: true,
            enable_wall_proximity_scaling: true,
            enable_neighbor_avoidance: true,
        }
    }
}

impl AdaptiveSerpentineConfig {
    /// Create a new adaptive serpentine configuration with validation
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        node_distance_normalization: f64,
        plateau_width_factor: f64,
        horizontal_ratio_threshold: f64,
        middle_section_amplitude_factor: f64,
        plateau_amplitude_factor: f64,
        enable_distance_based_scaling: bool,
        enable_wall_proximity_scaling: bool,
        enable_neighbor_avoidance: bool,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            node_distance_normalization,
            plateau_width_factor,
            horizontal_ratio_threshold,
            middle_section_amplitude_factor,
            plateau_amplitude_factor,
            enable_distance_based_scaling,
            enable_wall_proximity_scaling,
            enable_neighbor_avoidance,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the configuration parameters
    pub fn validate(&self) -> ConfigurationResult<()> {
        if self.node_distance_normalization < constants::MIN_NODE_DISTANCE_NORMALIZATION
            || self.node_distance_normalization > constants::MAX_NODE_DISTANCE_NORMALIZATION
        {
            return Err(ConfigurationError::invalid_generation_config(
                "node_distance_normalization",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_NODE_DISTANCE_NORMALIZATION,
                    constants::MAX_NODE_DISTANCE_NORMALIZATION
                ),
            ));
        }

        if self.plateau_width_factor < constants::MIN_PLATEAU_WIDTH_FACTOR
            || self.plateau_width_factor > constants::MAX_PLATEAU_WIDTH_FACTOR
        {
            return Err(ConfigurationError::invalid_generation_config(
                "plateau_width_factor",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_PLATEAU_WIDTH_FACTOR,
                    constants::MAX_PLATEAU_WIDTH_FACTOR
                ),
            ));
        }

        if self.horizontal_ratio_threshold < constants::MIN_HORIZONTAL_RATIO_THRESHOLD
            || self.horizontal_ratio_threshold > constants::MAX_HORIZONTAL_RATIO_THRESHOLD
        {
            return Err(ConfigurationError::invalid_generation_config(
                "horizontal_ratio_threshold",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_HORIZONTAL_RATIO_THRESHOLD,
                    constants::MAX_HORIZONTAL_RATIO_THRESHOLD
                ),
            ));
        }

        if self.middle_section_amplitude_factor < constants::MIN_MIDDLE_SECTION_AMPLITUDE_FACTOR
            || self.middle_section_amplitude_factor > constants::MAX_MIDDLE_SECTION_AMPLITUDE_FACTOR
        {
            return Err(ConfigurationError::invalid_generation_config(
                "middle_section_amplitude_factor",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_MIDDLE_SECTION_AMPLITUDE_FACTOR,
                    constants::MAX_MIDDLE_SECTION_AMPLITUDE_FACTOR
                ),
            ));
        }

        if self.plateau_amplitude_factor < constants::MIN_PLATEAU_AMPLITUDE_FACTOR
            || self.plateau_amplitude_factor > constants::MAX_PLATEAU_AMPLITUDE_FACTOR
        {
            return Err(ConfigurationError::invalid_generation_config(
                "plateau_amplitude_factor",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_PLATEAU_AMPLITUDE_FACTOR,
                    constants::MAX_PLATEAU_AMPLITUDE_FACTOR
                ),
            ));
        }

        Ok(())
    }

    /// Create a conservative configuration with minimal adaptive behavior
    #[must_use]
    pub const fn conservative() -> Self {
        Self {
            node_distance_normalization: 20.0,
            plateau_width_factor: 0.2,
            horizontal_ratio_threshold: 0.9,
            middle_section_amplitude_factor: 0.2,
            plateau_amplitude_factor: 0.6,
            enable_distance_based_scaling: true,
            enable_wall_proximity_scaling: true,
            enable_neighbor_avoidance: true,
        }
    }

    /// Create an aggressive configuration with strong adaptive behavior
    #[must_use]
    pub const fn aggressive() -> Self {
        Self {
            node_distance_normalization: 5.0,
            plateau_width_factor: 0.6,
            horizontal_ratio_threshold: 0.6,
            middle_section_amplitude_factor: 0.6,
            plateau_amplitude_factor: 0.95,
            enable_distance_based_scaling: true,
            enable_wall_proximity_scaling: true,
            enable_neighbor_avoidance: true,
        }
    }

    /// Create a configuration with adaptive behavior disabled
    #[must_use]
    pub fn disabled() -> Self {
        Self {
            enable_distance_based_scaling: false,
            enable_wall_proximity_scaling: false,
            enable_neighbor_avoidance: false,
            ..Self::default()
        }
    }
}

/// Configuration for geometry generation parameters
///
/// This struct controls the quality and precision of geometry generation,
/// including point density for path generation and wave calculations.
///
/// # Examples
///
/// ```rust
/// use scheme::config::GeometryGenerationConfig;
///
/// // Create with default values
/// let config = GeometryGenerationConfig::default();
///
/// // Create with custom values for high-quality generation
/// let high_quality = GeometryGenerationConfig {
///     serpentine_points: 200,
///     optimization_points: 100,
///     smooth_straight_middle_points: 20,
///     transition_wave_multiplier: 2.0,
/// };
/// ```
#[derive(Clone, Copy, Debug)]
pub struct GeometryGenerationConfig {
    /// Number of points to generate for serpentine paths (10-1000)
    pub serpentine_points: usize,
    /// Number of points to generate for optimization paths (10-200)
    pub optimization_points: usize,
    /// Number of middle points for smooth straight channels (2-50)
    pub smooth_straight_middle_points: usize,
    /// Wave multiplier for smooth transitions (0.5-10.0, where 2.0 = one complete wave)
    pub transition_wave_multiplier: f64,
}

impl Default for GeometryGenerationConfig {
    fn default() -> Self {
        Self {
            serpentine_points: constants::DEFAULT_SERPENTINE_POINTS,
            optimization_points: constants::DEFAULT_OPTIMIZATION_POINTS,
            smooth_straight_middle_points: constants::DEFAULT_SMOOTH_STRAIGHT_MIDDLE_POINTS,
            transition_wave_multiplier: constants::DEFAULT_TRANSITION_WAVE_MULTIPLIER,
        }
    }
}

impl GeometryGenerationConfig {
    /// Create a new geometry generation configuration with validation
    pub fn new(
        serpentine_points: usize,
        optimization_points: usize,
        smooth_straight_middle_points: usize,
        transition_wave_multiplier: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            serpentine_points,
            optimization_points,
            smooth_straight_middle_points,
            transition_wave_multiplier,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the configuration parameters
    pub fn validate(&self) -> ConfigurationResult<()> {
        if self.serpentine_points < constants::MIN_SERPENTINE_POINTS
            || self.serpentine_points > constants::MAX_SERPENTINE_POINTS
        {
            return Err(ConfigurationError::invalid_generation_config(
                "serpentine_points",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_SERPENTINE_POINTS,
                    constants::MAX_SERPENTINE_POINTS
                ),
            ));
        }

        if self.optimization_points < constants::MIN_OPTIMIZATION_POINTS
            || self.optimization_points > constants::MAX_OPTIMIZATION_POINTS
        {
            return Err(ConfigurationError::invalid_generation_config(
                "optimization_points",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_OPTIMIZATION_POINTS,
                    constants::MAX_OPTIMIZATION_POINTS
                ),
            ));
        }

        if self.smooth_straight_middle_points < constants::MIN_SMOOTH_STRAIGHT_MIDDLE_POINTS
            || self.smooth_straight_middle_points > constants::MAX_SMOOTH_STRAIGHT_MIDDLE_POINTS
        {
            return Err(ConfigurationError::invalid_generation_config(
                "smooth_straight_middle_points",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_SMOOTH_STRAIGHT_MIDDLE_POINTS,
                    constants::MAX_SMOOTH_STRAIGHT_MIDDLE_POINTS
                ),
            ));
        }

        if self.transition_wave_multiplier < constants::MIN_TRANSITION_WAVE_MULTIPLIER
            || self.transition_wave_multiplier > constants::MAX_TRANSITION_WAVE_MULTIPLIER
        {
            return Err(ConfigurationError::invalid_generation_config(
                "transition_wave_multiplier",
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_TRANSITION_WAVE_MULTIPLIER,
                    constants::MAX_TRANSITION_WAVE_MULTIPLIER
                ),
            ));
        }

        Ok(())
    }

    /// Create a high-quality configuration for detailed generation
    #[must_use]
    pub const fn high_quality() -> Self {
        Self {
            serpentine_points: 200,
            optimization_points: 100,
            smooth_straight_middle_points: 20,
            transition_wave_multiplier: 2.0,
        }
    }

    /// Create a fast configuration for quick generation
    #[must_use]
    pub const fn fast() -> Self {
        Self {
            serpentine_points: 50,
            optimization_points: 25,
            smooth_straight_middle_points: 5,
            transition_wave_multiplier: 2.0,
        }
    }
}

/// Configuration for basic geometry parameters
///
/// This struct defines the fundamental geometric constraints for microfluidic
/// channel generation, including wall clearances and channel dimensions.
///
/// # Examples
///
/// ```rust
/// use scheme::config::GeometryConfig;
///
/// // Create with default values
/// let config = GeometryConfig::default();
///
/// // Create with custom values
/// let custom_config = GeometryConfig::new(0.5, 1.0, 0.5).unwrap();
/// ```
#[derive(Clone, Copy, Debug)]
pub struct GeometryConfig {
    /// Minimum distance between channels and walls (mm)
    pub wall_clearance: f64,
    /// Width of channels (mm)
    pub channel_width: f64,
    /// Height of channels (mm) - used for 3D-aware calculations
    pub channel_height: f64,
    /// Configuration for geometry generation quality and precision
    pub generation: GeometryGenerationConfig,
}

impl GeometryConfig {
    /// Create a new geometry configuration with validation
    pub fn new(
        wall_clearance: f64,
        channel_width: f64,
        channel_height: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            wall_clearance,
            channel_width,
            channel_height,
            generation: GeometryGenerationConfig::default(),
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a new geometry configuration with custom generation settings
    pub fn with_generation(
        wall_clearance: f64,
        channel_width: f64,
        channel_height: f64,
        generation: GeometryGenerationConfig,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            wall_clearance,
            channel_width,
            channel_height,
            generation,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the geometry configuration
    pub fn validate(&self) -> ConfigurationResult<()> {
        if self.wall_clearance < constants::MIN_WALL_CLEARANCE
            || self.wall_clearance > constants::MAX_WALL_CLEARANCE
        {
            return Err(ConfigurationError::invalid_geometry_config(
                "wall_clearance",
                self.wall_clearance,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_WALL_CLEARANCE,
                    constants::MAX_WALL_CLEARANCE
                ),
            ));
        }

        if self.channel_width < constants::MIN_CHANNEL_WIDTH
            || self.channel_width > constants::MAX_CHANNEL_WIDTH
        {
            return Err(ConfigurationError::invalid_geometry_config(
                "channel_width",
                self.channel_width,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_CHANNEL_WIDTH,
                    constants::MAX_CHANNEL_WIDTH
                ),
            ));
        }

        if self.channel_height < constants::MIN_CHANNEL_HEIGHT
            || self.channel_height > constants::MAX_CHANNEL_HEIGHT
        {
            return Err(ConfigurationError::invalid_geometry_config(
                "channel_height",
                self.channel_height,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_CHANNEL_HEIGHT,
                    constants::MAX_CHANNEL_HEIGHT
                ),
            ));
        }

        // Validate generation configuration
        self.generation.validate()?;

        Ok(())
    }
}

impl Default for GeometryConfig {
    fn default() -> Self {
        Self {
            wall_clearance: constants::DEFAULT_WALL_CLEARANCE,
            channel_width: constants::DEFAULT_CHANNEL_WIDTH,
            channel_height: constants::DEFAULT_CHANNEL_HEIGHT,
            generation: GeometryGenerationConfig::default(),
        }
    }
}

/// Optimization profile for serpentine channel optimization
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd)]
pub enum OptimizationProfile {
    /// Fast optimization with limited parameter exploration (5-10x slower)
    Fast,
    /// Balanced optimization with moderate exploration (20-50x slower)
    Balanced,
    /// Thorough optimization with extensive exploration (100-500x slower)
    Thorough,
}

/// Wave shape types for serpentine channels
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum WaveShape {
    /// Smooth sine wave (default) - provides natural, flowing curves
    #[default]
    Sine,
    /// Square wave - provides sharp, angular transitions
    Square,
}

/// Taper profile types for frustum channels
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum TaperProfile {
    /// Linear taper - constant rate of width change
    #[default]
    Linear,
    /// Exponential taper - smooth exponential transition
    Exponential,
    /// Smooth taper - uses cosine function for very smooth transitions
    Smooth,
}

/// Configuration for serpentine (S-shaped) channels
#[derive(Debug, Clone, Copy)]
pub struct SerpentineConfig {
    /// Fraction of available vertical space to fill (0.1 to 0.95)
    pub fill_factor: f64,
    /// Multiplier for channel width to determine wavelength (1.0 to 10.0)
    pub wavelength_factor: f64,
    /// Controls width of Gaussian envelope - sigma = length / `gaussian_width_factor` (2.0 to 20.0)
    pub gaussian_width_factor: f64,
    /// Controls wave density relative to channel length - higher = more waves (0.5 to 10.0)
    pub wave_density_factor: f64,
    /// Controls wave phase direction for symmetry: -1.0=force inward, 1.0=force outward, 0.0=auto-symmetric
    pub wave_phase_direction: f64,
    /// Wave shape type - sine for smooth curves, square for angular transitions
    pub wave_shape: WaveShape,
    /// Enable length optimization algorithm (default: false for backward compatibility)
    pub optimization_enabled: bool,
    /// Target fill ratio for optimization - fraction of maximum possible length to achieve (0.8 to 0.99)
    pub target_fill_ratio: f64,
    /// Optimization profile controlling speed vs quality tradeoff
    pub optimization_profile: OptimizationProfile,
    /// Adaptive behavior configuration for dynamic channel properties
    pub adaptive_config: AdaptiveSerpentineConfig,
}

impl SerpentineConfig {
    /// Create a new serpentine configuration with validation
    pub fn new(
        fill_factor: f64,
        wavelength_factor: f64,
        gaussian_width_factor: f64,
        wave_density_factor: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            fill_factor,
            wavelength_factor,
            gaussian_width_factor,
            wave_density_factor,
            wave_phase_direction: 0.0, // Auto-determine for perfect symmetry
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false, // Disabled by default for backward compatibility
            target_fill_ratio: 0.9,    // Default target for optimization
            optimization_profile: OptimizationProfile::Balanced, // Default profile
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a new serpentine configuration with optimization enabled
    pub fn new_with_optimization(
        fill_factor: f64,
        wavelength_factor: f64,
        gaussian_width_factor: f64,
        wave_density_factor: f64,
        target_fill_ratio: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            fill_factor,
            wavelength_factor,
            gaussian_width_factor,
            wave_density_factor,
            wave_phase_direction: 0.0, // Auto-determine for perfect symmetry
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: true,
            target_fill_ratio,
            optimization_profile: OptimizationProfile::Balanced, // Default profile
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a new serpentine configuration with optimization and profile
    pub fn new_with_optimization_profile(
        fill_factor: f64,
        wavelength_factor: f64,
        gaussian_width_factor: f64,
        wave_density_factor: f64,
        target_fill_ratio: f64,
        optimization_profile: OptimizationProfile,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            fill_factor,
            wavelength_factor,
            gaussian_width_factor,
            wave_density_factor,
            wave_phase_direction: 0.0, // Auto-determine for perfect symmetry
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: true,
            target_fill_ratio,
            optimization_profile,
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a new serpentine configuration with explicit wave phase direction control
    pub fn new_with_phase_direction(
        fill_factor: f64,
        wavelength_factor: f64,
        gaussian_width_factor: f64,
        wave_density_factor: f64,
        wave_phase_direction: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            fill_factor,
            wavelength_factor,
            gaussian_width_factor,
            wave_density_factor,
            wave_phase_direction,
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false,      // Disabled by default for backward compatibility
            target_fill_ratio: 0.9,           // Default target for optimization
            optimization_profile: OptimizationProfile::Balanced, // Default profile
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the serpentine configuration
    pub fn validate(&self) -> ConfigurationResult<()> {
        if self.fill_factor < constants::MIN_FILL_FACTOR
            || self.fill_factor > constants::MAX_FILL_FACTOR
        {
            return Err(ConfigurationError::invalid_serpentine_config(
                "fill_factor",
                self.fill_factor,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_FILL_FACTOR,
                    constants::MAX_FILL_FACTOR
                ),
            ));
        }

        if self.wavelength_factor < constants::MIN_WAVELENGTH_FACTOR
            || self.wavelength_factor > constants::MAX_WAVELENGTH_FACTOR
        {
            return Err(ConfigurationError::invalid_serpentine_config(
                "wavelength_factor",
                self.wavelength_factor,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_WAVELENGTH_FACTOR,
                    constants::MAX_WAVELENGTH_FACTOR
                ),
            ));
        }

        if self.gaussian_width_factor < constants::MIN_GAUSSIAN_WIDTH_FACTOR
            || self.gaussian_width_factor > constants::MAX_GAUSSIAN_WIDTH_FACTOR
        {
            return Err(ConfigurationError::invalid_serpentine_config(
                "gaussian_width_factor",
                self.gaussian_width_factor,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_GAUSSIAN_WIDTH_FACTOR,
                    constants::MAX_GAUSSIAN_WIDTH_FACTOR
                ),
            ));
        }

        if self.wave_density_factor < constants::MIN_WAVE_DENSITY_FACTOR
            || self.wave_density_factor > constants::MAX_WAVE_DENSITY_FACTOR
        {
            return Err(ConfigurationError::invalid_serpentine_config(
                "wave_density_factor",
                self.wave_density_factor,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_WAVE_DENSITY_FACTOR,
                    constants::MAX_WAVE_DENSITY_FACTOR
                ),
            ));
        }

        if self.wave_phase_direction.abs() > 1.0 {
            return Err(ConfigurationError::invalid_serpentine_config(
                "wave_phase_direction",
                self.wave_phase_direction,
                "Must be between -1.0 and 1.0",
            ));
        }

        if self.target_fill_ratio < 0.8 || self.target_fill_ratio > 0.99 {
            return Err(ConfigurationError::invalid_serpentine_config(
                "target_fill_ratio",
                self.target_fill_ratio,
                "Must be between 0.8 and 0.99",
            ));
        }

        Ok(())
    }

    /// Convert this configuration to use square wave shape
    #[must_use]
    pub const fn with_square_wave(mut self) -> Self {
        self.wave_shape = WaveShape::Square;
        self
    }

    /// Convert this configuration to use sine wave shape
    #[must_use]
    pub const fn with_sine_wave(mut self) -> Self {
        self.wave_shape = WaveShape::Sine;
        self
    }

    /// Convert this configuration to use the specified wave shape
    #[must_use]
    pub const fn with_wave_shape(mut self, wave_shape: WaveShape) -> Self {
        self.wave_shape = wave_shape;
        self
    }
}

/// Configuration for arc (curved) channels
#[derive(Debug, Clone, Copy)]
pub struct ArcConfig {
    /// Controls how curved the arc is - 0.0 = straight, 1.0 = semicircle (0.0 to 2.0)
    pub curvature_factor: f64,
    /// Number of points to generate along the arc - higher = smoother (3 to 1000)
    pub smoothness: usize,
    /// Controls the direction of arc curvature: 1.0 = outward, -1.0 = inward, 0.0 = auto
    pub curvature_direction: f64,
    /// Minimum separation distance between arc channels (0.1 to 10.0)
    pub min_separation_distance: f64,
    /// Enable collision detection and prevention (default: true)
    pub enable_collision_prevention: bool,
    /// Maximum allowed curvature reduction factor for collision prevention (0.1 to 1.0)
    pub max_curvature_reduction: f64,
    /// Enable adaptive curvature based on neighbor proximity (default: true)
    pub enable_adaptive_curvature: bool,
}

impl ArcConfig {
    /// Create a new arc configuration with validation
    pub fn new(curvature_factor: f64, smoothness: usize) -> ConfigurationResult<Self> {
        let config = Self {
            curvature_factor,
            smoothness,
            curvature_direction: 0.0, // Auto-determine direction
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a new arc configuration with explicit curvature direction
    pub fn new_with_direction(
        curvature_factor: f64,
        smoothness: usize,
        curvature_direction: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            curvature_factor,
            smoothness,
            curvature_direction,
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a new arc configuration with full proximity control
    pub fn new_with_proximity_control(
        curvature_factor: f64,
        smoothness: usize,
        curvature_direction: f64,
        min_separation_distance: f64,
        enable_collision_prevention: bool,
        max_curvature_reduction: f64,
        enable_adaptive_curvature: bool,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            curvature_factor,
            smoothness,
            curvature_direction,
            min_separation_distance,
            enable_collision_prevention,
            max_curvature_reduction,
            enable_adaptive_curvature,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the arc configuration
    pub fn validate(&self) -> ConfigurationResult<()> {
        if self.curvature_factor < constants::MIN_CURVATURE_FACTOR
            || self.curvature_factor > constants::MAX_CURVATURE_FACTOR
        {
            return Err(ConfigurationError::invalid_arc_config(
                "curvature_factor",
                self.curvature_factor,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_CURVATURE_FACTOR,
                    constants::MAX_CURVATURE_FACTOR
                ),
            ));
        }

        if self.smoothness < constants::MIN_SMOOTHNESS
            || self.smoothness > constants::MAX_SMOOTHNESS
        {
            return Err(ConfigurationError::InvalidArcConfig {
                field: "smoothness".to_string(),
                value: self.smoothness as f64,
                constraint: format!(
                    "Must be between {} and {}",
                    constants::MIN_SMOOTHNESS,
                    constants::MAX_SMOOTHNESS
                ),
            });
        }

        if self.curvature_direction.abs() > 1.0 {
            return Err(ConfigurationError::InvalidArcConfig {
                field: "curvature_direction".to_string(),
                value: self.curvature_direction,
                constraint: "Must be between -1.0 and 1.0".to_string(),
            });
        }

        if self.min_separation_distance < constants::MIN_SEPARATION_DISTANCE
            || self.min_separation_distance > constants::MAX_SEPARATION_DISTANCE
        {
            return Err(ConfigurationError::invalid_arc_config(
                "min_separation_distance",
                self.min_separation_distance,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_SEPARATION_DISTANCE,
                    constants::MAX_SEPARATION_DISTANCE
                ),
            ));
        }

        if self.max_curvature_reduction < constants::MIN_CURVATURE_REDUCTION
            || self.max_curvature_reduction > constants::MAX_CURVATURE_REDUCTION_LIMIT
        {
            return Err(ConfigurationError::invalid_arc_config(
                "max_curvature_reduction",
                self.max_curvature_reduction,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_CURVATURE_REDUCTION,
                    constants::MAX_CURVATURE_REDUCTION_LIMIT
                ),
            ));
        }

        Ok(())
    }
}

impl Default for SerpentineConfig {
    fn default() -> Self {
        Self {
            fill_factor: constants::DEFAULT_FILL_FACTOR,
            wavelength_factor: constants::DEFAULT_WAVELENGTH_FACTOR,
            gaussian_width_factor: constants::DEFAULT_GAUSSIAN_WIDTH_FACTOR,
            wave_density_factor: constants::DEFAULT_WAVE_DENSITY_FACTOR,
            wave_phase_direction: 0.0, // Auto-determine for perfect symmetry
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false, // Disabled by default for backward compatibility
            target_fill_ratio: 0.9,    // Default target for optimization
            optimization_profile: OptimizationProfile::Balanced, // Default profile
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        }
    }
}

impl Default for ArcConfig {
    fn default() -> Self {
        Self {
            curvature_factor: constants::DEFAULT_CURVATURE_FACTOR,
            smoothness: constants::DEFAULT_SMOOTHNESS,
            curvature_direction: 0.0, // Auto-determine direction
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        }
    }
}

/// Configuration for frustum (tapered) channels with venturi throat functionality
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FrustumConfig {
    /// Inlet width (starting width) - must be positive (0.1 to 50.0)
    pub inlet_width: f64,
    /// Throat width (minimum width at center) - must be less than inlet/outlet widths (0.05 to 25.0)
    pub throat_width: f64,
    /// Outlet width (ending width) - must be positive (0.1 to 50.0)
    pub outlet_width: f64,
    /// Taper profile type - controls how width changes along the channel
    pub taper_profile: TaperProfile,
    /// Number of points to generate along the path - higher = smoother (10 to 1000)
    pub smoothness: usize,
    /// Throat position factor - 0.5 = center, 0.0 = at inlet, 1.0 = at outlet (0.1 to 0.9)
    pub throat_position: f64,
}

impl FrustumConfig {
    /// Create a new `FrustumConfig` with validation
    ///
    /// # Arguments
    /// * `inlet_width` - Starting width of the channel (must be positive)
    /// * `throat_width` - Minimum width at the throat (must be less than inlet/outlet)
    /// * `outlet_width` - Ending width of the channel (must be positive)
    /// * `taper_profile` - Type of taper profile to use
    /// * `smoothness` - Number of points for path generation
    /// * `throat_position` - Position of throat along channel (0.5 = center)
    ///
    /// # Returns
    /// * `ConfigurationResult<Self>` - The validated configuration or an error
    pub fn new(
        inlet_width: f64,
        throat_width: f64,
        outlet_width: f64,
        taper_profile: TaperProfile,
        smoothness: usize,
        throat_position: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            inlet_width,
            throat_width,
            outlet_width,
            taper_profile,
            smoothness,
            throat_position,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the frustum configuration
    pub fn validate(&self) -> ConfigurationResult<()> {
        // Validate inlet width
        if self.inlet_width <= 0.0 || self.inlet_width > 50.0 {
            return Err(ConfigurationError::invalid_frustum_config(
                "inlet_width",
                self.inlet_width,
                "Inlet width must be positive and <= 50.0",
            ));
        }

        // Validate outlet width
        if self.outlet_width <= 0.0 || self.outlet_width > 50.0 {
            return Err(ConfigurationError::invalid_frustum_config(
                "outlet_width",
                self.outlet_width,
                "Outlet width must be positive and <= 50.0",
            ));
        }

        // Validate throat width
        if self.throat_width <= 0.0 || self.throat_width > 25.0 {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_width",
                self.throat_width,
                "Throat width must be positive and <= 25.0",
            ));
        }

        // Throat must be narrower than both inlet and outlet
        if self.throat_width >= self.inlet_width {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_width",
                self.throat_width,
                "Throat width must be less than inlet width",
            ));
        }

        if self.throat_width >= self.outlet_width {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_width",
                self.throat_width,
                "Throat width must be less than outlet width",
            ));
        }

        // Validate smoothness
        if self.smoothness < 10 || self.smoothness > 1000 {
            return Err(ConfigurationError::invalid_frustum_config(
                "smoothness",
                self.smoothness as f64,
                "Smoothness must be between 10 and 1000",
            ));
        }

        // Validate throat position
        if self.throat_position < 0.1 || self.throat_position > 0.9 {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_position",
                self.throat_position,
                "Throat position must be between 0.1 and 0.9",
            ));
        }

        Ok(())
    }

    /// Calculate width at a given position along the channel (0.0 to 1.0)
    #[must_use]
    pub fn width_at_position(&self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);

        match self.taper_profile {
            TaperProfile::Linear => self.linear_taper(t),
            TaperProfile::Exponential => self.exponential_taper(t),
            TaperProfile::Smooth => self.smooth_taper(t),
        }
    }

    /// Linear taper profile
    fn linear_taper(&self, t: f64) -> f64 {
        if t <= self.throat_position {
            // Inlet to throat
            let local_t = t / self.throat_position;
            (self.throat_width - self.inlet_width).mul_add(local_t, self.inlet_width)
        } else {
            // Throat to outlet
            let local_t = (t - self.throat_position) / (1.0 - self.throat_position);
            (self.outlet_width - self.throat_width).mul_add(local_t, self.throat_width)
        }
    }

    /// Exponential taper profile
    fn exponential_taper(&self, t: f64) -> f64 {
        if t <= self.throat_position {
            // Inlet to throat - exponential decay
            let local_t = t / self.throat_position;
            // Use a normalized exponential that goes from 1 to 0
            let exp_factor = 2.0;
            let factor =
                ((-exp_factor * local_t).exp() - (-exp_factor).exp()) / (1.0 - (-exp_factor).exp());
            (self.throat_width - self.inlet_width).mul_add(1.0 - factor, self.inlet_width)
        } else {
            // Throat to outlet - exponential growth
            let local_t = (t - self.throat_position) / (1.0 - self.throat_position);
            // Use a normalized exponential that goes from 0 to 1
            let exp_factor = 2.0;
            let factor = (1.0 - (-exp_factor * local_t).exp()) / (1.0 - (-exp_factor).exp());
            (self.outlet_width - self.throat_width).mul_add(factor, self.throat_width)
        }
    }

    /// Smooth taper profile using cosine function
    fn smooth_taper(&self, t: f64) -> f64 {
        if t <= self.throat_position {
            // Inlet to throat - smooth cosine transition
            let local_t = t / self.throat_position;
            let factor = 0.5 * (1.0 + (std::f64::consts::PI * local_t).cos());
            (self.inlet_width - self.throat_width).mul_add(factor, self.throat_width)
        } else {
            // Throat to outlet - smooth cosine transition
            let local_t = (t - self.throat_position) / (1.0 - self.throat_position);
            let factor = 0.5 * (1.0 - (std::f64::consts::PI * local_t).cos());
            (self.outlet_width - self.throat_width).mul_add(factor, self.throat_width)
        }
    }
}

impl Default for FrustumConfig {
    fn default() -> Self {
        Self {
            inlet_width: 2.0,
            throat_width: 0.5,
            outlet_width: 2.0,
            taper_profile: TaperProfile::Linear,
            smoothness: 50,
            throat_position: 0.5, // Center
        }
    }
}

/// Configuration for selecting channel types in microfluidic schematics
///
/// This enum provides different strategies for determining what type of channel
/// (straight, serpentine, or arc) to use for each connection in the system.
///
/// # Examples
///
/// ```rust
/// use scheme::config::{ChannelTypeConfig, SerpentineConfig, ArcConfig, FrustumConfig};
///
/// // All channels will be straight lines
/// let straight_config = ChannelTypeConfig::AllStraight;
///
/// // All channels will be serpentine with default parameters
/// let serpentine_config = ChannelTypeConfig::AllSerpentine(SerpentineConfig::default());
///
/// // All channels will be frustum (tapered) with default parameters
/// let frustum_config = ChannelTypeConfig::AllFrustum(FrustumConfig::default());
///
/// // Adaptive selection based on channel characteristics
/// let adaptive_config = ChannelTypeConfig::Adaptive {
///     serpentine_config: SerpentineConfig::default(),
///     arc_config: ArcConfig::default(),
///     frustum_config: FrustumConfig::default(),
/// };
/// ```
#[derive(Debug, Clone, Copy)]
pub enum ChannelTypeConfig {
    /// All channels will be straight lines
    AllStraight,
    /// All channels will be smooth straight lines with transition zones
    AllSmoothStraight(SmoothTransitionConfig),
    /// All channels will be serpentine with the specified configuration
    AllSerpentine(SerpentineConfig),
    /// All channels will be arcs with the specified configuration
    AllArcs(ArcConfig),
    /// All channels will be frustum (tapered) with the specified configuration
    AllFrustum(FrustumConfig),
    /// Channels are selected based on their position in the layout
    MixedByPosition {
        /// Fraction of the box width that defines the middle zone for serpentine channels (0.0 to 1.0)
        middle_zone_fraction: f64,
        /// Configuration for serpentine channels in the middle zone
        serpentine_config: SerpentineConfig,
        /// Configuration for arc channels outside the middle zone
        arc_config: ArcConfig,
    },
    /// Adaptive channel type selection based on channel characteristics
    Adaptive {
        /// Configuration for serpentine channels when selected by adaptive algorithm
        serpentine_config: SerpentineConfig,
        /// Configuration for arc channels when selected by adaptive algorithm
        arc_config: ArcConfig,
        /// Configuration for frustum channels when selected by adaptive algorithm
        frustum_config: FrustumConfig,
    },
    /// Smooth serpentine channels with smooth straight junction connectors
    SmoothSerpentineWithTransitions {
        /// Configuration for serpentine channels in branches
        serpentine_config: SerpentineConfig,
        /// Configuration for smooth straight channels in junction connectors
        smooth_straight_config: SmoothTransitionConfig,
    },
    /// Custom function for determining channel type based on endpoints and box dimensions
    #[allow(clippy::type_complexity)]
    Custom(fn(from: (f64, f64), to: (f64, f64), box_dims: (f64, f64)) -> ChannelType),
}

impl Default for ChannelTypeConfig {
    fn default() -> Self {
        Self::Adaptive {
            serpentine_config: SerpentineConfig::default(),
            arc_config: ArcConfig::default(),
            frustum_config: FrustumConfig::default(),
        }
    }
}

/// Configuration presets for common use cases
pub mod presets {
    use super::{
        constants, AdaptiveSerpentineConfig, ArcConfig, GeometryConfig, GeometryGenerationConfig,
        OptimizationProfile, SerpentineConfig, WaveShape,
    };

    /// Preset for microfluidic devices with fine features
    #[must_use]
    pub const fn fine_features() -> GeometryConfig {
        GeometryConfig {
            wall_clearance: 0.2,
            channel_width: 0.5,
            channel_height: 0.3,
            generation: GeometryGenerationConfig::high_quality(),
        }
    }

    /// Preset for standard microfluidic devices
    #[must_use]
    pub fn standard() -> GeometryConfig {
        GeometryConfig::default()
    }

    /// Preset for large-scale microfluidic devices
    #[must_use]
    pub const fn large_scale() -> GeometryConfig {
        GeometryConfig {
            wall_clearance: 2.0,
            channel_width: 5.0,
            channel_height: 2.0,
            generation: GeometryGenerationConfig::fast(),
        }
    }

    /// Preset for high-density serpentine channels
    #[must_use]
    pub fn high_density_serpentine() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.9,
            wavelength_factor: 2.0,
            gaussian_width_factor: 8.0,
            wave_density_factor: 4.0,
            wave_phase_direction: 0.0,        // Auto-symmetric
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false,
            target_fill_ratio: 0.9,
            optimization_profile: OptimizationProfile::Balanced,
            adaptive_config: AdaptiveSerpentineConfig::aggressive(), // High-density needs aggressive adaptation
        }
    }

    /// Preset for smooth serpentine channels
    #[must_use]
    pub fn smooth_serpentine() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.6,
            wavelength_factor: 4.0,
            gaussian_width_factor: 10.0,
            wave_density_factor: 1.0, // Lower than default (1.5) for smoother, less dense waves
            wave_phase_direction: 0.0, // Auto-symmetric
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false,
            target_fill_ratio: 0.9,
            optimization_profile: OptimizationProfile::Balanced,
            adaptive_config: AdaptiveSerpentineConfig::conservative(), // Smooth channels need conservative adaptation
        }
    }

    /// Preset for inward-phase serpentine channels (all waves start inward)
    #[must_use]
    pub fn inward_serpentines() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.8,
            wavelength_factor: 3.0,
            gaussian_width_factor: 6.0,
            wave_density_factor: 2.0,
            wave_phase_direction: -1.0,       // Force inward phase
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false,
            target_fill_ratio: 0.9,
            optimization_profile: OptimizationProfile::Balanced,
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        }
    }

    /// Preset for outward-phase serpentine channels (all waves start outward)
    #[must_use]
    pub fn outward_serpentines() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.8,
            wavelength_factor: 3.0,
            gaussian_width_factor: 6.0,
            wave_density_factor: 2.0,
            wave_phase_direction: 1.0,        // Force outward phase
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: false,
            target_fill_ratio: 0.9,
            optimization_profile: OptimizationProfile::Balanced,
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        }
    }

    /// Preset for length-optimized serpentine channels
    #[must_use]
    pub fn optimized_serpentine() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.8,
            wavelength_factor: 3.0,
            gaussian_width_factor: 6.0,
            wave_density_factor: 2.0,
            wave_phase_direction: 0.0,        // Auto-symmetric
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: true,
            target_fill_ratio: 0.95, // Aggressive optimization target
            optimization_profile: OptimizationProfile::Balanced,
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        }
    }

    /// Preset for fast-optimized serpentine channels
    #[must_use]
    pub fn fast_optimized_serpentine() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.8,
            wavelength_factor: 3.0,
            gaussian_width_factor: 6.0,
            wave_density_factor: 2.0,
            wave_phase_direction: 0.0,        // Auto-symmetric
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: true,
            target_fill_ratio: 0.9, // Moderate optimization target
            optimization_profile: OptimizationProfile::Fast,
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        }
    }

    /// Preset for thorough-optimized serpentine channels
    #[must_use]
    pub fn thorough_optimized_serpentine() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.8,
            wavelength_factor: 3.0,
            gaussian_width_factor: 6.0,
            wave_density_factor: 2.0,
            wave_phase_direction: 0.0,        // Auto-symmetric
            wave_shape: WaveShape::default(), // Default to sine wave
            optimization_enabled: true,
            target_fill_ratio: 0.98, // Very aggressive optimization target
            optimization_profile: OptimizationProfile::Thorough,
            adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
        }
    }

    /// Preset for square wave serpentine channels (angular transitions)
    #[must_use]
    pub fn square_wave_serpentine() -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: 0.8,
            wavelength_factor: 3.0,
            gaussian_width_factor: 6.0,
            wave_density_factor: 2.0,
            wave_phase_direction: 0.0,     // Auto-symmetric
            wave_shape: WaveShape::Square, // Sharp square wave transitions
            optimization_enabled: false,
            target_fill_ratio: 0.9,
            optimization_profile: OptimizationProfile::Balanced,
            adaptive_config: AdaptiveSerpentineConfig::default(),
        }
    }

    /// Preset for subtle arc channels
    #[must_use]
    pub const fn subtle_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 0.2,
            smoothness: 30,
            curvature_direction: 0.0, // Auto-determine
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        }
    }

    /// Preset for high-quality geometry generation
    #[must_use]
    pub const fn high_quality_generation() -> GeometryConfig {
        GeometryConfig {
            wall_clearance: constants::DEFAULT_WALL_CLEARANCE,
            channel_width: constants::DEFAULT_CHANNEL_WIDTH,
            channel_height: constants::DEFAULT_CHANNEL_HEIGHT,
            generation: GeometryGenerationConfig::high_quality(),
        }
    }

    /// Preset for fast geometry generation
    #[must_use]
    pub const fn fast_generation() -> GeometryConfig {
        GeometryConfig {
            wall_clearance: constants::DEFAULT_WALL_CLEARANCE,
            channel_width: constants::DEFAULT_CHANNEL_WIDTH,
            channel_height: constants::DEFAULT_CHANNEL_HEIGHT,
            generation: GeometryGenerationConfig::fast(),
        }
    }

    /// Preset for microfluidic research applications
    #[must_use]
    pub const fn research_grade() -> GeometryConfig {
        GeometryConfig {
            wall_clearance: 0.3,
            channel_width: 0.8,
            channel_height: 0.5,
            generation: GeometryGenerationConfig::high_quality(),
        }
    }

    /// Preset for industrial manufacturing applications
    #[must_use]
    pub const fn manufacturing_grade() -> GeometryConfig {
        GeometryConfig {
            wall_clearance: 1.0,
            channel_width: 2.0,
            channel_height: 1.5,
            generation: GeometryGenerationConfig::fast(),
        }
    }

    /// Preset for pronounced arc channels
    #[must_use]
    pub const fn pronounced_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 0.8,
            smoothness: 50,
            curvature_direction: 0.0, // Auto-determine
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        }
    }

    /// Preset for inward-curving arc channels (all arcs curve toward center)
    #[must_use]
    pub const fn inward_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 0.5,
            smoothness: 30,
            curvature_direction: -1.0, // Force inward curvature
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        }
    }

    /// Preset for outward-curving arc channels (all arcs curve away from center)
    #[must_use]
    pub const fn outward_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 0.5,
            smoothness: 30,
            curvature_direction: 1.0, // Force outward curvature
            min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
            enable_collision_prevention: true,
            max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
            enable_adaptive_curvature: true,
        }
    }

    /// Preset for safe high-curvature arcs with collision prevention
    #[must_use]
    pub const fn safe_high_curvature_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 1.5, // High curvature but with safety measures
            smoothness: 50,
            curvature_direction: 0.0,     // Auto-determine
            min_separation_distance: 2.0, // Increased separation for safety
            enable_collision_prevention: true,
            max_curvature_reduction: 0.3, // Allow significant reduction if needed
            enable_adaptive_curvature: true,
        }
    }

    /// Preset for maximum curvature with aggressive collision prevention
    #[must_use]
    pub const fn maximum_safe_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 2.0, // Maximum curvature
            smoothness: 100,
            curvature_direction: 0.0,     // Auto-determine
            min_separation_distance: 3.0, // Maximum separation for safety
            enable_collision_prevention: true,
            max_curvature_reduction: 0.1, // Allow maximum reduction if needed
            enable_adaptive_curvature: true,
        }
    }

    /// Preset for dense layouts with conservative curvature
    #[must_use]
    pub const fn dense_layout_arcs() -> ArcConfig {
        ArcConfig {
            curvature_factor: 0.3, // Conservative curvature
            smoothness: 25,
            curvature_direction: 0.0,     // Auto-determine
            min_separation_distance: 0.5, // Tight separation for dense layouts
            enable_collision_prevention: true,
            max_curvature_reduction: 0.7, // Allow significant reduction
            enable_adaptive_curvature: true,
        }
    }
}

/// Builder for complex channel type configurations
pub struct ChannelTypeConfigBuilder {
    serpentine_config: SerpentineConfig,
    arc_config: ArcConfig,
    frustum_config: FrustumConfig,
    middle_zone_fraction: f64,
}

impl ChannelTypeConfigBuilder {
    /// Create a new builder with default configurations
    #[must_use]
    pub fn new() -> Self {
        Self {
            serpentine_config: SerpentineConfig::default(),
            arc_config: ArcConfig::default(),
            frustum_config: FrustumConfig::default(),
            middle_zone_fraction: constants::strategy_thresholds::DEFAULT_MIDDLE_ZONE_FRACTION,
        }
    }

    /// Set the serpentine configuration
    #[must_use]
    pub const fn with_serpentine_config(mut self, config: SerpentineConfig) -> Self {
        self.serpentine_config = config;
        self
    }

    /// Set the arc configuration
    #[must_use]
    pub const fn with_arc_config(mut self, config: ArcConfig) -> Self {
        self.arc_config = config;
        self
    }

    /// Set the frustum configuration
    #[must_use]
    pub const fn with_frustum_config(mut self, config: FrustumConfig) -> Self {
        self.frustum_config = config;
        self
    }

    /// Set the middle zone fraction for mixed by position
    #[must_use]
    pub const fn with_middle_zone_fraction(mut self, fraction: f64) -> Self {
        self.middle_zone_fraction = fraction;
        self
    }

    /// Build an adaptive channel type configuration
    #[must_use]
    pub const fn build_adaptive(self) -> ChannelTypeConfig {
        ChannelTypeConfig::Adaptive {
            serpentine_config: self.serpentine_config,
            arc_config: self.arc_config,
            frustum_config: self.frustum_config,
        }
    }

    /// Build a mixed by position channel type configuration
    #[must_use]
    pub const fn build_mixed_by_position(self) -> ChannelTypeConfig {
        ChannelTypeConfig::MixedByPosition {
            middle_zone_fraction: self.middle_zone_fraction,
            serpentine_config: self.serpentine_config,
            arc_config: self.arc_config,
        }
    }
}

impl Default for ChannelTypeConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}
