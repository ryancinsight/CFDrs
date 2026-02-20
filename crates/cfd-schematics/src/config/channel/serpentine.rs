use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};
use crate::config::adaptive::AdaptiveSerpentineConfig;

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
