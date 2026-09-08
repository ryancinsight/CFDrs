use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};

/// Configuration for geometry generation parameters
///
/// This struct controls the quality and precision of geometry generation,
/// including point density for path generation and wave calculations.
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::config::GeometryGenerationConfig;
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
