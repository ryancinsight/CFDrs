use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};
use super::generation::GeometryGenerationConfig;

/// Main configuration for converting schematics to geometry
///
/// This struct aggregates all settings required to generate the physical layout
/// of the microfluidic device, including dimensions, clearances, and quality settings.
///
/// # Examples
///
/// ```rust
/// use scheme::config::{GeometryConfig, GeometryGenerationConfig};
///
/// // Create with default settings
/// let config = GeometryConfig::default();
///
/// // Create with custom dimensions
/// let custom_config = GeometryConfig {
///     wall_clearance: 1.0,
///     channel_width: 2.0,
///     channel_height: 1.5,
///     generation: GeometryGenerationConfig::high_quality(),
/// };
/// ```
#[derive(Clone, Copy, Debug)]
pub struct GeometryConfig {
    /// Minimum clearance between channels and walls/other channels (0.1-100.0 mm)
    pub wall_clearance: f64,
    /// Width of the fluidic channels (0.1-50.0 mm)
    pub channel_width: f64,
    /// Height (depth) of the fluidic channels (0.1-50.0 mm)
    pub channel_height: f64,
    /// Settings for geometry generation quality and precision
    pub generation: GeometryGenerationConfig,
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

impl GeometryConfig {
    /// Create a new geometry configuration with validation
    pub fn new(
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

    /// Validate the configuration parameters
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

        if self.channel_width < 0.1 || self.channel_width > 50.0 {
            return Err(ConfigurationError::invalid_geometry_config(
                "channel_width",
                self.channel_width,
                "Must be between 0.1 and 50.0",
            ));
        }

        if self.channel_height < 0.1 || self.channel_height > 50.0 {
            return Err(ConfigurationError::invalid_geometry_config(
                "channel_height",
                self.channel_height,
                "Must be between 0.1 and 50.0",
            ));
        }

        // Validate nested configuration
        self.generation.validate()?;

        Ok(())
    }
}
