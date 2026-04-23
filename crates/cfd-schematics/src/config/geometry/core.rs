use super::generation::GeometryGenerationConfig;
use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};

/// Main configuration for converting schematics to geometry
///
/// This struct aggregates all settings required to generate the physical layout
/// of the microfluidic device, including dimensions, clearances, and quality settings.
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::config::{GeometryConfig, GeometryGenerationConfig};
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
    /// Minimum clearance between channels and walls/other channels
    /// (`constants::MIN_WALL_CLEARANCE..=constants::MAX_WALL_CLEARANCE` mm).
    pub wall_clearance: f64,
    /// Width of the fluidic channels
    /// (`constants::MIN_CHANNEL_WIDTH..=constants::MAX_CHANNEL_WIDTH` mm).
    pub channel_width: f64,
    /// Height (depth) of the fluidic channels
    /// (`constants::MIN_CHANNEL_HEIGHT..=constants::MAX_CHANNEL_HEIGHT` mm).
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

    /// Validate the configuration parameters.
    ///
    /// # Theorem - Canonical Geometry Contract
    ///
    /// A geometry configuration is admissible iff each scalar lies within its
    /// canonical range and `wall_clearance < channel_width`.
    ///
    /// **Proof sketch**: The range checks keep every scalar within the bounded
    /// physical domain used by the generator, while the strict ordering keeps
    /// the planar layout from degenerating to zero available spacing.
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

        if self.wall_clearance >= self.channel_width {
            return Err(ConfigurationError::invalid_geometry_config(
                "wall_clearance",
                self.wall_clearance,
                "Must be less than channel_width",
            ));
        }

        // Validate nested configuration
        self.generation.validate()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::constants::primitives as constants;

    #[test]
    fn geometry_config_accepts_canonical_bounds() {
        let config = GeometryConfig::new(
            constants::MIN_WALL_CLEARANCE,
            constants::MIN_CHANNEL_WIDTH * 2.0,
            constants::MAX_CHANNEL_HEIGHT,
            GeometryGenerationConfig::default(),
        )
        .expect("canonical bounds should validate");

        assert_eq!(config.wall_clearance, constants::MIN_WALL_CLEARANCE);
        assert_eq!(config.channel_width, constants::MIN_CHANNEL_WIDTH * 2.0);
        assert_eq!(config.channel_height, constants::MAX_CHANNEL_HEIGHT);
    }

    #[test]
    fn geometry_config_rejects_degenerate_clearance_width_relation() {
        let err = GeometryConfig::new(
            constants::MIN_WALL_CLEARANCE,
            constants::MIN_WALL_CLEARANCE,
            constants::DEFAULT_CHANNEL_HEIGHT,
            GeometryGenerationConfig::default(),
        )
        .expect_err("wall_clearance must be strictly less than channel_width");

        match err {
            ConfigurationError::InvalidGeometryConfig {
                field,
                value,
                constraint,
            } => {
                assert_eq!(field, "wall_clearance");
                assert_eq!(value, constants::MIN_WALL_CLEARANCE);
                assert!(constraint.contains("channel_width"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
