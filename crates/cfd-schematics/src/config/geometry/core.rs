use super::generation::GeometryGenerationConfig;
use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};
use aequitas::systems::si::quantities::Length;

/// Main configuration for converting schematics to geometry
///
/// This struct aggregates all settings required to generate the physical layout
/// of the microfluidic device, including dimensions, clearances, and quality settings.
/// Physical dimensions are stored as Eunomia lengths in base metres. The
/// explicit `*_mm` methods are layout/formula projections for the legacy
/// millimetre coordinate space used by the generator.
///
/// # Examples
///
/// ```rust
/// use aequitas::systems::si::quantities::Length;
/// use cfd_schematics::config::{GeometryConfig, GeometryGenerationConfig};
///
/// // Create with default settings
/// let config = GeometryConfig::default();
///
/// // Create with custom dimensions
/// let custom_config = GeometryConfig {
///     wall_clearance: Length::from_base(1.0e-3),
///     channel_width: Length::from_base(2.0e-3),
///     channel_height: Length::from_base(1.5e-3),
///     generation: GeometryGenerationConfig::high_quality(),
/// };
/// ```
#[derive(Clone, Copy, Debug)]
pub struct GeometryConfig {
    /// Minimum clearance between channels and walls/other channels
    /// (`constants::MIN_WALL_CLEARANCE..=constants::MAX_WALL_CLEARANCE` mm).
    pub wall_clearance: Length<f64>,
    /// Width of the fluidic channels
    /// (`constants::MIN_CHANNEL_WIDTH..=constants::MAX_CHANNEL_WIDTH` mm).
    pub channel_width: Length<f64>,
    /// Height (depth) of the fluidic channels
    /// (`constants::MIN_CHANNEL_HEIGHT..=constants::MAX_CHANNEL_HEIGHT` mm).
    pub channel_height: Length<f64>,
    /// Settings for geometry generation quality and precision
    pub generation: GeometryGenerationConfig,
}

impl Default for GeometryConfig {
    fn default() -> Self {
        Self {
            wall_clearance: Length::from_base(constants::DEFAULT_WALL_CLEARANCE * 1.0e-3),
            channel_width: Length::from_base(constants::DEFAULT_CHANNEL_WIDTH * 1.0e-3),
            channel_height: Length::from_base(constants::DEFAULT_CHANNEL_HEIGHT * 1.0e-3),
            generation: GeometryGenerationConfig::default(),
        }
    }
}

impl GeometryConfig {
    /// Create a new geometry configuration with validation
    pub fn new(
        wall_clearance: Length<f64>,
        channel_width: Length<f64>,
        channel_height: Length<f64>,
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
    /// A geometry configuration is admissible iff each projected dimension lies
    /// within its canonical range and `wall_clearance < channel_width`.
    ///
    /// **Proof sketch**: The range checks keep every scalar within the bounded
    /// physical domain used by the generator, while the strict ordering keeps
    /// the planar layout from degenerating to zero available spacing.
    pub fn validate(&self) -> ConfigurationResult<()> {
        let wall_clearance_mm = self.wall_clearance_mm();
        let channel_width_mm = self.channel_width_mm();
        let channel_height_mm = self.channel_height_mm();

        if !(constants::MIN_WALL_CLEARANCE..=constants::MAX_WALL_CLEARANCE)
            .contains(&wall_clearance_mm)
        {
            return Err(ConfigurationError::invalid_geometry_config(
                "wall_clearance",
                wall_clearance_mm,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_WALL_CLEARANCE,
                    constants::MAX_WALL_CLEARANCE
                ),
            ));
        }

        if !(constants::MIN_CHANNEL_WIDTH..=constants::MAX_CHANNEL_WIDTH)
            .contains(&channel_width_mm)
        {
            return Err(ConfigurationError::invalid_geometry_config(
                "channel_width",
                channel_width_mm,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_CHANNEL_WIDTH,
                    constants::MAX_CHANNEL_WIDTH
                ),
            ));
        }

        if !(constants::MIN_CHANNEL_HEIGHT..=constants::MAX_CHANNEL_HEIGHT)
            .contains(&channel_height_mm)
        {
            return Err(ConfigurationError::invalid_geometry_config(
                "channel_height",
                channel_height_mm,
                &format!(
                    "Must be between {} and {}",
                    constants::MIN_CHANNEL_HEIGHT,
                    constants::MAX_CHANNEL_HEIGHT
                ),
            ));
        }

        if wall_clearance_mm >= channel_width_mm {
            return Err(ConfigurationError::invalid_geometry_config(
                "wall_clearance",
                wall_clearance_mm,
                "Must be less than channel_width",
            ));
        }

        // Validate nested configuration
        self.generation.validate()?;

        Ok(())
    }

    /// Project the authored wall clearance to millimetres for layout formulas.
    #[must_use]
    pub fn wall_clearance_mm(&self) -> f64 {
        self.wall_clearance.into_base() * 1.0e3
    }

    /// Project the authored channel width to millimetres for layout formulas.
    #[must_use]
    pub fn channel_width_mm(&self) -> f64 {
        self.channel_width.into_base() * 1.0e3
    }

    /// Project the authored channel height to millimetres for layout formulas.
    #[must_use]
    pub fn channel_height_mm(&self) -> f64 {
        self.channel_height.into_base() * 1.0e3
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::constants::primitives as constants;
    use aequitas::systems::si::quantities::Length;

    #[test]
    fn geometry_config_accepts_canonical_bounds() {
        let config = GeometryConfig::new(
            Length::from_base(constants::MIN_WALL_CLEARANCE * 1.0e-3),
            Length::from_base(constants::MIN_CHANNEL_WIDTH * 2.0e-3),
            Length::from_base(constants::MAX_CHANNEL_HEIGHT * 1.0e-3),
            GeometryGenerationConfig::default(),
        )
        .expect("canonical bounds should validate");

        assert_eq!(
            config.wall_clearance.into_base().to_bits(),
            (constants::MIN_WALL_CLEARANCE * 1.0e-3).to_bits()
        );
        assert_eq!(
            config.channel_width.into_base().to_bits(),
            (constants::MIN_CHANNEL_WIDTH * 2.0e-3).to_bits()
        );
        assert_eq!(
            config.channel_height.into_base().to_bits(),
            (constants::MAX_CHANNEL_HEIGHT * 1.0e-3).to_bits()
        );
    }

    #[test]
    fn geometry_config_rejects_degenerate_clearance_width_relation() {
        let err = GeometryConfig::new(
            Length::from_base(constants::MIN_WALL_CLEARANCE * 1.0e-3),
            Length::from_base(constants::MIN_WALL_CLEARANCE * 1.0e-3),
            Length::from_base(constants::DEFAULT_CHANNEL_HEIGHT * 1.0e-3),
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
                assert_eq!(value.to_bits(), constants::MIN_WALL_CLEARANCE.to_bits());
                assert!(constraint.contains("channel_width"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
