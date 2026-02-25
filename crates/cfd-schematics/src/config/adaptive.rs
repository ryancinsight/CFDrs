use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};

/// Configuration for adaptive serpentine channel behavior
///
/// This configuration controls how serpentine channels adapt their wave properties
/// based on geometric constraints such as distance from nodes, walls, and other channels.
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::config::AdaptiveSerpentineConfig;
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
