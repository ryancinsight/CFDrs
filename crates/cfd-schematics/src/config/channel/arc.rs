use crate::config::constants::primitives as constants;
use crate::error::{ConfigurationError, ConfigurationResult};

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
