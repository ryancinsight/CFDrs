//! Straight and smooth-straight channel strategies.

use crate::config::{ConstantsRegistry, GeometryConfig};
use crate::geometry::{ChannelType, Point2D};

use super::ChannelTypeStrategy;

/// Strategy for creating straight channels
#[derive(Debug, Clone)]
pub struct StraightChannelStrategy;

impl ChannelTypeStrategy for StraightChannelStrategy {
    fn create_channel(
        &self,
        _from: Point2D,
        _to: Point2D,
        _geometry_config: &GeometryConfig,
        _box_dims: (f64, f64),
        _total_branches: usize,
        _neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        ChannelType::Straight
    }
}

/// Strategy for creating smooth straight channels with transition zones
#[derive(Debug, Clone)]
pub struct SmoothStraightChannelStrategy {
    /// Configuration for transition zones
    pub transition_config: SmoothTransitionConfig,
}

/// Configuration for smooth transition zones in straight channels
///
/// This configuration controls the smooth transition zones at the endpoints
/// of straight channels to eliminate sharp corners when connecting to other
/// channel types.
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::geometry::strategies::SmoothTransitionConfig;
///
/// // Create with default values
/// let config = SmoothTransitionConfig::default();
///
/// // Create with custom values for subtle transitions
/// let subtle = SmoothTransitionConfig {
///     transition_length_factor: 0.1,
///     transition_amplitude_factor: 0.2,
///     transition_smoothness: 15,
///     wave_multiplier: 1.5,
/// };
/// ```
#[derive(Debug, Clone, Copy)]
pub struct SmoothTransitionConfig {
    /// Length of transition zone as fraction of channel length (0.0 to 0.5)
    pub transition_length_factor: f64,
    /// Maximum amplitude of transition waves relative to channel width (0.0 to 1.0)
    pub transition_amplitude_factor: f64,
    /// Number of points to use for transition smoothing (5 to 100)
    pub transition_smoothness: usize,
    /// Wave multiplier for transition waves (0.5 to 10.0, where 2.0 = one complete wave)
    pub wave_multiplier: f64,
}

impl Default for SmoothTransitionConfig {
    fn default() -> Self {
        let constants = ConstantsRegistry::new();
        Self {
            transition_length_factor: constants.get_default_transition_length_factor(),
            transition_amplitude_factor: constants.get_default_transition_amplitude_factor(),
            transition_smoothness: constants.get_default_transition_smoothness(),
            wave_multiplier: constants.get_default_wave_multiplier(),
        }
    }
}

impl SmoothTransitionConfig {
    /// Create a new smooth transition configuration with validation
    ///
    /// # Errors
    ///
    /// Returns an error if any parameter is outside its valid range.
    pub fn new(
        transition_length_factor: f64,
        transition_amplitude_factor: f64,
        transition_smoothness: usize,
        wave_multiplier: f64,
    ) -> Result<Self, String> {
        let config = Self {
            transition_length_factor,
            transition_amplitude_factor,
            transition_smoothness,
            wave_multiplier,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the configuration parameters
    ///
    /// # Errors
    ///
    /// Returns an error if any parameter is outside its valid range.
    pub fn validate(&self) -> Result<(), String> {
        if self.transition_length_factor < 0.0 || self.transition_length_factor > 0.5 {
            return Err("transition_length_factor must be between 0.0 and 0.5".to_string());
        }

        if self.transition_amplitude_factor < 0.0 || self.transition_amplitude_factor > 1.0 {
            return Err("transition_amplitude_factor must be between 0.0 and 1.0".to_string());
        }

        if self.transition_smoothness < 5 || self.transition_smoothness > 100 {
            return Err("transition_smoothness must be between 5 and 100".to_string());
        }

        if self.wave_multiplier < 0.5 || self.wave_multiplier > 10.0 {
            return Err("wave_multiplier must be between 0.5 and 10.0".to_string());
        }

        Ok(())
    }

    /// Create a subtle transition configuration
    #[must_use]
    pub const fn subtle() -> Self {
        Self {
            transition_length_factor: 0.1,
            transition_amplitude_factor: 0.2,
            transition_smoothness: 15,
            wave_multiplier: 1.5,
        }
    }

    /// Create a pronounced transition configuration
    #[must_use]
    pub const fn pronounced() -> Self {
        Self {
            transition_length_factor: 0.25,
            transition_amplitude_factor: 0.5,
            transition_smoothness: 30,
            wave_multiplier: 3.0,
        }
    }

    /// Create a high-quality transition configuration for detailed work
    #[must_use]
    pub const fn high_quality() -> Self {
        Self {
            transition_length_factor: 0.2,
            transition_amplitude_factor: 0.4,
            transition_smoothness: 50,
            wave_multiplier: 2.0,
        }
    }

    /// Create a fast transition configuration for quick generation
    #[must_use]
    pub const fn fast() -> Self {
        Self {
            transition_length_factor: 0.15,
            transition_amplitude_factor: 0.3,
            transition_smoothness: 10,
            wave_multiplier: 2.0,
        }
    }
}

impl SmoothStraightChannelStrategy {
    #[must_use]
    pub const fn new(config: SmoothTransitionConfig) -> Self {
        Self {
            transition_config: config,
        }
    }
}

impl ChannelTypeStrategy for SmoothStraightChannelStrategy {
    fn create_channel(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        _box_dims: (f64, f64),
        _total_branches: usize,
        _neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        let path = self.generate_smooth_straight_path(from, to, geometry_config);
        ChannelType::SmoothStraight { path }
    }
}

impl SmoothStraightChannelStrategy {
    /// Generate a smooth straight path with transition zones at endpoints
    fn generate_smooth_straight_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
    ) -> Vec<Point2D> {
        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let channel_length = dx.hypot(dy);

        let constants = ConstantsRegistry::new();
        // For very short channels, just return straight line
        if channel_length
            < geometry_config.channel_width * constants.get_short_channel_width_multiplier()
        {
            return vec![p1, p2];
        }

        let transition_length = channel_length * self.transition_config.transition_length_factor;
        let max_amplitude =
            geometry_config.channel_width * self.transition_config.transition_amplitude_factor;

        // Calculate total points: transition + middle + transition
        let transition_points = self.transition_config.transition_smoothness;
        let middle_points = geometry_config.generation.smooth_straight_middle_points;
        let total_points = transition_points * 2 + middle_points;

        let mut path = Vec::with_capacity(total_points);

        // Perpendicular direction for wave displacement
        let perp_x = -dy / channel_length;
        let perp_y = dx / channel_length;

        for i in 0..total_points {
            let t = i as f64 / (total_points - 1) as f64;

            // Base position along the line
            let base_x = t.mul_add(dx, p1.0);
            let base_y = t.mul_add(dy, p1.1);

            // Calculate smooth transition amplitude
            let amplitude = self.calculate_transition_amplitude(
                t,
                transition_length / channel_length,
                max_amplitude,
            );

            // Apply small wave for smooth transition
            let wave_phase = std::f64::consts::PI * self.transition_config.wave_multiplier * t;
            let wave_amplitude = amplitude * wave_phase.sin();

            let x = wave_amplitude.mul_add(perp_x, base_x);
            let y = wave_amplitude.mul_add(perp_y, base_y);

            // Ensure exact endpoints
            if i == 0 {
                path.push(p1);
            } else if i == total_points - 1 {
                path.push(p2);
            } else {
                path.push((x, y));
            }
        }

        path
    }

    /// Calculate transition amplitude that smoothly goes to zero at endpoints
    fn calculate_transition_amplitude(
        &self,
        t: f64,
        transition_factor: f64,
        max_amplitude: f64,
    ) -> f64 {
        // Create smooth transitions at both ends
        let start_transition = if t < transition_factor {
            // Smooth ramp up from 0 to 1 using smoothstep
            let local_t = t / transition_factor;
            local_t * local_t * 2.0f64.mul_add(-local_t, 3.0)
        } else {
            1.0
        };

        let end_transition = if t > (1.0 - transition_factor) {
            // Smooth ramp down from 1 to 0 using smoothstep
            let local_t = (1.0 - t) / transition_factor;
            local_t * local_t * 2.0f64.mul_add(-local_t, 3.0)
        } else {
            1.0
        };

        max_amplitude * start_transition * end_transition
    }
}
