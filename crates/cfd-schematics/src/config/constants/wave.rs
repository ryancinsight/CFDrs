//! Wave generation constants for serpentine channel shaping

/// Wave generation constants previously hardcoded in strategies
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaveGenerationConstants {
    /// Sharpness factor for square wave generation
    pub square_wave_sharpness: f64,

    /// Transition zone factor for smooth endpoints
    pub transition_zone_factor: f64,

    /// Gaussian envelope scaling factor
    pub gaussian_envelope_scale: f64,

    /// Phase direction calculation threshold
    pub phase_direction_threshold: f64,

    /// Wave amplitude safety margin
    pub amplitude_safety_margin: f64,

    /// Smooth endpoint transition start threshold
    pub smooth_endpoint_start_threshold: f64,

    /// Smooth endpoint transition end threshold
    pub smooth_endpoint_end_threshold: f64,

    /// Default transition length factor for smooth transitions
    pub default_transition_length_factor: f64,

    /// Default transition amplitude factor
    pub default_transition_amplitude_factor: f64,

    /// Default transition smoothness points
    pub default_transition_smoothness: usize,

    /// Default wave multiplier for transitions
    pub default_wave_multiplier: f64,

    /// Wall proximity scaling factor
    pub wall_proximity_scaling_factor: f64,

    /// Neighbor avoidance scaling factor
    pub neighbor_avoidance_scaling_factor: f64,

    /// Geometric tolerance for distance comparisons
    pub geometric_tolerance: f64,
}

impl WaveGenerationConstants {
    /// Canonical default wave generation constants.
    pub const DEFAULT: Self = Self {
        square_wave_sharpness: 5.0,
        transition_zone_factor: 0.1,
        gaussian_envelope_scale: 1.0,
        phase_direction_threshold: 0.5,
        amplitude_safety_margin: 0.8,
        smooth_endpoint_start_threshold: 0.1,
        smooth_endpoint_end_threshold: 0.9,
        default_transition_length_factor: 0.15,
        default_transition_amplitude_factor: 0.3,
        default_transition_smoothness: 20,
        default_wave_multiplier: 2.0,
        wall_proximity_scaling_factor: 0.8,
        neighbor_avoidance_scaling_factor: 0.8,
        geometric_tolerance: 1e-6,
    };
}

impl Default for WaveGenerationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
