//! Wave generation constants for serpentine channel shaping

use super::primitives;

/// Wave generation constants previously hardcoded in strategies
///
/// The values themselves live in [`primitives`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaveGenerationConstants {
    /// Sharpness factor for square wave generation
    pub square_wave_sharpness: f64,

    /// Transition zone factor for smooth endpoints
    pub transition_zone_factor: f64,

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

    /// Neighbor avoidance scaling factor
    pub neighbor_avoidance_scaling_factor: f64,

    /// Geometric tolerance for distance comparisons
    pub geometric_tolerance: f64,
}

impl WaveGenerationConstants {
    /// Canonical default wave generation constants.
    pub const DEFAULT: Self = Self {
        square_wave_sharpness: primitives::SQUARE_WAVE_SHARPNESS,
        transition_zone_factor: primitives::TRANSITION_ZONE_FACTOR,
        smooth_endpoint_start_threshold: primitives::SMOOTH_ENDPOINT_START_THRESHOLD,
        smooth_endpoint_end_threshold: primitives::SMOOTH_ENDPOINT_END_THRESHOLD,
        default_transition_length_factor: primitives::DEFAULT_TRANSITION_LENGTH_FACTOR,
        default_transition_amplitude_factor: primitives::DEFAULT_TRANSITION_AMPLITUDE_FACTOR,
        default_transition_smoothness: primitives::DEFAULT_TRANSITION_SMOOTHNESS,
        default_wave_multiplier: primitives::DEFAULT_TRANSITION_WAVE_MULTIPLIER,
        neighbor_avoidance_scaling_factor: primitives::NEIGHBOR_AVOIDANCE_SCALING_FACTOR,
        geometric_tolerance: primitives::GEOMETRIC_TOLERANCE,
    };
}

impl Default for WaveGenerationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
