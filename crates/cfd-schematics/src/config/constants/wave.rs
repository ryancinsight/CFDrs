//! Wave generation constants for serpentine channel shaping

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

/// Wave generation constants previously hardcoded in strategies
pub struct WaveGenerationConstants {
    /// Sharpness factor for square wave generation
    pub square_wave_sharpness: ConfigurableParameter<f64>,

    /// Transition zone factor for smooth endpoints
    pub transition_zone_factor: ConfigurableParameter<f64>,

    /// Gaussian envelope scaling factor
    pub gaussian_envelope_scale: ConfigurableParameter<f64>,

    /// Phase direction calculation threshold
    pub phase_direction_threshold: ConfigurableParameter<f64>,

    /// Wave amplitude safety margin
    pub amplitude_safety_margin: ConfigurableParameter<f64>,

    /// Smooth endpoint transition start threshold
    pub smooth_endpoint_start_threshold: ConfigurableParameter<f64>,

    /// Smooth endpoint transition end threshold
    pub smooth_endpoint_end_threshold: ConfigurableParameter<f64>,

    /// Default transition length factor for smooth transitions
    pub default_transition_length_factor: ConfigurableParameter<f64>,

    /// Default transition amplitude factor
    pub default_transition_amplitude_factor: ConfigurableParameter<f64>,

    /// Default transition smoothness points
    pub default_transition_smoothness: ConfigurableParameter<usize>,

    /// Default wave multiplier for transitions
    pub default_wave_multiplier: ConfigurableParameter<f64>,

    /// Wall proximity scaling factor
    pub wall_proximity_scaling_factor: ConfigurableParameter<f64>,

    /// Neighbor avoidance scaling factor
    pub neighbor_avoidance_scaling_factor: ConfigurableParameter<f64>,

    /// Geometric tolerance for distance comparisons
    pub geometric_tolerance: ConfigurableParameter<f64>,
}

impl WaveGenerationConstants {
    /// Create default wave generation constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            square_wave_sharpness: ConfigurableParameter::new(
                5.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 20.0),
                ]),
                ParameterMetadata::new(
                    "square_wave_sharpness",
                    "Sharpness factor for square wave generation using tanh",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            transition_zone_factor: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 0.5),
                ]),
                ParameterMetadata::new(
                    "transition_zone_factor",
                    "Factor for smooth transition zones at wave endpoints",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            gaussian_envelope_scale: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 5.0),
                ]),
                ParameterMetadata::new(
                    "gaussian_envelope_scale",
                    "Scaling factor for Gaussian envelope calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            phase_direction_threshold: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "phase_direction_threshold",
                    "Threshold for phase direction calculation in bilateral symmetry",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            amplitude_safety_margin: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "amplitude_safety_margin",
                    "Safety margin factor for amplitude calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            smooth_endpoint_start_threshold: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 0.5),
                ]),
                ParameterMetadata::new(
                    "smooth_endpoint_start_threshold",
                    "Threshold for smooth endpoint transition start",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            smooth_endpoint_end_threshold: ConfigurableParameter::new(
                0.9,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.5, 0.99),
                ]),
                ParameterMetadata::new(
                    "smooth_endpoint_end_threshold",
                    "Threshold for smooth endpoint transition end",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            default_transition_length_factor: ConfigurableParameter::new(
                0.15,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.05, 0.5),
                ]),
                ParameterMetadata::new(
                    "default_transition_length_factor",
                    "Default length factor for smooth transitions",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            default_transition_amplitude_factor: ConfigurableParameter::new(
                0.3,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "default_transition_amplitude_factor",
                    "Default amplitude factor for smooth transitions",
                    "wave_generation",
                )
                .with_units("ratio"),
            ),

            default_transition_smoothness: ConfigurableParameter::new(
                20,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(5, 100),
                ]),
                ParameterMetadata::new(
                    "default_transition_smoothness",
                    "Default number of points for transition smoothing",
                    "wave_generation",
                )
                .with_units("points"),
            ),

            default_wave_multiplier: ConfigurableParameter::new(
                2.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.5, 10.0),
                ]),
                ParameterMetadata::new(
                    "default_wave_multiplier",
                    "Default wave multiplier for transitions",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            wall_proximity_scaling_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "wall_proximity_scaling_factor",
                    "Scaling factor for wall proximity calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            neighbor_avoidance_scaling_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "neighbor_avoidance_scaling_factor",
                    "Scaling factor for neighbor avoidance calculations",
                    "wave_generation",
                )
                .with_units("factor"),
            ),

            geometric_tolerance: ConfigurableParameter::new(
                1e-6,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1e-12, 1e-3),
                ]),
                ParameterMetadata::new(
                    "geometric_tolerance",
                    "Tolerance for geometric distance comparisons",
                    "wave_generation",
                )
                .with_units("units"),
            ),
        }
    }
}

impl Default for WaveGenerationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}
