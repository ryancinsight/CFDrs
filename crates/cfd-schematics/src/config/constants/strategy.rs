//! Strategy selection thresholds and parameters

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

/// Strategy selection thresholds and parameters
pub struct StrategyThresholds {
    /// Minimum curvature factor to use arc strategy instead of straight
    pub arc_curvature_threshold: ConfigurableParameter<f64>,

    /// Maximum fill factor to use serpentine strategy
    pub serpentine_fill_threshold: ConfigurableParameter<f64>,

    /// Minimum channel length for complex strategies
    pub min_complex_strategy_length: ConfigurableParameter<f64>,

    /// Branch count threshold for adaptive behavior
    pub adaptive_branch_threshold: ConfigurableParameter<usize>,

    /// Minimum length threshold for frustum channel selection in smart mode
    pub frustum_min_length_threshold: ConfigurableParameter<f64>,

    /// Maximum length threshold for frustum channel selection in smart mode
    pub frustum_max_length_threshold: ConfigurableParameter<f64>,

    /// Maximum angle threshold for frustum channel selection (horizontal preference)
    pub frustum_angle_threshold: ConfigurableParameter<f64>,
}

impl StrategyThresholds {
    /// Create default strategy thresholds
    #[must_use]
    fn make_default() -> Self {
        Self {
            arc_curvature_threshold: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::non_negative(),
                    ParameterConstraints::range(0.0, 2.0),
                ]),
                ParameterMetadata::new(
                    "arc_curvature_threshold",
                    "Minimum curvature factor to trigger arc strategy selection",
                    "strategy_selection"
                ).with_units("factor").affects_others()
            ),
            
            serpentine_fill_threshold: ConfigurableParameter::new(
                0.95,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::normalized(),
                ]),
                ParameterMetadata::new(
                    "serpentine_fill_threshold",
                    "Maximum fill factor for serpentine strategy selection",
                    "strategy_selection"
                ).affects_others()
            ),
            
            min_complex_strategy_length: ConfigurableParameter::new(
                10.0,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 1000.0),
                ]),
                ParameterMetadata::new(
                    "min_complex_strategy_length",
                    "Minimum channel length to use complex strategies",
                    "strategy_selection"
                ).with_units("mm")
            ),
            
            adaptive_branch_threshold: ConfigurableParameter::new(
                4usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(1, 100),
                ]),
                ParameterMetadata::new(
                    "adaptive_branch_threshold",
                    "Branch count threshold to enable adaptive parameter behavior",
                    "strategy_selection"
                )
            ),

            frustum_min_length_threshold: ConfigurableParameter::new(
                0.3,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 1.0),
                ]),
                ParameterMetadata::new(
                    "frustum_min_length_threshold",
                    "Minimum length threshold (as fraction of box width) for frustum channel selection",
                    "strategy_selection"
                ).with_units("fraction")
            ),

            frustum_max_length_threshold: ConfigurableParameter::new(
                0.7,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.2, 1.0),
                ]),
                ParameterMetadata::new(
                    "frustum_max_length_threshold",
                    "Maximum length threshold (as fraction of box width) for frustum channel selection",
                    "strategy_selection"
                ).with_units("fraction")
            ),

            frustum_angle_threshold: ConfigurableParameter::new(
                0.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 2.0),
                ]),
                ParameterMetadata::new(
                    "frustum_angle_threshold",
                    "Maximum angle threshold for frustum channel selection (dy/dx ratio)",
                    "strategy_selection"
                ).with_units("ratio")
            ),
        }
    }
}

impl Default for StrategyThresholds {
    fn default() -> Self {
        Self::make_default()
    }
}
