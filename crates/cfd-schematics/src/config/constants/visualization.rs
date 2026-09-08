//! Visualization constants for chart rendering

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

/// Visualization constants previously hardcoded
pub struct VisualizationConstants {
    /// Default margin for chart rendering
    pub default_chart_margin: ConfigurableParameter<u32>,

    /// Default right margin for chart rendering
    pub default_chart_right_margin: ConfigurableParameter<u32>,

    /// Default label area size for x-axis
    pub default_x_label_area_size: ConfigurableParameter<u32>,

    /// Default label area size for y-axis
    pub default_y_label_area_size: ConfigurableParameter<u32>,

    /// Default buffer factor for chart boundaries
    pub default_boundary_buffer_factor: ConfigurableParameter<f64>,
}

impl VisualizationConstants {
    /// Create default visualization constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            default_chart_margin: ConfigurableParameter::new(
                20u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(5, 100),
                ]),
                ParameterMetadata::new(
                    "default_chart_margin",
                    "Default margin for chart rendering",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_chart_right_margin: ConfigurableParameter::new(
                150u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(50, 500),
                ]),
                ParameterMetadata::new(
                    "default_chart_right_margin",
                    "Default right margin for chart rendering",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_x_label_area_size: ConfigurableParameter::new(
                30u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(10, 100),
                ]),
                ParameterMetadata::new(
                    "default_x_label_area_size",
                    "Default label area size for x-axis",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_y_label_area_size: ConfigurableParameter::new(
                30u32,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<u32>::positive(),
                    ParameterConstraints::range(10, 100),
                ]),
                ParameterMetadata::new(
                    "default_y_label_area_size",
                    "Default label area size for y-axis",
                    "visualization",
                )
                .with_units("pixels"),
            ),

            default_boundary_buffer_factor: ConfigurableParameter::new(
                0.1,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.01, 0.5),
                ]),
                ParameterMetadata::new(
                    "default_boundary_buffer_factor",
                    "Default buffer factor for chart boundaries",
                    "visualization",
                )
                .with_units("ratio"),
            ),
        }
    }
}

impl Default for VisualizationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}
