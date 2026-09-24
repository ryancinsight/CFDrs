//! Visualization constants for chart rendering

/// Visualization constants previously hardcoded
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VisualizationConstants {
    /// Default margin for chart rendering
    pub default_chart_margin: u32,

    /// Default right margin for chart rendering
    pub default_chart_right_margin: u32,

    /// Default label area size for x-axis
    pub default_x_label_area_size: u32,

    /// Default label area size for y-axis
    pub default_y_label_area_size: u32,

    /// Default buffer factor for chart boundaries
    pub default_boundary_buffer_factor: f64,
}

impl VisualizationConstants {
    /// Canonical default visualization constants.
    pub const DEFAULT: Self = Self {
        default_chart_margin: 20,
        default_chart_right_margin: 150,
        default_x_label_area_size: 30,
        default_y_label_area_size: 30,
        default_boundary_buffer_factor: 0.1,
    };
}

impl Default for VisualizationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
