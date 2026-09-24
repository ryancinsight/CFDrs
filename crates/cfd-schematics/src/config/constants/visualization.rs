//! Visualization constants for chart rendering

use super::primitives;

/// Visualization constants previously hardcoded
///
/// The values themselves live in [`primitives`].
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
}

impl VisualizationConstants {
    /// Canonical default visualization constants.
    pub const DEFAULT: Self = Self {
        default_chart_margin: primitives::DEFAULT_CHART_MARGIN,
        default_chart_right_margin: primitives::DEFAULT_CHART_RIGHT_MARGIN,
        default_x_label_area_size: primitives::DEFAULT_X_LABEL_AREA_SIZE,
        default_y_label_area_size: primitives::DEFAULT_Y_LABEL_AREA_SIZE,
    };
}

impl Default for VisualizationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
