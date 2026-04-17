use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;

use crate::domain::model::NetworkBlueprint;
use crate::error::VisualizationResult;
use crate::visualizations::analysis_field::{colorize, AnalysisOverlay};

use super::convert_color;

pub(super) fn draw_colorbar<DB: DrawingBackend>(
    chart: &mut ChartContext<'_, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    system: &NetworkBlueprint,
    x_buffer: f64,
    y_buffer: f64,
    overlay: &AnalysisOverlay,
) -> VisualizationResult<()> {
    let (length, width) = system.box_dims;
    let has_edge_data = !overlay.edge_data.is_empty();

    let (min_val, max_val) = if has_edge_data {
        overlay.edge_range()
    } else {
        overlay.node_range()
    };

    let bar_steps = 20usize;
    let bar_x_start = length + x_buffer * 0.2;
    let bar_x_end = length + x_buffer * 0.7;
    let bar_y_start = y_buffer;
    let bar_y_end = width - y_buffer;
    let step_height = (bar_y_end - bar_y_start) / bar_steps as f64;

    for idx in 0..bar_steps {
        let t = 1.0 - (idx as f64) / (bar_steps as f64);
        let color = colorize(t, overlay.colormap);
        let y_top = bar_y_start + (idx as f64) * step_height;
        let y_bot = y_top + step_height;
        let rect_pts = vec![
            (bar_x_start, y_top),
            (bar_x_end, y_top),
            (bar_x_end, y_bot),
            (bar_x_start, y_bot),
            (bar_x_start, y_top),
        ];
        chart
            .draw_series(std::iter::once(PathElement::new(
                rect_pts,
                convert_color(&color).stroke_width(4),
            )))
            .map_err(|e| crate::error::VisualizationError::rendering_error(&e.to_string()))?;
    }

    chart
        .draw_series(std::iter::once(Text::new(
            format!("{max_val:.2e}"),
            (bar_x_start, bar_y_start),
            ("sans-serif", 10).into_font(),
        )))
        .map_err(|e| crate::error::VisualizationError::rendering_error(&e.to_string()))?;

    chart
        .draw_series(std::iter::once(Text::new(
            format!("{min_val:.2e}"),
            (bar_x_start, bar_y_end),
            ("sans-serif", 10).into_font(),
        )))
        .map_err(|e| crate::error::VisualizationError::rendering_error(&e.to_string()))?;

    chart
        .draw_series(std::iter::once(Text::new(
            overlay.field.label().to_string(),
            (bar_x_start, f64::midpoint(bar_y_start, bar_y_end)),
            ("sans-serif", 9).into_font(),
        )))
        .map_err(|e| crate::error::VisualizationError::rendering_error(&e.to_string()))?;

    Ok(())
}
