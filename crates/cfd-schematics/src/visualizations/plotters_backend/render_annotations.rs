use crate::domain::model::NetworkBlueprint;
use crate::error::{VisualizationError, VisualizationResult};
use crate::visualizations::annotations::{should_render_label, MarkerRole, SchematicAnnotations};
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;

use super::convert_color;

pub(super) fn draw_annotation_overlay<DB: DrawingBackend>(
    chart: &mut ChartContext<'_, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    system: &NetworkBlueprint,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    let label_dx = system.box_dims.0 * 0.008;
    let label_dy = system.box_dims.1 * 0.008;

    for marker in &annotations.markers {
        let marker_color = annotations.style.color_for_role(marker.role);
        let style = convert_color(&marker_color).filled();

        chart
            .draw_series(std::iter::once(Circle::new(
                marker.point,
                annotations.style.dot_radius_px,
                style,
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        if should_render_label(marker, annotations.label_density) {
            let label_text = marker.label.as_deref().unwrap_or_default();
            chart
                .draw_series(std::iter::once(Text::new(
                    label_text.to_string(),
                    (marker.point.0 + label_dx, marker.point.1 + label_dy),
                    ("sans-serif", annotations.style.label_font_size_pt.max(8))
                        .into_font()
                        .color(&convert_color(&marker_color)),
                )))
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }
    }

    if annotations.style.show_legend {
        draw_annotation_legend(chart, system, annotations)?;
    }

    Ok(())
}

fn draw_annotation_legend<DB: DrawingBackend>(
    chart: &mut ChartContext<'_, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    system: &NetworkBlueprint,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    use std::collections::BTreeSet;

    let present_roles: BTreeSet<MarkerRole> = annotations.markers.iter().map(|m| m.role).collect();
    if present_roles.is_empty() && annotations.legend_note.is_none() {
        return Ok(());
    }

    let legend_x = system.box_dims.0 * 0.02;
    let mut legend_y = system.box_dims.1 * 0.95;
    let row_step = system.box_dims.1 * 0.035;

    if !present_roles.is_empty() {
        chart
            .draw_series(std::iter::once(Text::new(
                "Legend".to_string(),
                (legend_x, legend_y),
                ("sans-serif", 12).into_font().color(&RGBColor(20, 20, 20)),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        legend_y -= row_step;
    }

    for role in present_roles {
        let color = annotations.style.color_for_role(role);
        chart
            .draw_series(std::iter::once(Circle::new(
                (legend_x, legend_y),
                (annotations.style.dot_radius_px - 1).max(3),
                convert_color(&color).filled(),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        chart
            .draw_series(std::iter::once(Text::new(
                role.legend_label().to_string(),
                (legend_x + system.box_dims.0 * 0.018, legend_y),
                ("sans-serif", 10).into_font().color(&RGBColor(50, 50, 50)),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        legend_y -= row_step;
    }

    if let Some(note) = annotations.legend_note.as_deref() {
        chart
            .draw_series(std::iter::once(Text::new(
                note.to_string(),
                (legend_x, legend_y),
                ("sans-serif", 10).into_font().color(&RGBColor(30, 30, 30)),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
    }

    Ok(())
}
