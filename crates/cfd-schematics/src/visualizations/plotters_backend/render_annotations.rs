use crate::domain::model::NetworkBlueprint;
use crate::error::{VisualizationError, VisualizationResult};
use crate::visualizations::annotations::{
    should_render_label, AnnotationMarker, MarkerRole, SchematicAnnotations,
};
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;

use super::convert_color;

pub(super) fn draw_annotation_overlay<DB: DrawingBackend>(
    chart: &mut ChartContext<'_, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    system: &NetworkBlueprint,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
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
            let (label_dx, label_dy) = label_offset_for_marker(marker, system);
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

fn label_offset_for_marker(marker: &AnnotationMarker, system: &NetworkBlueprint) -> (f64, f64) {
    let width = system.box_dims.0.max(1.0);
    let height = system.box_dims.1.max(1.0);
    let label_dx = width * 0.020;
    let label_dy = height * 0.018;
    let center_x = width * 0.5;
    let center_y = height * 0.5;

    let horizontal = match marker.role {
        MarkerRole::Inlet => label_dx,
        MarkerRole::Outlet => -label_dx,
        _ if marker.point.0 < center_x => label_dx,
        _ => -label_dx,
    };

    let vertical_scale = match marker.role {
        MarkerRole::Inlet | MarkerRole::Outlet => 0.55,
        MarkerRole::VenturiThroat | MarkerRole::TherapyTarget => 1.25,
        MarkerRole::Split | MarkerRole::Merge => 1.15,
        MarkerRole::Bypass | MarkerRole::Internal => 0.95,
    };
    let vertical = if marker.point.1 <= center_y {
        -label_dy * vertical_scale
    } else {
        label_dy * vertical_scale
    };

    (horizontal, vertical)
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

#[cfg(test)]
mod tests {
    use super::label_offset_for_marker;
    use crate::domain::model::{NetworkBlueprint, NodeId, NodeKind, NodeSpec};
    use crate::visualizations::annotations::{AnnotationMarker, MarkerRole};

    fn synthetic_system() -> NetworkBlueprint {
        NetworkBlueprint {
            name: "test".to_string(),
            box_dims: (127.76, 85.47),
            nodes: vec![NodeSpec {
                id: NodeId::new("0"),
                kind: NodeKind::Inlet,
                point: (0.0, 42.735),
                layout: None,
                junction_geometry: None,
                metadata: None,
            }],
            channels: Vec::new(),
            box_outline: vec![
                ((0.0, 0.0), (127.76, 0.0)),
                ((127.76, 0.0), (127.76, 85.47)),
                ((127.76, 85.47), (0.0, 85.47)),
                ((0.0, 85.47), (0.0, 0.0)),
            ],
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: false,
        }
    }

    #[test]
    fn label_offsets_move_edge_labels_outward_and_center_labels_off_the_lane() {
        let system = synthetic_system();
        let inlet = AnnotationMarker::new((0.0, 42.735), MarkerRole::Inlet).with_label("IN1", true);
        let outlet =
            AnnotationMarker::new((127.76, 42.735), MarkerRole::Outlet).with_label("OUT1", true);
        let upper_split =
            AnnotationMarker::new((32.0, 18.0), MarkerRole::Split).with_label("S1", true);
        let lower_split =
            AnnotationMarker::new((32.0, 68.0), MarkerRole::Split).with_label("S2", true);
        let venturi = AnnotationMarker::new((60.0, 42.735), MarkerRole::VenturiThroat)
            .with_label("TH1", true);

        let (inlet_dx, inlet_dy) = label_offset_for_marker(&inlet, &system);
        let (outlet_dx, outlet_dy) = label_offset_for_marker(&outlet, &system);
        let (upper_dx, upper_dy) = label_offset_for_marker(&upper_split, &system);
        let (lower_dx, lower_dy) = label_offset_for_marker(&lower_split, &system);
        let (venturi_dx, venturi_dy) = label_offset_for_marker(&venturi, &system);

        assert!(
            inlet_dx > 0.0,
            "inlet labels should be shifted into the figure"
        );
        assert!(
            outlet_dx < 0.0,
            "outlet labels should be shifted away from the outlet edge"
        );
        assert!(
            inlet_dy < 0.0,
            "inlet labels should sit slightly above the node"
        );
        assert!(
            outlet_dy < 0.0,
            "outlet labels should sit slightly above the node"
        );
        assert!(
            upper_dy < 0.0,
            "upper split labels should move above the main lane"
        );
        assert!(
            lower_dy > 0.0,
            "lower split labels should move below the main lane"
        );
        assert!(
            venturi_dy.abs() > 0.0,
            "venturi labels need a vertical offset"
        );
        assert!(
            upper_dx > 0.0,
            "left-half split labels should remain to the right of the node"
        );
        assert!(
            lower_dx > 0.0,
            "left-half split labels should remain to the right of the node"
        );
        assert!(
            venturi_dx > 0.0,
            "centerline venturi labels should avoid starting on the marker"
        );
    }
}
