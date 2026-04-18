use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::BlueprintRenderHints;
use crate::visualizations::annotations::{
    classify_node_roles, project_markers_along_path, throat_count_from_blueprint_metadata,
    venturi_marker_points_from_blueprint, AnnotationMarker, MarkerRole, SchematicAnnotations,
};

pub(super) fn build_auto_annotations(blueprint: &NetworkBlueprint) -> SchematicAnnotations {
    let mut annotations = SchematicAnnotations::report_default();
    let mut roles = classify_node_roles(blueprint);

    let has_target =
        blueprint.length_in_zone(crate::domain::therapy_metadata::TherapyZone::CancerTarget) > 0.0;
    let has_bypass =
        blueprint.length_in_zone(crate::domain::therapy_metadata::TherapyZone::HealthyBypass) > 0.0;

    if has_target || has_bypass {
        let box_dims = blueprint.box_dims;
        let y_mid = box_dims.1 * 0.5;
        let tx_min = box_dims.0 * 0.35;
        let tx_max = box_dims.0 * 0.65;
        for (node_idx, node) in blueprint.nodes.iter().enumerate() {
            let base = roles
                .get(&node_idx)
                .copied()
                .unwrap_or(MarkerRole::Internal);
            if matches!(
                base,
                MarkerRole::Inlet | MarkerRole::Outlet | MarkerRole::Split | MarkerRole::Merge
            ) {
                continue;
            }
            let near_cl = (node.point.1 - y_mid).abs() <= box_dims.1 * 0.18;
            let in_win = node.point.0 >= tx_min && node.point.0 <= tx_max;
            if has_target && near_cl && in_win {
                roles.insert(node_idx, MarkerRole::TherapyTarget);
            } else if has_bypass && !near_cl {
                roles.insert(node_idx, MarkerRole::Bypass);
            }
        }
    }

    let mut split_idx = 1usize;
    let mut merge_idx = 1usize;
    let mut inlet_idx = 1usize;
    let mut outlet_idx = 1usize;
    for (node_idx, node) in blueprint.nodes.iter().enumerate() {
        let role = roles
            .get(&node_idx)
            .copied()
            .unwrap_or(MarkerRole::Internal);
        let marker = match role {
            MarkerRole::Inlet => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("IN{inlet_idx}"), true);
                inlet_idx += 1;
                m
            }
            MarkerRole::Outlet => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("OUT{outlet_idx}"), true);
                outlet_idx += 1;
                m
            }
            MarkerRole::Split => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("S{split_idx}"), true);
                split_idx += 1;
                m
            }
            MarkerRole::Merge => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("M{merge_idx}"), true);
                merge_idx += 1;
                m
            }
            _ => AnnotationMarker::new(node.point, role),
        };
        annotations.markers.push(marker);
    }

    let hints: Option<&BlueprintRenderHints> = blueprint.render_hints();
    let mut throat_count = throat_count_from_blueprint_metadata(blueprint);
    if throat_count == 0 {
        if let Some(h) = hints {
            throat_count = h.throat_count_hint;
        }
    }

    if throat_count > 0 {
        let actual_points = venturi_marker_points_from_blueprint(blueprint);
        let mut th_idx = 1usize;
        if actual_points.is_empty() {
            let main_path = crate::visualizations::annotations::center_biased_main_path(blueprint);
            let box_dims = blueprint.box_dims;
            let zone = (box_dims.0 * 0.35, box_dims.0 * 0.65);
            for point in project_markers_along_path(&main_path, throat_count, zone) {
                annotations.markers.push(
                    AnnotationMarker::new(point, MarkerRole::VenturiThroat)
                        .with_label(format!("TH{th_idx}"), true),
                );
                th_idx += 1;
            }
        } else {
            for point in actual_points {
                annotations.markers.push(
                    AnnotationMarker::new(point, MarkerRole::VenturiThroat)
                        .with_label(format!("TH{th_idx}"), true),
                );
                th_idx += 1;
            }
        }
    }

    if let Some(h) = hints {
        let min_w = blueprint
            .channels
            .iter()
            .map(|ch| ch.cross_section.dims().0)
            .fold(f64::INFINITY, f64::min);
        let max_w = blueprint
            .channels
            .iter()
            .map(|ch| ch.cross_section.dims().0)
            .fold(f64::NEG_INFINITY, f64::max);
        annotations.legend_note = Some(format!(
            "{}  |  seq: {}  |  layers: {}  |  throats: {}  |  width: {:.2}-{:.2} mm  |  line thickness ∝ channel width",
            h.treatment_label,
            h.stage_sequence,
            h.split_layers,
            throat_count,
            min_w * 1e3,
            max_w * 1e3
        ));
        annotations.markers.push(
            AnnotationMarker::new(
                (blueprint.box_dims.0 * 0.50, blueprint.box_dims.1 * 0.82),
                MarkerRole::Internal,
            )
            .with_label(
                format!(
                    "{}  |  {} split layers  |  {} treatment",
                    h.stage_sequence, h.split_layers, h.treatment_label
                ),
                true,
            ),
        );
    }

    let volume_note = blueprint.fluid_volume_summary().display_label;
    let show_volume_marker = hints.is_none();
    let y_fraction = 0.88;
    annotations.legend_note = Some(match annotations.legend_note.take() {
        Some(note) => format!("{note} | {volume_note}"),
        None => volume_note.clone(),
    });
    if show_volume_marker {
        annotations.markers.push(
            AnnotationMarker::new(
                (
                    blueprint.box_dims.0 * 0.50,
                    blueprint.box_dims.1 * y_fraction,
                ),
                MarkerRole::Internal,
            )
            .with_label(volume_note, true),
        );
    }

    annotations
}
