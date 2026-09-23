use super::path_geometry::{path_intersects_any, polyline_length_mm};
use super::routing::route_monotone_treatment_path;
use super::{
    CenterSerpentinePathSpec, PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest,
    create_primitive_selective_tree_geometry,
};
use aequitas::systems::si::quantities::Length;

#[test]
fn primitive_selective_tree_annotation_preserves_positive_channel_lengths() {
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "primitive-selective-lengths".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(8.0e-3),
        throat_width_m: Length::from_base(55.0e-6),
        throat_length_m: Length::from_base(110.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: false,
        treatment_branch_throat_count: 1,
        center_serpentine: None,
    });

    let invalid_lengths: Vec<String> = blueprint
        .channels
        .iter()
        .filter(|channel| {
            let length_m = channel.length_m.into_base();
            !(length_m.is_finite() && length_m > 0.0)
        })
        .map(|channel| format!("{}={}", channel.id.as_str(), channel.length_m.into_base()))
        .collect();

    assert!(
        invalid_lengths.is_empty(),
        "primitive selective annotation should preserve positive channel lengths: {}",
        invalid_lengths.join(", ")
    );
}

#[test]
fn primitive_selective_annotation_lengths_match_materialized_paths() {
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "primitive-selective-path-lengths".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(4.0e-3),
        throat_width_m: Length::from_base(55.0e-6),
        throat_length_m: Length::from_base(110.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: false,
        treatment_branch_throat_count: 1,
        center_serpentine: Some(CenterSerpentinePathSpec {
            segments: 5,
            bend_radius_m: Length::from_base(3.0e-3),
            wave_type: crate::SerpentineWaveType::default(),
        }),
    });

    let mismatches: Vec<String> = blueprint
        .channels
        .iter()
        .filter(|channel| channel.path.len() >= 2)
        .filter_map(|channel| {
            let expected = polyline_length_mm(&channel.path) * 1.0e-3;
            (channel.length_m.into_base() != expected).then(|| {
                format!(
                    "{}: expected {expected}, got {}",
                    channel.id.as_str(),
                    channel.length_m.into_base()
                )
            })
        })
        .collect();

    assert!(
        mismatches.is_empty(),
        "annotated channel lengths must match their materialized paths: {}",
        mismatches.join(", ")
    );
}

#[test]
fn primitive_selective_venturi_paths_have_no_unresolved_crossings() {
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "primitive-selective-no-crossings".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(8.0e-3),
        throat_width_m: Length::from_base(55.0e-6),
        throat_length_m: Length::from_base(110.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: true,
        treatment_branch_throat_count: 1,
        center_serpentine: None,
    });

    assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
    blueprint
        .validate()
        .expect("venturi treatment paths must remain planar");

    for channel in blueprint.venturi_channels() {
        if let (Some(start), Some(end)) = (channel.path.first(), channel.path.last())
            && (start.1 - end.1).abs() < 1e-9
        {
            assert!(
                channel
                    .path
                    .iter()
                    .all(|point| (point.1 - start.1).abs() < 1e-9),
                "equal-y treatment channel {} must stay on its own lane",
                channel.id.as_str()
            );
        }
    }
}

#[test]
fn primitive_selective_penta_quad_tri_is_geometry_authorable() {
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "primitive-selective-penta-quad-tri".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Penta,
            PrimitiveSelectiveSplitKind::Quad,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(4.0e-3),
        throat_width_m: Length::from_base(35.0e-6),
        throat_length_m: Length::from_base(180.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.45,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: true,
        treatment_branch_throat_count: 1,
        center_serpentine: None,
    });

    let envelope_bound = 8.0 * f64::EPSILON * 128.0;
    assert!((blueprint.box_dims.0 - 127.76).abs() <= envelope_bound);
    assert!((blueprint.box_dims.1 - 85.47).abs() <= envelope_bound);
    assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
    blueprint
        .validate()
        .expect("Penta->Quad->Tri primitive selective tree should remain planar");
}

#[test]
fn primitive_selective_tri_tri_retains_multiple_treatment_window_lanes() {
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "primitive-selective-tritri-lanes".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(8.0e-3),
        throat_width_m: Length::from_base(55.0e-6),
        throat_length_m: Length::from_base(110.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: false,
        treatment_branch_throat_count: 1,
        center_serpentine: None,
    });

    let mid_x = blueprint.box_dims.0 * 0.5;
    let mut lane_keys = std::collections::BTreeSet::new();
    for channel in &blueprint.channels {
        if channel.therapy_zone
            != Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
        {
            continue;
        }
        let points = if channel.path.is_empty() {
            let start = blueprint
                .nodes
                .iter()
                .find(|node| node.id == channel.from)
                .map(|node| node.point);
            let end = blueprint
                .nodes
                .iter()
                .find(|node| node.id == channel.to)
                .map(|node| node.point);
            [start, end].into_iter().flatten().collect::<Vec<_>>()
        } else {
            channel.path.clone()
        };
        if points.is_empty() {
            continue;
        }
        let min_x = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
        let max_x = points
            .iter()
            .map(|(x, _)| *x)
            .fold(f64::NEG_INFINITY, f64::max);
        if max_x <= mid_x + 1.0e-6 || min_x < mid_x - 1.0e-6 {
            continue;
        }
        let mean_y = points.iter().map(|(_, y)| *y).sum::<f64>() / points.len() as f64;
        lane_keys.insert((mean_y * 100.0).round() as i64);
    }

    assert!(
        lane_keys.len() >= 3,
        "Tri->Tri selective trees must preserve multiple treatment-window lanes instead of collapsing to a centerline surrogate, got {lane_keys:?}"
    );
}

#[test]
fn monotone_treatment_routing_preserves_equal_y_lane() {
    let routed = route_monotone_treatment_path(
        Some((10.0, 24.0)),
        Some((40.0, 24.0)),
        24.0,
        42.0,
        |_| false,
    );
    assert_eq!(routed, vec![(10.0, 24.0), (40.0, 24.0)]);
}

#[test]
fn monotone_treatment_routing_doglegs_around_existing_branch() {
    let existing_paths = vec![vec![(3.0, 0.0), (7.0, 4.0)]];
    let routed = route_monotone_treatment_path(
        Some((0.0, 2.0)),
        Some((10.0, 2.0)),
        0.0,
        2.0,
        |candidate| path_intersects_any(candidate, existing_paths.iter()),
    );

    assert!(
        routed.len() > 2,
        "crossing direct path should reroute with a dogleg"
    );
    assert!(
        !path_intersects_any(&routed, &existing_paths),
        "rerouted path must not cross the existing treatment lane"
    );
}

#[test]
fn monotone_treatment_routing_falls_back_to_direct_segment() {
    let routed =
        route_monotone_treatment_path(Some((0.0, 2.0)), Some((10.0, 2.0)), 0.0, 2.0, |_| true);

    assert_eq!(routed, vec![(0.0, 2.0), (10.0, 2.0)]);
}

#[test]
fn primitive_selective_serpentine_mirrors_on_both_sides_of_midline() {
    use crate::domain::model::ChannelShape;
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "serp-mirror-check".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(4.0e-3),
        throat_width_m: Length::from_base(55.0e-6),
        throat_length_m: Length::from_base(110.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: false,
        treatment_branch_throat_count: 1,
        center_serpentine: Some(CenterSerpentinePathSpec {
            segments: 5,
            bend_radius_m: Length::from_base(3.0e-3),
            wave_type: crate::SerpentineWaveType::default(),
        }),
    });

    let mid_x = blueprint.box_dims.0 * 0.5;
    let serpentine_channels: Vec<_> = blueprint
        .channels
        .iter()
        .filter(|ch| matches!(ch.channel_shape, ChannelShape::Serpentine { .. }))
        .collect();

    assert!(
        serpentine_channels.len() >= 2,
        "at least one serpentine channel on each side expected, got {}",
        serpentine_channels.len()
    );

    let node_pts: std::collections::HashMap<_, _> = blueprint
        .nodes
        .iter()
        .map(|n| (n.id.clone(), n.point))
        .collect();
    let left_count = serpentine_channels
        .iter()
        .filter(|ch| {
            let max_x = [
                node_pts.get(&ch.from).map(|p| p.0),
                node_pts.get(&ch.to).map(|p| p.0),
            ]
            .into_iter()
            .flatten()
            .fold(f64::NEG_INFINITY, f64::max);
            max_x <= mid_x + 1e-6
        })
        .count();
    let right_count = serpentine_channels
        .iter()
        .filter(|ch| {
            let min_x = [
                node_pts.get(&ch.from).map(|p| p.0),
                node_pts.get(&ch.to).map(|p| p.0),
            ]
            .into_iter()
            .flatten()
            .fold(f64::INFINITY, f64::min);
            min_x >= mid_x - 1e-6
        })
        .count();

    assert!(
        left_count >= 1,
        "split-side (left) treatment channels must have serpentine overlays, got {left_count}"
    );
    assert!(
        right_count >= 1,
        "merge-side (right) treatment channels must have serpentine overlays, got {right_count}"
    );
}

#[test]
fn treatment_channels_retain_center_treatment_role_on_merge_side() {
    use crate::geometry::metadata::ChannelVisualRole;
    let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "role-mirror-check".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(4.0e-3),
        throat_width_m: Length::from_base(55.0e-6),
        throat_length_m: Length::from_base(110.0e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.45,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: false,
        treatment_branch_throat_count: 1,
        center_serpentine: Some(CenterSerpentinePathSpec {
            segments: 5,
            bend_radius_m: Length::from_base(3.0e-3),
            wave_type: crate::SerpentineWaveType::default(),
        }),
    });

    let mid_x = blueprint.box_dims.0 * 0.5;
    let node_pts: std::collections::HashMap<_, _> = blueprint
        .nodes
        .iter()
        .map(|n| (n.id.clone(), n.point))
        .collect();

    let merge_side_treatment: Vec<_> = blueprint
        .channels
        .iter()
        .filter(|ch| {
            ch.therapy_zone == Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
        })
        .filter(|ch| {
            let max_x = [
                node_pts.get(&ch.from).map(|p| p.0),
                node_pts.get(&ch.to).map(|p| p.0),
            ]
            .into_iter()
            .flatten()
            .fold(f64::NEG_INFINITY, f64::max);
            max_x > mid_x + 1e-6
        })
        .collect();

    assert!(
        !merge_side_treatment.is_empty(),
        "should have treatment channels on the merge side"
    );
    for ch in &merge_side_treatment {
        assert_eq!(
            ch.visual_role,
            Some(ChannelVisualRole::CenterTreatment),
            "merge-side treatment channel {} must retain CenterTreatment role",
            ch.id.as_str()
        );
    }
}
