use cfd_schematics::domain::model::ChannelShape;
use cfd_schematics::geometry::metadata::JunctionFamily;
use cfd_schematics::interface::presets::{
    cascade_center_trifurcation_rect, double_trifurcation_cif_venturi_rect,
    incremental_filtration_tri_bi_rect_staged_remerge, CenterSerpentineSpec,
};
use cfd_schematics::visualizations::throat_count_from_blueprint_metadata;
use cfd_schematics::NetworkBlueprint;

fn node_has_point(bp: &NetworkBlueprint, node_id: &str) -> bool {
    bp.nodes.iter().any(|node| node.id.0 == node_id)
}

fn channel_path(bp: &NetworkBlueprint, channel_id: &str) -> Vec<(f64, f64)> {
    bp.channels
        .iter()
        .find(|channel| channel.id.0 == channel_id)
        .map(|channel| channel.path.clone())
        .expect("channel path must exist")
}

fn node_y(bp: &NetworkBlueprint, node_id: &str) -> f64 {
    bp.nodes
        .iter()
        .find(|node| node.id.0 == node_id)
        .map(|node| node.point.1)
        .expect("node must exist")
}

#[test]
fn staged_cif_blueprint_uses_native_geometry() {
    let bp = incremental_filtration_tri_bi_rect_staged_remerge(
        "cif-render",
        12e-3,
        8e-3,
        6e-3,
        2,
        2.0e-3,
        0.45,
        0.55,
        0.68,
        100e-6,
        300e-6,
        1.5e-3,
        1.0e-3,
        true,
        None,
    );

    assert!(node_has_point(&bp, "split_lv0"));
    assert!(node_has_point(&bp, "hy_tri"));

    let throat = bp
        .channels
        .iter()
        .find(|channel| channel.venturi_geometry.is_some())
        .expect("venturi blueprint should render a throat channel");
    match &throat.cross_section {
        cfd_schematics::domain::model::CrossSectionSpec::Rectangular { width_m, .. } => {
            assert!(*width_m > 0.0);
        }
        _ => panic!("Expected rectangular"),
    }
    assert_eq!(throat_count_from_blueprint_metadata(&bp), 1);
}

#[test]
fn staged_cif_layout_keeps_upper_lower_bypass_paths_mirrored() {
    let bp = incremental_filtration_tri_bi_rect_staged_remerge(
        "cif-symmetric",
        12e-3,
        8e-3,
        6e-3,
        2,
        2.0e-3,
        0.45,
        0.55,
        0.68,
        100e-6,
        300e-6,
        1.5e-3,
        1.0e-3,
        true,
        None,
    );

    let upper = channel_path(&bp, "L_lv0");
    let lower = channel_path(&bp, "R_lv0");
    let y_mid = node_y(&bp, "split_lv0");

    assert!(upper.len() >= 4);
    assert!(lower.len() >= 4);
    assert!((upper[1].0 - lower[1].0).abs() < 1e-9);
    assert!(((upper[1].1 + lower[1].1) - 2.0 * y_mid).abs() < 1e-9);
}

#[test]
fn blueprint_render_adds_wall_port_stubs_for_inlet_and_outlet() {
    let bp = incremental_filtration_tri_bi_rect_staged_remerge(
        "cif-ports",
        12e-3,
        8e-3,
        6e-3,
        2,
        2.0e-3,
        0.45,
        0.55,
        0.68,
        100e-6,
        300e-6,
        1.5e-3,
        1.0e-3,
        false,
        None,
    );

    let inlet = bp
        .nodes
        .iter()
        .find(|node| node.id.0 == "inlet")
        .expect("inlet node must exist");
    let outlet = bp
        .nodes
        .iter()
        .find(|node| node.id.0 == "outlet")
        .expect("outlet node must exist");

    assert_eq!(inlet.point.0, 0.0);
    assert!(outlet.point.0 > 10.0); // Right edge side
}

#[test]
fn cct_blueprint_render_preserves_center_serpentine_lane() {
    let bp = cascade_center_trifurcation_rect(
        "cct-serp",
        12e-3,
        8e-3,
        3,
        2.0e-3,
        0.45,
        100e-6,
        300e-6,
        1.0e-3,
        false,
        Some(CenterSerpentineSpec {
            segments: 5,
            bend_radius_m: 3.0e-3,
            segment_length_m: 7.5e-3,
        }),
    );

    let center_lane = bp
        .channels
        .iter()
        .find(|channel| matches!(channel.channel_shape, ChannelShape::Serpentine { .. }))
        .expect("center treatment lane should stay serpentine in render path");

    assert!(center_lane.path.len() >= 2);
}

#[test]
fn cct_serpentine_lane_materializes_mirrored_rotated_s_offsets() {
    let bp = cascade_center_trifurcation_rect(
        "cct-serp-mirror",
        12e-3,
        8e-3,
        3,
        2.0e-3,
        0.45,
        100e-6,
        300e-6,
        1.0e-3,
        false,
        Some(CenterSerpentineSpec {
            segments: 3,
            bend_radius_m: 3.0e-3,
            segment_length_m: 7.5e-3,
        }),
    );

    let center_lane = bp
        .channels
        .iter()
        .find(|channel| matches!(channel.channel_shape, ChannelShape::Serpentine { .. }))
        .expect("center treatment lane should stay serpentine");
    let start = center_lane.path.first().copied().expect("path start");
    let end = center_lane.path.last().copied().expect("path end");
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    assert!(length > 0.0, "serpentine chord length must be non-zero");
    let nx = -dy / length;
    let ny = dx / length;
    let signed_offsets: Vec<f64> = center_lane
        .path
        .iter()
        .skip(1)
        .take(center_lane.path.len().saturating_sub(2))
        .map(|point| ((point.0 - start.0) * nx) + ((point.1 - start.1) * ny))
        .filter(|offset| offset.abs() > 1.0e-6)
        .collect();

    assert!(
        signed_offsets.iter().any(|offset| *offset > 0.0),
        "serpentine path must lobe to one side of the centerline"
    );
    assert!(
        signed_offsets.iter().any(|offset| *offset < 0.0),
        "serpentine path must mirror to the opposite side of the centerline"
    );
}

#[test]
fn dtcv_blueprint_render_builds_four_true_trifurcation_nodes() {
    let bp = double_trifurcation_cif_venturi_rect(
        "dtcv-full-tree",
        12e-3,
        8e-3,
        5.0e-3,
        0.54,
        0.45,
        100e-6,
        300e-6,
        1.2e-3,
        2,
        1.0e-3,
    );

    assert!(node_has_point(&bp, "split1"));
    assert!(node_has_point(&bp, "split2_upper"));
    assert!(node_has_point(&bp, "split2"));
    assert!(node_has_point(&bp, "split2_lower"));

    let tri_nodes = bp
        .nodes
        .iter()
        .filter(|node| {
            node.junction_geometry
                .as_ref()
                .is_some_and(|g| g.junction_family == JunctionFamily::Trifurcation)
        })
        .count();
    assert_eq!(tri_nodes, 4);
}

#[test]
fn dtcv_acoustic_render_keeps_full_tree_without_frusta() {
    let bp = double_trifurcation_cif_venturi_rect(
        "dtcv-acoustic",
        12e-3,
        8e-3,
        5.0e-3,
        0.54,
        0.45,
        100e-6,
        300e-6,
        1.2e-3,
        0,
        1.0e-3,
    );

    assert!(node_has_point(&bp, "split2_upper"));
    assert!(channel_path(&bp, "upper_L").len() >= 2);
    assert!(!bp
        .channels
        .iter()
        .any(|channel| channel.venturi_geometry.is_some()));
}
