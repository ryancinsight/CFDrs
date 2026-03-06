use cfd_schematics::geometry::{
    metadata::{ChannelPathMetadata, JunctionFamily, JunctionGeometryMetadata, NodeLayoutMetadata},
    ChannelType,
};
use cfd_schematics::interface::presets::{
    cascade_center_trifurcation_rect, double_trifurcation_cif_venturi_rect,
    incremental_filtration_tri_bi_rect_staged_remerge, CenterSerpentineSpec,
};
use cfd_schematics::visualizations::throat_count_from_blueprint_metadata;
use cfd_schematics::{channel_system_from_blueprint, NetworkBlueprint};

fn node_has_layout(bp: &NetworkBlueprint, node_id: &str) -> bool {
    bp.nodes
        .iter()
        .find(|node| node.id.as_str() == node_id)
        .and_then(|node| {
            node.layout.or_else(|| {
                node.metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<NodeLayoutMetadata>())
                    .copied()
            })
        })
        .is_some()
}

fn channel_path(bp: &NetworkBlueprint, channel_id: &str) -> Vec<(f64, f64)> {
    bp.channels
        .iter()
        .find(|channel| channel.id.as_str() == channel_id)
        .and_then(|channel| {
            channel.path.clone().or_else(|| {
                channel
                    .metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<ChannelPathMetadata>())
                    .cloned()
            })
        })
        .map(|path| path.polyline_mm.clone())
        .expect("channel path metadata must exist")
}

fn node_y(bp: &NetworkBlueprint, node_id: &str) -> f64 {
    bp.nodes
        .iter()
        .find(|node| node.id.as_str() == node_id)
        .and_then(|node| {
            node.layout.or_else(|| {
                node.metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<NodeLayoutMetadata>())
                    .copied()
            })
        })
        .map(|layout| layout.y_mm)
        .expect("node layout metadata must exist")
}

#[test]
fn staged_cif_blueprint_render_uses_layout_metadata_and_throat_geometry() {
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

    assert!(node_has_layout(&bp, "split_lv0"));
    assert!(node_has_layout(&bp, "hy_tri"));

    let system = channel_system_from_blueprint(&bp, Some((127.76, 85.47)))
        .expect("blueprint-native schematic conversion must succeed");
    let throat = system
        .channels
        .iter()
        .find(|channel| matches!(channel.channel_type, ChannelType::Frustum { .. }))
        .expect("venturi blueprint should render a frustum channel");
    assert!(throat.width > 0.0);
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

    assert!(upper.len() >= 5);
    assert!(lower.len() >= 5);
    assert!((upper[1].0 - lower[1].0).abs() < 1e-9);
    assert!((upper[2].0 - lower[2].0).abs() < 1e-9);
    assert!(((upper[1].1 + lower[1].1) - 2.0 * y_mid).abs() < 1e-9);
    assert!(((upper[2].1 + lower[2].1) - 2.0 * y_mid).abs() < 1e-9);
    assert!(
        upper
            .windows(2)
            .any(|segment| (segment[0].1 - segment[1].1).abs() < 1e-9),
        "upper bypass lane should retain horizontal tree segments across split columns"
    );
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

    let system = channel_system_from_blueprint(&bp, Some((127.76, 85.47)))
        .expect("blueprint-native schematic conversion must succeed");
    let inlet = bp
        .nodes
        .iter()
        .find(|node| node.id.as_str() == "inlet")
        .and_then(|node| {
            node.layout.or_else(|| {
                node.metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<NodeLayoutMetadata>())
                    .copied()
            })
        })
        .expect("inlet node layout must exist");
    let outlet = bp
        .nodes
        .iter()
        .find(|node| node.id.as_str() == "outlet")
        .and_then(|node| {
            node.layout.or_else(|| {
                node.metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<NodeLayoutMetadata>())
                    .copied()
            })
        })
        .expect("outlet node layout must exist");
    assert_eq!(inlet.x_mm, 0.0);
    assert_eq!(outlet.x_mm, 127.76);
    assert!(system
        .box_outline
        .iter()
        .any(|(p0, p1)| p0.0 == 0.0 && p1.0 > 0.0));
    assert!(system
        .box_outline
        .iter()
        .any(|(p0, p1)| p0.0 < 127.76 && (p1.0 - 127.76).abs() < 1e-9));
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

    let system = channel_system_from_blueprint(&bp, Some((127.76, 85.47)))
        .expect("cct blueprint render conversion must succeed");
    let center_lane = system
        .channels
        .iter()
        .find(|channel| matches!(channel.channel_type, ChannelType::Serpentine { .. }))
        .expect("center treatment lane should stay serpentine in render path");
    match &center_lane.channel_type {
        ChannelType::Serpentine { path } => assert!(
            path.len() > 2,
            "serpentine centerline should contain explicit intermediate bends"
        ),
        _ => unreachable!("channel was filtered as serpentine"),
    }
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

    assert!(node_has_layout(&bp, "split1"));
    assert!(node_has_layout(&bp, "split2_upper"));
    assert!(node_has_layout(&bp, "split2"));
    assert!(node_has_layout(&bp, "split2_lower"));

    let upper = channel_path(&bp, "upper_C");
    let lower = channel_path(&bp, "lower_C");
    let y_mid = node_y(&bp, "split1");
    assert!(upper.len() >= 4);
    assert!(lower.len() >= 4);
    assert!((upper[1].0 - lower[1].0).abs() < 1e-9);
    assert!(((upper[1].1 + lower[1].1) - 2.0 * y_mid).abs() < 1e-9);

    let tri_nodes = bp
        .nodes
        .iter()
        .filter(|node| {
            node.junction_geometry.as_ref().or_else(|| {
                node.metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<JunctionGeometryMetadata>())
            })
                .is_some_and(|meta| meta.junction_family == JunctionFamily::Trifurcation)
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

    let system = channel_system_from_blueprint(&bp, Some((127.76, 85.47)))
        .expect("acoustic dtcv render conversion must succeed");
    assert!(node_has_layout(&bp, "split2_upper"));
    assert!(channel_path(&bp, "upper_L").len() >= 4);
    assert!(!system
        .channels
        .iter()
        .any(|channel| matches!(channel.channel_type, ChannelType::Frustum { .. })));
}
