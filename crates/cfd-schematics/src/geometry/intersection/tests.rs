use super::detection::{channel_centerline_points, segment_intersection};
use super::*;
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::geometry::SplitType;
use crate::geometry::generator::create_geometry;
use std::borrow::Cow;
use std::collections::HashMap;

#[test]
fn segment_intersection_detects_crossing() {
    let result = segment_intersection((0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0));
    let (t, u, point) = result.expect("structural invariant");
    assert!((t - 0.5).abs() < 1e-10);
    assert!((u - 0.5).abs() < 1e-10);
    assert!((point.0 - 0.5).abs() < 1e-10);
    assert!((point.1 - 0.5).abs() < 1e-10);
}

#[test]
fn segment_intersection_detects_no_crossing_for_parallel() {
    let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0));
    assert!(result.is_none());
}

#[test]
fn segment_intersection_excludes_shared_endpoints() {
    let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (2.0, 1.0));
    assert!(result.is_none());
}

#[test]
fn centerline_storage_matches_path_shape() {
    let node_points = HashMap::from([("a".to_string(), (0.0, 0.0)), ("b".to_string(), (1.0, 0.0))]);

    let mut complete = ChannelSpec::new_pipe_rect("complete", "a", "b", 1.0, 0.1, 0.05, 0.0, 0.0);
    complete.path = vec![(0.25, 0.0), (0.75, 0.0)];
    let complete_centerline = channel_centerline_points(&complete, &node_points);
    assert!(matches!(&complete_centerline, Cow::Borrowed(_)));
    assert_eq!(complete_centerline.as_ref(), complete.path.as_slice());

    let empty = ChannelSpec::new_pipe_rect("empty", "a", "b", 1.0, 0.1, 0.05, 0.0, 0.0);
    let empty_centerline = channel_centerline_points(&empty, &node_points);
    assert!(matches!(&empty_centerline, Cow::Owned(_)));
    assert_eq!(
        empty_centerline.as_ref(),
        [(0.0, 0.0), (1.0, 0.0)].as_slice()
    );

    let mut singleton = ChannelSpec::new_pipe_rect("singleton", "a", "b", 1.0, 0.1, 0.05, 0.0, 0.0);
    singleton.path = vec![(0.5, 0.25)];
    let singleton_centerline = channel_centerline_points(&singleton, &node_points);
    assert!(matches!(&singleton_centerline, Cow::Owned(_)));
    assert_eq!(
        singleton_centerline.as_ref(),
        [(0.0, 0.0), (0.5, 0.25), (1.0, 0.0)].as_slice()
    );
}

#[test]
fn adaptive_box_dims_scales_with_branches() {
    let (w1, h1) = adaptive_box_dims(45.0, 45.0, 1, 2.0, 2.0);
    let (w2, h2) = adaptive_box_dims(45.0, 45.0, 4, 2.0, 2.0);
    let (_w8, h8) = adaptive_box_dims(45.0, 45.0, 8, 2.0, 2.0);

    assert!((w1 - 45.0).abs() < 1e-10);
    assert!((w2 - 45.0).abs() < 1e-10);
    assert!(h2 > h1, "4 branches should use more height than 1");
    assert!(h8 > h2, "8 branches should use more height than 4");
    assert!(h8 <= 45.0 + 1e-10);
}

#[test]
fn no_intersections_for_simple_bifurcation() {
    let mut system = create_geometry(
        (100.0, 50.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllStraight,
    );
    let result = insert_intersection_nodes(&mut system);
    assert_eq!(result.intersection_count, 0);
}

#[test]
fn intersection_detection_finds_manual_crossing() {
    let mut system = NetworkBlueprint::new_with_explicit_positions("crossing");
    system.box_dims = (1.0, 1.0);
    system.nodes = vec![
        NodeSpec::new_at("0".to_string(), NodeKind::Junction, (0.0, 0.5)),
        NodeSpec::new_at("1".to_string(), NodeKind::Junction, (1.0, 0.5)),
        NodeSpec::new_at("2".to_string(), NodeKind::Junction, (0.5, 0.0)),
        NodeSpec::new_at("3".to_string(), NodeKind::Junction, (0.5, 1.0)),
    ];

    let mut ch0 = ChannelSpec::new_pipe_rect(
        "ch0".to_string(),
        "0".to_string(),
        "1".to_string(),
        1.0,
        0.1,
        0.05,
        0.0,
        0.0,
    );
    ch0.path = vec![(0.0, 0.5), (1.0, 0.5)];
    let mut ch1 = ChannelSpec::new_pipe_rect(
        "ch1".to_string(),
        "2".to_string(),
        "3".to_string(),
        1.0,
        0.1,
        0.05,
        0.0,
        0.0,
    );
    ch1.path = vec![(0.5, 0.0), (0.5, 1.0)];
    system.channels = vec![ch0, ch1];

    let result = insert_intersection_nodes(&mut system);
    assert_eq!(result.intersection_count, 1);
    assert_eq!(result.junction_node_ids.len(), 1);
    assert_eq!(
        system.channels.len(),
        4,
        "2 channels × 1 crossing = 4 sub-channels"
    );

    let junction = &system.nodes[result.junction_node_ids[0]];
    assert!((junction.point.0 - 0.5).abs() < 1e-6);
    assert!((junction.point.1 - 0.5).abs() < 1e-6);
}
