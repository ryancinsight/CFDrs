#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; all nodes use NodeSpec::new_at() with explicit positions.

use cfd_schematics::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

fn crossing_system() -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new("crossing");
    bp.box_dims = (1.0, 1.0);

    bp.add_node(NodeSpec::new_at("0", NodeKind::Inlet, (0.0, 0.5)));
    bp.add_node(NodeSpec::new_at("1", NodeKind::Outlet, (1.0, 0.5)));
    bp.add_node(NodeSpec::new_at("2", NodeKind::Inlet, (0.5, 0.0)));
    bp.add_node(NodeSpec::new_at("3", NodeKind::Outlet, (0.5, 1.0)));

    // Explicit paths are required for crossing detection.
    let mut ch0 = ChannelSpec::new_pipe_rect("ch0", "0", "1", 1.0, 0.1, 0.05, 0.0, 0.0);
    ch0.path = vec![(0.0, 0.5), (1.0, 0.5)];
    let mut ch1 = ChannelSpec::new_pipe_rect("ch1", "2", "3", 1.0, 0.1, 0.05, 0.0, 0.0);
    ch1.path = vec![(0.5, 0.0), (0.5, 1.0)];
    bp.add_channel(ch0);
    bp.add_channel(ch1);
    bp
}

#[test]
fn overlap_resolution_survives_json_roundtrip() {
    let original = crossing_system();
    assert_eq!(original.unresolved_channel_overlap_count(), 1);
    assert!(
        original.validate().is_err(),
        "unresolved crossing systems must fail final blueprint validation"
    );
    let json = original.to_json().expect("serialize original system");
    let mut imported = NetworkBlueprint::from_json(&json).expect("deserialize original system");

    let result = imported.resolve_channel_overlaps();
    assert_eq!(result.intersection_count, 1);
    assert_eq!(result.junction_node_ids.len(), 1);
    assert_eq!(imported.nodes.len(), 5);
    assert_eq!(imported.channels.len(), 4);

    // junction_node_ids[0] is the index into the nodes Vec.
    let jn = &imported.nodes[result.junction_node_ids[0]];
    assert!((jn.point.0 - 0.5).abs() < 1e-6);
    assert!((jn.point.1 - 0.5).abs() < 1e-6);

    imported.validate().expect("resolved system must be valid");

    let json2 = imported.to_json().expect("serialize resolved system");
    let imported2 = NetworkBlueprint::from_json(&json2).expect("deserialize resolved system");
    assert_eq!(imported2.nodes.len(), 5);
    assert_eq!(imported2.channels.len(), 4);
}

#[test]
fn resolved_interchange_contains_inserted_junction_nodes() {
    let system = crossing_system();
    let unresolved = system.to_interchange();
    let resolved = system.to_interchange_resolved();

    assert_eq!(unresolved.nodes.len(), 4);
    assert_eq!(unresolved.channels.len(), 2);
    assert_eq!(resolved.nodes.len(), 5);
    assert_eq!(resolved.channels.len(), 4);
}
