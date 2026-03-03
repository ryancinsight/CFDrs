use cfd_schematics::geometry::{Channel, ChannelSystem, ChannelType, Node};

fn crossing_system() -> ChannelSystem {
    ChannelSystem {
        box_dims: (1.0, 1.0),
        nodes: vec![
            Node {
                id: 0,
                point: (0.0, 0.5),
                metadata: None,
            },
            Node {
                id: 1,
                point: (1.0, 0.5),
                metadata: None,
            },
            Node {
                id: 2,
                point: (0.5, 0.0),
                metadata: None,
            },
            Node {
                id: 3,
                point: (0.5, 1.0),
                metadata: None,
            },
        ],
        channels: vec![
            Channel {
                id: 0,
                from_node: 0,
                to_node: 1,
                width: 0.1,
                height: 0.05,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 1,
                from_node: 2,
                to_node: 3,
                width: 0.1,
                height: 0.05,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
        ],
        box_outline: vec![],
    }
}

#[test]
fn overlap_resolution_survives_json_roundtrip() {
    let original = crossing_system();
    let json = original.to_json().expect("serialize original system");
    let mut imported = ChannelSystem::from_json(&json).expect("deserialize original system");

    let result = imported.resolve_channel_overlaps();
    assert_eq!(result.intersection_count, 1);
    assert_eq!(result.junction_node_ids.len(), 1);
    assert_eq!(imported.nodes.len(), 5);
    assert_eq!(imported.channels.len(), 4);

    let jn = &imported.nodes[result.junction_node_ids[0]];
    assert!((jn.point.0 - 0.5).abs() < 1e-6);
    assert!((jn.point.1 - 0.5).abs() < 1e-6);

    imported.validate().expect("resolved system must be valid");

    let json2 = imported.to_json().expect("serialize resolved system");
    let imported2 = ChannelSystem::from_json(&json2).expect("deserialize resolved system");
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
