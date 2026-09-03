use super::*;
use crate::config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig};
use crate::domain::model::ChannelShape;
use crate::geometry::builders::{ChannelExt, NodeExt};
use crate::geometry::metadata::{ChannelGeometryMetadata, PerformanceMetadata};
use crate::geometry::types::{polyline_length, ChannelType};
use crate::geometry::SplitType;
use crate::topology::presets::single_venturi_series_spec;
use aequitas::systems::si::quantities::Length;

#[test]
fn test_generator_with_performance_metadata() {
    let metadata_config = MetadataConfig {
        track_performance: true,
        track_optimization: false,
        channel_diameter_m: None,
    };

    let system = create_geometry_with_metadata(
        (100.0, 50.0),
        &[],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllStraight,
        &metadata_config,
    );

    for channel in &system.channels {
        assert!(channel.has_metadata::<PerformanceMetadata>());
        let perf_data = channel
            .get_metadata::<PerformanceMetadata>()
            .expect("Performance metadata should be available after creation");
        assert!(perf_data.generation_time_us > 0);
        assert!(perf_data.memory_usage_bytes > 0);
    }

    for node in &system.nodes {
        assert!(node.has_metadata::<PerformanceMetadata>());
    }
}

#[test]
fn test_channel_diameter_metadata_adjusts_split_spacing() {
    let box_dims = (200.0, 100.0);
    let splits = [SplitType::Bifurcation];
    let config = GeometryConfig::default();
    let channel_type_config = ChannelTypeConfig::AllStraight;

    let baseline_metadata = MetadataConfig::default();
    let baseline_system = create_geometry_with_metadata(
        box_dims,
        &splits,
        &config,
        &channel_type_config,
        &baseline_metadata,
    );

    let diameter_metadata =
        MetadataConfig::default().with_channel_diameter_m(Length::from_base(40.0e-3));
    let diameter_system = create_geometry_with_metadata(
        box_dims,
        &splits,
        &config,
        &channel_type_config,
        &diameter_metadata,
    );

    let baseline_min_y = baseline_system
        .nodes
        .iter()
        .map(|node| node.point.1)
        .fold(f64::INFINITY, f64::min);
    let diameter_min_y = diameter_system
        .nodes
        .iter()
        .map(|node| node.point.1)
        .fold(f64::INFINITY, f64::min);

    assert!(
        diameter_min_y > baseline_min_y + 1.0,
        "Larger diameter metadata should push branches farther from walls (baseline_min_y={baseline_min_y:.3}, diameter_min_y={diameter_min_y:.3})"
    );

    for channel in &diameter_system.channels {
        let metadata = channel
            .get_metadata::<ChannelGeometryMetadata>()
            .expect("channel should include ChannelGeometryMetadata when diameter is configured");
        assert!((metadata.channel_diameter_m.into_base() - 40.0e-3).abs() < 1e-15);
    }
}

#[test]
fn generated_serpentine_channels_persist_physical_length_and_shape() {
    let system = create_geometry(
        (200.0, 100.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
    );

    let serpentine_channels: Vec<_> = system
        .channels
        .iter()
        .filter(|channel| channel.path.len() > 2)
        .collect();
    assert!(
        !serpentine_channels.is_empty(),
        "expected generated geometry to include routed serpentine channels"
    );

    for channel in serpentine_channels {
        let expected_length_m = polyline_length(&channel.path) * 1.0e-3;
        assert!(channel.length_m.into_base() > 0.0);
        assert!(
            (channel.length_m.into_base() - expected_length_m).abs() < 1e-9,
            "channel {:?} length should follow stored polyline length",
            channel.id
        );

        match channel.channel_shape {
            ChannelShape::Serpentine {
                segments,
                bend_radius_m,
                wave_type: _,
            } => {
                assert!(
                    segments >= 2,
                    "channel {:?} should expose at least one serpentine bend",
                    channel.id
                );
                assert!(
                    bend_radius_m.into_base() > 0.0,
                    "channel {:?} should expose a positive bend radius",
                    channel.id
                );
            }
            shape @ ChannelShape::Straight => panic!(
                "channel {:?} should be marked serpentine for 1D modeling, got {:?}",
                channel.id, shape
            ),
        }
    }
}

#[test]
fn explicit_channel_paths_survive_owned_construction_modes() {
    let path = vec![(0.0, 0.0), (20.0, 3.0), (40.0, 0.0)];

    let mut generator = GeometryGenerator::new(
        (100.0, 50.0),
        GeometryConfig::default(),
        ChannelTypeConfig::AllStraight,
        1,
    );
    generator.add_channel_with_type(
        (0.0, 0.0),
        (40.0, 0.0),
        Some(ChannelType::Serpentine { path: path.clone() }),
        None,
    );

    let mut metadata_generator = GeometryGenerator::new_with_metadata(
        (100.0, 50.0),
        GeometryConfig::default(),
        ChannelTypeConfig::AllStraight,
        1,
        MetadataConfig::default(),
    );
    metadata_generator.add_channel_with_type(
        (0.0, 0.0),
        (40.0, 0.0),
        Some(ChannelType::Serpentine { path: path.clone() }),
        None,
    );

    assert_eq!(generator.channels[0].path, path);
    assert_eq!(metadata_generator.channels[0].path, path);
}

#[test]
fn series_geometry_uses_adjacent_node_coordinates() {
    let spec = single_venturi_series_spec("series-index", 1.0e-3, 100.0e-6, 500.0e-6, 200.0e-6);
    let blueprint = create_series_geometry_from_spec(&spec);

    assert_eq!(blueprint.channels.len(), 3);
    assert_eq!(blueprint.nodes.len(), blueprint.channels.len() + 1);
    assert_eq!(
        blueprint
            .channels
            .iter()
            .map(|channel| (channel.from.0.as_str(), channel.to.0.as_str()))
            .collect::<Vec<_>>(),
        vec![
            ("inlet", "junction_0"),
            ("junction_0", "junction_1"),
            ("junction_1", "outlet"),
        ]
    );

    for (idx, channel) in blueprint.channels.iter().enumerate() {
        assert_eq!(
            channel.path.first().copied(),
            Some(blueprint.nodes[idx].point)
        );
        assert_eq!(
            channel.path.last().copied(),
            Some(blueprint.nodes[idx + 1].point)
        );
    }
}
