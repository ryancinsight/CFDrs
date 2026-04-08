use super::*;
use crate::config::SerpentineConfig;
use crate::domain::model::ChannelShape;
use crate::geometry::builders::{ChannelExt, NodeExt};
use crate::geometry::metadata::{ChannelGeometryMetadata, PerformanceMetadata};
use crate::geometry::SplitType;

#[test]
fn test_generator_with_performance_metadata() {
    let metadata_config = MetadataConfig {
        track_performance: true,
        track_optimization: false,
        channel_diameter_mm: None,
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

    let diameter_metadata = MetadataConfig::default().with_channel_diameter_mm(40.0);
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
        let metadata = channel.get_metadata::<ChannelGeometryMetadata>().expect(
            "channel should include ChannelGeometryMetadata when diameter is configured",
        );
        assert!((metadata.channel_diameter_mm - 40.0).abs() < 1e-10);
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
        assert!(channel.length_m > 0.0);
        assert!(
            (channel.length_m - expected_length_m).abs() < 1e-9,
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
                    bend_radius_m > 0.0,
                    "channel {:?} should expose a positive bend radius",
                    channel.id
                );
            }
            ref shape => panic!(
                "channel {:?} should be marked serpentine for 1D modeling, got {:?}",
                channel.id, shape
            ),
        }
    }
}