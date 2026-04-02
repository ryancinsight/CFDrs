use cfd_schematics::{
    config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig},
    geometry::{
        builders::ChannelExt,
        generator::{create_geometry, create_geometry_with_metadata, MetadataConfig},
        metadata::PerformanceMetadata,
        SplitType,
    },
};

#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    tracing::info!("Unified Generator Demo");
    tracing::info!("====================");

    let config = GeometryConfig::default();
    let serpentine_config = SerpentineConfig::default();

    tracing::info!("1. Standard Generation (No Metadata)");
    let start_time = std::time::Instant::now();
    let standard_system = create_geometry(
        (200.0, 100.0),
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllSerpentine(serpentine_config),
    );
    let standard_time = start_time.elapsed();

    tracing::info!(
        "   Generated {} channels and {} nodes in {:?}",
        standard_system.channels.len(),
        standard_system.nodes.len(),
        standard_time
    );

    let has_metadata = standard_system.channels[0].has_metadata::<PerformanceMetadata>();
    tracing::info!("   Has performance metadata: {}", has_metadata);

    tracing::info!("2. Generation with Metadata");
    let metadata_config = MetadataConfig {
        track_performance: true,
        track_optimization: true,
        channel_diameter_mm: None,
    };

    let start_time = std::time::Instant::now();
    let metadata_system = create_geometry_with_metadata(
        (200.0, 100.0),
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllSerpentine(serpentine_config),
        &metadata_config,
    );
    let metadata_time = start_time.elapsed();

    tracing::info!(
        "   Generated {} channels and {} nodes in {:?}",
        metadata_system.channels.len(),
        metadata_system.nodes.len(),
        metadata_time
    );

    let has_perf_metadata = metadata_system.channels[0].has_metadata::<PerformanceMetadata>();
    tracing::info!("   Has performance metadata: {}", has_perf_metadata);

    if has_perf_metadata {
        match metadata_system.channels[0].get_metadata::<PerformanceMetadata>() {
            Some(perf_data) => {
                tracing::info!(
                    "   Generation time: {}μs, Memory: {} bytes, Path points: {}",
                    perf_data.generation_time_us,
                    perf_data.memory_usage_bytes,
                    perf_data.path_points_count
                );
            }
            None => {
                tracing::info!(
                    "   Warning: Performance metadata exists but could not be retrieved"
                );
            }
        }
    }

    tracing::info!("3. Performance Comparison");
    tracing::info!("   Standard generation: {:?}", standard_time);
    tracing::info!("   Metadata generation: {:?}", metadata_time);
    let overhead = metadata_time.as_nanos() as f64 / standard_time.as_nanos() as f64;
    tracing::info!("   Metadata overhead: {:.1}x", overhead);

    tracing::info!("4. Generating Visualizations");

    shared::output::save_example_output_with_name(
        &standard_system,
        "unified_generator_demo",
        "standard_system",
    );
    shared::output::save_example_output_with_name(
        &metadata_system,
        "unified_generator_demo",
        "metadata_system",
    );

    tracing::info!("5. Unified Generator API Summary");
    tracing::info!("   • create_geometry() - Fast generation without metadata");
    tracing::info!(
        "   • create_geometry_with_metadata() - Generation with optional metadata tracking"
    );
    tracing::info!("   • MetadataConfig - Configure what metadata to track");
    tracing::info!("   • Zero overhead when metadata is disabled");
    tracing::info!("   • Canonical API preserved across generator variants");

    tracing::info!("Demo complete! The unified generator provides:");
    tracing::info!("• Single source of truth for geometry generation");
    tracing::info!("• Optional metadata support with zero overhead when unused");
    tracing::info!("• Clean, consistent API");
    tracing::info!("• Canonical API consistency");

    Ok(())
}
