//! Comprehensive Split Patterns Demo
//!
//! This example demonstrates various split patterns (bifurcation and trifurcation)
//! with different levels of complexity. It replaces multiple individual examples
//! with a single comprehensive demonstration.

use cfd_schematics::{
    config::{ChannelTypeConfig, GeometryConfig},
    geometry::{generator::create_geometry, SplitType},
    visualizations::{schematic::plot_geometry_auto_annotated, traits::RenderConfig},
};
use std::path::PathBuf;

#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let config = GeometryConfig::default();
    let channel_config = ChannelTypeConfig::default();

    tracing::info!("🔬 Generating Comprehensive Split Pattern Examples");
    tracing::info!("================================================");

    tracing::info!("\n📊 Bifurcation Patterns:");

    let bifurcation_patterns = vec![
        ("single", vec![SplitType::Bifurcation], (200.0, 100.0)),
        (
            "double",
            vec![SplitType::Bifurcation, SplitType::Bifurcation],
            (300.0, 150.0),
        ),
        (
            "triple",
            vec![
                SplitType::Bifurcation,
                SplitType::Bifurcation,
                SplitType::Bifurcation,
            ],
            (400.0, 200.0),
        ),
        ("quadruple", vec![SplitType::Bifurcation; 4], (500.0, 250.0)),
    ];

    let render_config = RenderConfig::default();

    for (name, splits, box_dims) in bifurcation_patterns {
        let system = create_geometry(box_dims, &splits, &config, &channel_config);
        let filename = format!("bifurcation_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_split_patterns", &filename);

        tracing::info!(
            "   ✓ {}: {} channels, {} nodes",
            name,
            system.channels.len(),
            system.nodes.len()
        );
    }

    tracing::info!("\n🔱 Trifurcation Patterns:");

    let trifurcation_patterns = vec![
        ("single", vec![SplitType::Trifurcation], (250.0, 120.0)),
        (
            "double",
            vec![SplitType::Trifurcation, SplitType::Trifurcation],
            (400.0, 200.0),
        ),
        (
            "triple",
            vec![
                SplitType::Trifurcation,
                SplitType::Trifurcation,
                SplitType::Trifurcation,
            ],
            (600.0, 300.0),
        ),
    ];

    for (name, splits, box_dims) in trifurcation_patterns {
        let system = create_geometry(box_dims, &splits, &config, &channel_config);
        let filename = format!("trifurcation_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_split_patterns", &filename);

        tracing::info!(
            "   ✓ {}: {} channels, {} nodes",
            name,
            system.channels.len(),
            system.nodes.len()
        );
    }

    tracing::info!("\n🔀 Mixed Patterns:");

    let mixed_patterns = vec![
        (
            "bifurcation_trifurcation",
            vec![SplitType::Bifurcation, SplitType::Trifurcation],
            (300.0, 150.0),
        ),
        (
            "trifurcation_bifurcation",
            vec![SplitType::Trifurcation, SplitType::Bifurcation],
            (300.0, 150.0),
        ),
        (
            "alternating",
            vec![
                SplitType::Bifurcation,
                SplitType::Trifurcation,
                SplitType::Bifurcation,
            ],
            (400.0, 200.0),
        ),
        (
            "complex",
            vec![
                SplitType::Trifurcation,
                SplitType::Bifurcation,
                SplitType::Trifurcation,
                SplitType::Bifurcation,
            ],
            (500.0, 250.0),
        ),
    ];

    for (name, splits, box_dims) in mixed_patterns {
        let system = create_geometry(box_dims, &splits, &config, &channel_config);
        let filename = format!("mixed_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_split_patterns", &filename);

        tracing::info!(
            "   ✓ {}: {} channels, {} nodes",
            name,
            system.channels.len(),
            system.nodes.len()
        );
    }

    tracing::info!("\n📈 Summary:");
    tracing::info!("   • Bifurcation patterns: Each split creates 2 branches");
    tracing::info!("   • Trifurcation patterns: Each split creates 3 branches");
    tracing::info!("   • Mixed patterns: Combinations of both split types");
    tracing::info!("   • Smart channel selection: Automatically chooses optimal channel types");
    tracing::info!("   • All outputs saved to organized directory structure");

    tracing::info!("\n✅ Comprehensive split pattern generation complete!");

    Ok(())
}
