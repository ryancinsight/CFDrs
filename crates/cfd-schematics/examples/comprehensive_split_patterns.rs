//! Comprehensive Split Patterns Demo
//!
//! This example demonstrates various split patterns (bifurcation and trifurcation)
//! with different levels of complexity. It replaces multiple individual examples
//! with a single comprehensive demonstration.

use cfd_schematics::{
    config::{ChannelTypeConfig, GeometryConfig},
    geometry::{generator::create_geometry, SplitType},
    visualizations::schematic::plot_geometry,
};
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("split_patterns/bifurcation"))?;
    fs::create_dir_all(out.join("split_patterns/trifurcation"))?;
    fs::create_dir_all(out.join("split_patterns/mixed"))?;

    let config = GeometryConfig::default();
    let channel_config = ChannelTypeConfig::default();

    println!("🔬 Generating Comprehensive Split Pattern Examples");
    println!("================================================");

    println!("\n📊 Bifurcation Patterns:");

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

    for (name, splits, box_dims) in bifurcation_patterns {
        let system = create_geometry(box_dims, &splits, &config, &channel_config);
        let output_path = out.join(format!(
            "split_patterns/bifurcation/{}_bifurcation.png",
            name
        ));
        plot_geometry(&system, output_path.to_str().unwrap())?;

        println!(
            "   ✓ {}: {} channels, {} nodes -> {}",
            name,
            system.channels.len(),
            system.nodes.len(),
            output_path.display()
        );
    }

    println!("\n🔱 Trifurcation Patterns:");

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
        let output_path = out.join(format!(
            "split_patterns/trifurcation/{}_trifurcation.png",
            name
        ));
        plot_geometry(&system, output_path.to_str().unwrap())?;

        println!(
            "   ✓ {}: {} channels, {} nodes -> {}",
            name,
            system.channels.len(),
            system.nodes.len(),
            output_path.display()
        );
    }

    println!("\n🔀 Mixed Patterns:");

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
        let output_path = out.join(format!("split_patterns/mixed/{}_pattern.png", name));
        plot_geometry(&system, output_path.to_str().unwrap())?;

        println!(
            "   ✓ {}: {} channels, {} nodes -> {}",
            name,
            system.channels.len(),
            system.nodes.len(),
            output_path.display()
        );
    }

    println!("\n📈 Summary:");
    println!("   • Bifurcation patterns: Each split creates 2 branches");
    println!("   • Trifurcation patterns: Each split creates 3 branches");
    println!("   • Mixed patterns: Combinations of both split types");
    println!("   • Smart channel selection: Automatically chooses optimal channel types");
    println!("   • All outputs saved to organized directory structure");

    println!("\n✅ Comprehensive split pattern generation complete!");

    Ok(())
}
