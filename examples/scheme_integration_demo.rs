//! Example demonstrating scheme library integration for 2D microfluidic schematics

#[cfg(feature = "scheme-integration")]
use cfd_1d::prelude::*;

#[cfg(not(feature = "scheme-integration"))]
fn main() {
    println!("Scheme Integration Demo");
    println!("======================");
    println!();
    println!("❌ The 'scheme-integration' feature is not enabled.");
    println!();
    // TODO: This feature requires system dependencies (fontconfig) that may not be available
    // DEPENDENCIES: fontconfig library, pkg-config tool
    // BLOCKED BY: Missing system dependencies on some platforms
    // PRIORITY: Medium - Nice-to-have visualization feature
    println!("in all environments. To enable scheme integration:");
    println!();
    println!("1. Install system dependencies:");
    println!("   - Ubuntu/Debian: sudo apt install pkg-config libfontconfig1-dev");
    println!("   - CentOS/RHEL: sudo yum install pkgconfig fontconfig-devel");
    println!("   - macOS: brew install pkg-config");
    println!();
    println!("2. Run with the feature enabled:");
    println!("   cargo run --example scheme_integration_demo --features cfd-1d/scheme-integration");
    println!();
    println!("Alternatively, you can design microfluidic networks programmatically");
    println!("using the NetworkBuilder API without 2D schematic visualization.");
}

#[cfg(feature = "scheme-integration")]
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    use scheme::{
        config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig},
        geometry::{generator::create_geometry, SplitType},
        visualizations::schematic::plot_geometry,
    };

    println!("Scheme Integration Demo");
    println!("======================");

    // Create a simple bifurcation network
    println!("\n1. Creating a simple bifurcation network...");
    let bifurcation_system = helpers::create_bifurcation_schematic(200.0, 100.0, 2)?;
    println!(
        "   Created system with {} nodes and {} channels",
        bifurcation_system.nodes.len(),
        bifurcation_system.channels.len()
    );

    // Export to PNG
    println!("\n2. Exporting to PNG...");
    helpers::export_to_png(&bifurcation_system, "bifurcation.png")?;
    println!("   Saved to bifurcation.png");

    // Create a more complex system with mixed channel types
    println!("\n3. Creating a complex system with serpentine channels...");
    let config = GeometryConfig::default();
    let serpentine_config = SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        ..SerpentineConfig::default()
    };

    let complex_system = create_geometry(
        (300.0, 150.0),
        &[SplitType::Bifurcation, SplitType::Trifurcation],
        &config,
        &ChannelTypeConfig::AllSerpentine(serpentine_config),
    );

    println!(
        "   Created system with {} nodes and {} channels",
        complex_system.nodes.len(),
        complex_system.channels.len()
    );

    // Extract layout information
    println!("\n4. Extracting layout information...");
    let nodes = helpers::convert_nodes(&complex_system.nodes);
    let paths = helpers::extract_channel_paths(&complex_system.channels);

    println!("   Node positions:");
    for (id, pos) in nodes.iter().take(5) {
        println!("     Node {}: ({:.2}, {:.2})", id, pos.0, pos.1);
    }
    if nodes.len() > 5 {
        println!("     ... and {} more nodes", nodes.len() - 5);
    }

    println!("\n   Channel paths:");
    for (i, path) in paths.iter().enumerate().take(3) {
        println!(
            "     Channel {}: {} -> {} ({:?}, {} points)",
            i,
            path.source,
            path.target,
            path.channel_type,
            path.points.len()
        );
    }
    if paths.len() > 3 {
        println!("     ... and {} more channels", paths.len() - 3);
    }

    // Export complex system
    println!("\n5. Exporting complex system...");
    plot_geometry(&complex_system, "complex_serpentine.png")?;
    println!("   Saved to complex_serpentine.png");

    // Demonstrate JSON serialization
    println!("\n6. Serializing to JSON...");
    let json = serde_json::to_string_pretty(&complex_system)?;
    std::fs::write("complex_system.json", &json)?;
    println!("   Saved to complex_system.json ({} bytes)", json.len());

    println!("\nScheme integration demo completed successfully!");

    Ok(())
}
