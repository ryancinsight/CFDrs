//! Comprehensive Serpentine Channel Demo
//!
//! This example demonstrates all aspects of serpentine channel generation:
//! - Wave shapes (sine vs square)
//! - Phase directions (auto, inward, outward)
//! - Different configurations (smooth, high-density, optimized)
//! - Integration with split patterns
//! - Gaussian envelope improvements

use cfd_schematics::{
    config::{presets, ChannelTypeConfig, GeometryConfig, SerpentineConfig},
    geometry::{generator::create_geometry, SplitType},
    visualizations::schematic::plot_geometry,
};
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("serpentine/wave_shapes"))?;
    fs::create_dir_all(out.join("serpentine/configurations"))?;
    fs::create_dir_all(out.join("serpentine/phase_directions"))?;
    fs::create_dir_all(out.join("serpentine/optimization"))?;

    let config = GeometryConfig::default();
    let splits = vec![SplitType::Bifurcation, SplitType::Trifurcation];
    let box_dims = (300.0, 150.0);

    println!("🌊 Comprehensive Serpentine Channel Demonstration");
    println!("===============================================");

    println!("\n🎵 Wave Shape Comparison:");

    let base_config = SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 4.0,
        wave_density_factor: 2.5,
        ..SerpentineConfig::default()
    };

    let sine_config = base_config.with_sine_wave();
    let sine_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(sine_config),
    );
    let sine_output = out.join("serpentine/wave_shapes/sine_wave.png");
    plot_geometry(&sine_system, sine_output.to_str().unwrap())?;
    println!("   ✓ Sine wave: Smooth, natural curves -> {}", sine_output.display());

    let square_config = base_config.with_square_wave();
    let square_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(square_config),
    );
    let square_output = out.join("serpentine/wave_shapes/square_wave.png");
    plot_geometry(&square_system, square_output.to_str().unwrap())?;
    println!("   ✓ Square wave: Angular transitions -> {}", square_output.display());

    println!("\n🔄 Phase Direction Control:");

    let phase_configs = vec![
        ("auto_symmetric", 0.0, "Perfect bilateral mirror symmetry"),
        ("inward_phase", -1.0, "All waves start inward"),
        ("outward_phase", 1.0, "All waves start outward"),
    ];

    for (name, phase_direction, description) in phase_configs {
        let phase_config = SerpentineConfig {
            wave_phase_direction: phase_direction,
            ..base_config
        };

        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllSerpentine(phase_config),
        );
        let output = out.join(format!("serpentine/phase_directions/{}.png", name));
        plot_geometry(&system, output.to_str().unwrap())?;
        println!("   ✓ {}: {} -> {}", name, description, output.display());
    }

    println!("\n⚙️  Configuration Presets:");

    let preset_configs = vec![
        (
            "default",
            SerpentineConfig::default(),
            "Standard configuration",
        ),
        (
            "smooth",
            presets::smooth_serpentine(),
            "Smooth, low-density waves",
        ),
        (
            "high_density",
            presets::high_density_serpentine(),
            "High-density, detailed waves",
        ),
        (
            "square_wave",
            presets::square_wave_serpentine(),
            "Angular square wave preset",
        ),
    ];

    for (name, preset_config, description) in preset_configs {
        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllSerpentine(preset_config),
        );
        let output = out.join(format!("serpentine/configurations/{}.png", name));
        plot_geometry(&system, output.to_str().unwrap())?;
        println!("   ✓ {}: {} -> {}", name, description, output.display());
    }

    println!("\n🚀 Optimization Features:");

    let optimization_configs = vec![
        ("standard", SerpentineConfig::default(), "No optimization"),
        (
            "fast_optimized",
            presets::fast_optimized_serpentine(),
            "Fast optimization profile",
        ),
        (
            "thorough_optimized",
            presets::thorough_optimized_serpentine(),
            "Thorough optimization profile",
        ),
    ];

    for (name, opt_config, description) in optimization_configs {
        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllSerpentine(opt_config),
        );
        let output = out.join(format!("serpentine/optimization/{}.png", name));
        plot_geometry(&system, output.to_str().unwrap())?;

        let total_length: f64 = system
            .channels
            .iter()
            .map(|channel| match &channel.channel_type {
                cfd_schematics::geometry::ChannelType::Serpentine { path } => path
                    .windows(2)
                    .map(|w| {
                        let dx = w[1].0 - w[0].0;
                        let dy = w[1].1 - w[0].1;
                        (dx * dx + dy * dy).sqrt()
                    })
                    .sum::<f64>(),
                _ => 0.0,
            })
            .sum();

        println!(
            "   ✓ {}: {} (total length: {:.1}mm) -> {}",
            name, description, total_length, output.display()
        );
    }

    println!("\n🔗 Single Channel Examples:");

    let single_configs = vec![
        ("basic_sine", SerpentineConfig::default().with_sine_wave()),
        (
            "basic_square",
            SerpentineConfig::default().with_square_wave(),
        ),
        (
            "high_density_sine",
            presets::high_density_serpentine().with_sine_wave(),
        ),
        (
            "high_density_square",
            presets::high_density_serpentine().with_square_wave(),
        ),
    ];

    for (name, single_config) in single_configs {
        let system = create_geometry(
            (200.0, 100.0),
            &[],
            &config,
            &ChannelTypeConfig::AllSerpentine(single_config),
        );
        let output = out.join(format!("serpentine/wave_shapes/single_{}.png", name));
        plot_geometry(&system, output.to_str().unwrap())?;
        println!("   ✓ Single channel {}: -> {}", name, output.display());
    }

    println!("\n📊 Feature Summary:");
    println!("   • Wave Shapes: Sine (smooth) and Square (angular) with 200+ points for smoothness");
    println!("   • Phase Control: Auto-symmetric, inward, and outward phase directions");
    println!("   • Configurations: Default, smooth, high-density, and square wave presets");
    println!("   • Optimization: Fast and thorough optimization profiles for length maximization");
    println!("   • Envelope Functions: Improved Gaussian with distance-based normalization");
    println!(
        "   • Perfect Symmetry: Bilateral mirror symmetry maintained across all configurations"
    );

    println!("\n✅ Comprehensive serpentine demonstration complete!");
    println!("   All outputs organized in {}/serpentine/ subdirectories", out.display());

    Ok(())
}
