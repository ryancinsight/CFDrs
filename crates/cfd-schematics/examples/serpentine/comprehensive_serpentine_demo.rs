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
#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let config = GeometryConfig::default();
    let splits = vec![SplitType::Bifurcation, SplitType::Trifurcation];
    let box_dims = (300.0, 150.0);

    tracing::info!("🌊 Comprehensive Serpentine Channel Demonstration");
    tracing::info!("===============================================");

    tracing::info!("\n🎵 Wave Shape Comparison:");

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
    shared::save_example_output_with_name(&sine_system, "comprehensive_serpentine_demo", "sine_wave");
    tracing::info!(
        "   ✓ Sine wave: Smooth, natural curves"
    );

    let square_config = base_config.with_square_wave();
    let square_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(square_config),
    );
    shared::save_example_output_with_name(&square_system, "comprehensive_serpentine_demo", "square_wave");
    tracing::info!(
        "   ✓ Square wave: Angular transitions"
    );

    tracing::info!("\n🔄 Phase Direction Control:");

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
        let filename = format!("phase_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_serpentine_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n⚙️  Configuration Presets:");

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
        let filename = format!("config_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_serpentine_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n🚀 Optimization Features:");

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
        let filename = format!("opt_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_serpentine_demo", &filename);

        let total_length: f64 = system
            .channels
            .iter()
            .map(|channel| match &channel.channel_shape {
                cfd_schematics::domain::model::ChannelShape::Serpentine { .. } => channel
                    .path
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

        tracing::info!(
            "   ✓ {}: {} (total length: {:.1}mm)",
            name,
            description,
            total_length
        );
    }

    tracing::info!("\n🔗 Single Channel Examples:");

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
        let filename = format!("single_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_serpentine_demo", &filename);
        tracing::info!("   ✓ Single channel {}", name);
    }

    tracing::info!("\n📊 Feature Summary:");
    tracing::info!(
        "   • Wave Shapes: Sine (smooth) and Square (angular) with 200+ points for smoothness"
    );
    tracing::info!("   • Phase Control: Auto-symmetric, inward, and outward phase directions");
    tracing::info!("   • Configurations: Default, smooth, high-density, and square wave presets");
    tracing::info!(
        "   • Optimization: Fast and thorough optimization profiles for length maximization"
    );
    tracing::info!("   • Envelope Functions: Improved Gaussian with distance-based normalization");
    tracing::info!(
        "   • Perfect Symmetry: Bilateral mirror symmetry maintained across all configurations"
    );

    tracing::info!("\n✅ Comprehensive serpentine demonstration complete!");
    tracing::info!(
        "   All outputs organized in output/examples/comprehensive_serpentine_demo/"
    );

    Ok(())
}
