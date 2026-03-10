//! Comprehensive Arc Channel Demo
//!
//! This example demonstrates all aspects of arc channel generation:
//! - Different curvature factors
//! - Smoothness levels
//! - Directional control
//! - Integration with split patterns
//! - Smart mixed channel configurations

use cfd_schematics::{
    config::{presets, ArcConfig, ChannelTypeConfig, GeometryConfig},
    geometry::{generator::create_geometry, SplitType},
    visualizations::schematic::plot_geometry,
};
#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = GeometryConfig::default();
    let splits = vec![SplitType::Bifurcation, SplitType::Trifurcation];
    let box_dims = (300.0, 150.0);

    tracing::info!("🌈 Comprehensive Arc Channel Demonstration");
    tracing::info!("========================================");

    tracing::info!("\n📐 Curvature Factor Variations:");

    let curvature_factors = vec![
        (0.1, "subtle", "Very subtle curves"),
        (0.5, "moderate", "Moderate curvature"),
        (1.0, "pronounced", "Pronounced curves"),
        (1.5, "strong", "Strong curvature"),
        (2.0, "maximum", "Maximum curvature"),
    ];

    for (factor, name, description) in curvature_factors {
        let arc_config = ArcConfig {
            curvature_factor: factor,
            smoothness: 50,
            curvature_direction: 0.0,
            min_separation_distance: 1.0,
            enable_collision_prevention: true,
            max_curvature_reduction: 0.5,
            enable_adaptive_curvature: true,
        };

        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllArcs(arc_config),
        );
        let filename = format!("curvature_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ {}: {} (factor: {})", name, description, factor);
    }

    tracing::info!("\n✨ Smoothness Levels:");

    let smoothness_levels = vec![
        (10, "low", "Low resolution (10 points)"),
        (25, "medium", "Medium resolution (25 points)"),
        (50, "high", "High resolution (50 points)"),
        (100, "ultra", "Ultra-high resolution (100 points)"),
    ];

    for (smoothness, name, description) in smoothness_levels {
        let arc_config = ArcConfig {
            curvature_factor: 1.0,
            smoothness,
            curvature_direction: 0.0,
            min_separation_distance: 1.0,
            enable_collision_prevention: true,
            max_curvature_reduction: 0.5,
            enable_adaptive_curvature: true,
        };

        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllArcs(arc_config),
        );
        let filename = format!("smoothness_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n🧭 Directional Control:");

    let direction_configs = vec![
        ("auto", 0.0, "Auto-determine direction"),
        ("inward", -1.0, "Force inward curvature"),
        ("outward", 1.0, "Force outward curvature"),
    ];

    for (name, direction, description) in direction_configs {
        let arc_config = ArcConfig {
            curvature_factor: 1.0,
            smoothness: 50,
            curvature_direction: direction,
            min_separation_distance: 1.0,
            enable_collision_prevention: true,
            max_curvature_reduction: 0.5,
            enable_adaptive_curvature: true,
        };

        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllArcs(arc_config),
        );
        let filename = format!("direction_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n⚙️  Configuration Presets:");

    let preset_configs = vec![
        (
            "default",
            ArcConfig::default(),
            "Standard arc configuration",
        ),
        ("subtle", presets::subtle_arcs(), "Subtle, gentle curves"),
        (
            "pronounced",
            presets::pronounced_arcs(),
            "Pronounced, dramatic curves",
        ),
    ];

    for (name, preset_config, description) in preset_configs {
        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllArcs(preset_config),
        );
        let filename = format!("preset_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n🧠 Smart Mixed Channel Selection:");

    let adaptive_configs = vec![
        (
            "adaptive_default",
            ChannelTypeConfig::Adaptive {
                serpentine_config: cfd_schematics::config::SerpentineConfig::default(),
                arc_config: ArcConfig::default(),
                frustum_config: cfd_schematics::config::FrustumConfig::default(),
            },
            "Default adaptive selection",
        ),
        (
            "adaptive_pronounced",
            ChannelTypeConfig::Adaptive {
                serpentine_config: presets::smooth_serpentine(),
                arc_config: presets::pronounced_arcs(),
                frustum_config: cfd_schematics::config::FrustumConfig::default(),
            },
            "Pronounced curves with smooth serpentines",
        ),
        (
            "adaptive_mixed_position",
            ChannelTypeConfig::MixedByPosition {
                middle_zone_fraction: 0.4,
                serpentine_config: cfd_schematics::config::SerpentineConfig::default(),
                arc_config: ArcConfig::default(),
            },
            "Position-based channel selection",
        ),
    ];

    for (name, adaptive_config, description) in adaptive_configs {
        let system = create_geometry(box_dims, &splits, &config, &adaptive_config);
        let filename = format!("adaptive_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n🪞 Enhanced Bilateral Mirror Symmetry:");

    let symmetry_configs = vec![
        (
            "trifurcation_enhanced",
            vec![SplitType::Trifurcation],
            (300.0, 150.0),
            "Trifurcation with figure-8 style center weave and bilateral symmetry",
        ),
        (
            "double_trifurcation_enhanced",
            vec![SplitType::Trifurcation, SplitType::Trifurcation],
            (400.0, 200.0),
            "Double trifurcation with enhanced symmetry",
        ),
        (
            "mixed_bifurcation_trifurcation",
            vec![SplitType::Bifurcation, SplitType::Trifurcation],
            (350.0, 175.0),
            "Mixed pattern with enhanced symmetry",
        ),
    ];

    for (name, splits, box_dims, description) in symmetry_configs {
        let enhanced_config = ArcConfig {
            curvature_factor: 1.2,
            smoothness: 75,
            curvature_direction: 0.0,
            min_separation_distance: 1.0,
            enable_collision_prevention: true,
            max_curvature_reduction: 0.5,
            enable_adaptive_curvature: true,
        };

        let system = create_geometry(
            box_dims,
            &splits,
            &config,
            &ChannelTypeConfig::AllArcs(enhanced_config),
        );
        let filename = format!("symmetry_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ {}: {}", name, description);
    }

    tracing::info!("\n🔗 Single Arc Examples:");

    let single_configs = vec![
        ("basic", ArcConfig::default()),
        ("subtle", presets::subtle_arcs()),
        ("pronounced", presets::pronounced_arcs()),
        (
            "high_smoothness",
            ArcConfig {
                curvature_factor: 1.0,
                smoothness: 100,
                curvature_direction: 0.0,
                min_separation_distance: 1.0,
                enable_collision_prevention: true,
                max_curvature_reduction: 0.5,
                enable_adaptive_curvature: true,
            },
        ),
    ];

    for (name, single_config) in single_configs {
        let system = create_geometry(
            (200.0, 100.0),
            &[],
            &config,
            &ChannelTypeConfig::AllArcs(single_config),
        );
        let filename = format!("single_{}", name);
        shared::save_example_output_with_name(&system, "comprehensive_arc_demo", &filename);
        tracing::info!("   ✓ Single arc {}", name);
    }

    tracing::info!("\n📊 Feature Summary:");
    tracing::info!("   • Curvature Control: From subtle (0.1) to maximum (2.0) curvature factors");
    tracing::info!("   • Smoothness Options: 10-100+ points for different resolution needs");
    tracing::info!("   • Directional Control: Auto, inward, and outward curvature directions");
    tracing::info!("   • Enhanced Bilateral Mirror Symmetry: Perfect symmetry across vertical and horizontal centerlines");
    tracing::info!(
        "   • Straight Inlet/Outlet Channels: Inlet and outlet channels are automatically straight"
    );
    tracing::info!("   • Figure-8 Center Channels: Center channels in trifurcations use figure-8 style weave curvature without unresolved planar crossings");
    tracing::info!("   • Peripheral vs Internal Curvature: Peripheral arcs curve toward walls, internal arcs toward center");
    tracing::info!("   • Smart Selection: Automatic channel type selection based on geometry");
    tracing::info!("   • Mixed Configurations: Combine arcs and serpentines intelligently");
    tracing::info!("   • Performance Optimized: Efficient generation for all smoothness levels");

    tracing::info!("\n✅ Comprehensive arc demonstration complete!");
    tracing::info!("   All outputs organized in output/examples/comprehensive_arc_demo/");

    Ok(())
}
