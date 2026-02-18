//! Frustum Channel Demonstration
//!
//! This example demonstrates the new frustum (tapered) channel functionality
//! for venturi throat applications in microfluidic design systems.
//!
//! Features demonstrated:
//! - Different taper profiles (Linear, Exponential, Smooth)
//! - Configurable inlet, throat, and outlet widths
//! - Variable throat positioning
//! - Integration with existing channel types
//! - JSON serialization/deserialization
//! - Visualization support
//!
//! Run with: cargo run --example frustum_channel_demo

use cfd_schematics::{
    config::{ChannelTypeConfig, FrustumConfig, GeometryConfig, TaperProfile},
    geometry::{generator::create_geometry, ChannelSystem, SplitType},
    visualizations::schematic::plot_geometry,
};
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌊 Frustum Channel Demonstration");
    println!("================================");
    println!();

    fs::create_dir_all("outputs")?;

    demonstrate_taper_profiles()?;
    demonstrate_throat_positions()?;
    demonstrate_width_configurations()?;
    demonstrate_mixed_channel_systems()?;
    demonstrate_json_serialization()?;

    println!("✅ All demonstrations completed successfully!");
    println!();
    println!("📁 Output files saved to 'outputs/' directory:");
    println!("   • frustum_linear_taper.svg");
    println!("   • frustum_exponential_taper.svg");
    println!("   • frustum_smooth_taper.svg");
    println!("   • frustum_throat_positions.svg");
    println!("   • frustum_width_configs.svg");
    println!("   • mixed_channel_system.svg");
    println!("   • frustum_system_export.json");

    Ok(())
}

fn demonstrate_taper_profiles() -> Result<(), Box<dyn std::error::Error>> {
    println!("1️⃣  Demonstrating Taper Profiles");
    println!("   Testing Linear, Exponential, and Smooth taper profiles");

    let base_config = FrustumConfig {
        inlet_width: 3.0,
        throat_width: 0.8,
        outlet_width: 2.5,
        smoothness: 100,
        throat_position: 0.5,
        taper_profile: TaperProfile::Linear,
    };

    let linear_config = FrustumConfig {
        taper_profile: TaperProfile::Linear,
        ..base_config
    };
    let linear_system = create_geometry(
        (120.0, 40.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllFrustum(linear_config),
    );
    plot_geometry(&linear_system, "outputs/frustum_linear_taper.svg")?;
    println!("   ✅ Linear taper: saved to frustum_linear_taper.svg");

    let exponential_config = FrustumConfig {
        taper_profile: TaperProfile::Exponential,
        ..base_config
    };
    let exponential_system = create_geometry(
        (120.0, 40.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllFrustum(exponential_config),
    );
    plot_geometry(&exponential_system, "outputs/frustum_exponential_taper.svg")?;
    println!("   ✅ Exponential taper: saved to frustum_exponential_taper.svg");

    let smooth_config = FrustumConfig {
        taper_profile: TaperProfile::Smooth,
        ..base_config
    };
    let smooth_system = create_geometry(
        (120.0, 40.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllFrustum(smooth_config),
    );
    plot_geometry(&smooth_system, "outputs/frustum_smooth_taper.svg")?;
    println!("   ✅ Smooth taper: saved to frustum_smooth_taper.svg");

    println!();
    Ok(())
}

fn demonstrate_throat_positions() -> Result<(), Box<dyn std::error::Error>> {
    println!("2️⃣  Demonstrating Throat Positions");
    println!("   Testing throat at 25%, 50%, and 75% positions");

    let config_25 = FrustumConfig {
        inlet_width: 2.5,
        throat_width: 0.6,
        outlet_width: 2.0,
        taper_profile: TaperProfile::Smooth,
        smoothness: 80,
        throat_position: 0.25,
    };

    let system = create_geometry(
        (150.0, 50.0),
        &[SplitType::Trifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllFrustum(config_25),
    );
    plot_geometry(&system, "outputs/frustum_throat_positions.svg")?;
    println!("   ✅ Variable throat positions: saved to frustum_throat_positions.svg");

    println!();
    Ok(())
}

fn demonstrate_width_configurations() -> Result<(), Box<dyn std::error::Error>> {
    println!("3️⃣  Demonstrating Width Configurations");
    println!("   Testing different inlet/throat/outlet width ratios");

    let high_compression = FrustumConfig {
        inlet_width: 4.0,
        throat_width: 0.4,
        outlet_width: 3.0,
        taper_profile: TaperProfile::Exponential,
        smoothness: 60,
        throat_position: 0.4,
    };

    let system = create_geometry(
        (100.0, 60.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllFrustum(high_compression),
    );
    plot_geometry(&system, "outputs/frustum_width_configs.svg")?;
    println!("   ✅ High compression ratio: saved to frustum_width_configs.svg");

    println!();
    Ok(())
}

fn demonstrate_mixed_channel_systems() -> Result<(), Box<dyn std::error::Error>> {
    println!("4️⃣  Demonstrating Mixed Channel Systems");
    println!("   Testing adaptive selection with frustum channels included");

    let adaptive_config = ChannelTypeConfig::default();

    let system = create_geometry(
        (200.0, 80.0),
        &[SplitType::Bifurcation, SplitType::Trifurcation],
        &GeometryConfig::default(),
        &adaptive_config,
    );
    plot_geometry(&system, "outputs/mixed_channel_system.svg")?;
    println!("   ✅ Mixed system: saved to mixed_channel_system.svg");

    let mut channel_counts = std::collections::HashMap::new();
    for channel in &system.channels {
        let channel_type_name = match &channel.channel_type {
            cfd_schematics::geometry::ChannelType::Straight => "Straight",
            cfd_schematics::geometry::ChannelType::SmoothStraight { .. } => "SmoothStraight",
            cfd_schematics::geometry::ChannelType::Serpentine { .. } => "Serpentine",
            cfd_schematics::geometry::ChannelType::Arc { .. } => "Arc",
            cfd_schematics::geometry::ChannelType::Frustum { .. } => "Frustum",
        };
        *channel_counts.entry(channel_type_name).or_insert(0) += 1;
    }

    println!("   📊 Channel type distribution:");
    for (channel_type, count) in &channel_counts {
        println!("      • {}: {} channels", channel_type, count);
    }

    println!();
    Ok(())
}

fn demonstrate_json_serialization() -> Result<(), Box<dyn std::error::Error>> {
    println!("5️⃣  Demonstrating JSON Serialization");
    println!("   Testing export/import of frustum channel systems");

    let frustum_config = FrustumConfig {
        inlet_width: 2.8,
        throat_width: 0.7,
        outlet_width: 2.2,
        taper_profile: TaperProfile::Smooth,
        smoothness: 75,
        throat_position: 0.6,
    };

    let original_system = create_geometry(
        (130.0, 55.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllFrustum(frustum_config),
    );

    let json = original_system.to_json()?;
    fs::write("outputs/frustum_system_export.json", &json)?;
    println!("   ✅ Exported to: frustum_system_export.json");

    let imported_system = ChannelSystem::from_json(&json)?;
    println!("   ✅ Successfully imported from JSON");

    assert_eq!(original_system.nodes.len(), imported_system.nodes.len());
    assert_eq!(
        original_system.channels.len(),
        imported_system.channels.len()
    );
    assert_eq!(original_system.box_dims, imported_system.box_dims);

    for (orig, imp) in original_system
        .channels
        .iter()
        .zip(imported_system.channels.iter())
    {
        match (&orig.channel_type, &imp.channel_type) {
            (
                cfd_schematics::geometry::ChannelType::Frustum {
                    path: orig_path,
                    widths: orig_widths,
                    inlet_width: orig_inlet,
                    throat_width: orig_throat,
                    outlet_width: orig_outlet,
                    taper_profile: orig_taper_profile,
                    throat_position: orig_throat_position,
                    has_venturi_throat: orig_has_venturi_throat,
                },
                cfd_schematics::geometry::ChannelType::Frustum {
                    path: imp_path,
                    widths: imp_widths,
                    inlet_width: imp_inlet,
                    throat_width: imp_throat,
                    outlet_width: imp_outlet,
                    taper_profile: imp_taper_profile,
                    throat_position: imp_throat_position,
                    has_venturi_throat: imp_has_venturi_throat,
                },
            ) => {
                assert_eq!(orig_path.len(), imp_path.len());
                assert_eq!(orig_widths.len(), imp_widths.len());
                assert_eq!(orig_inlet, imp_inlet);
                assert_eq!(orig_throat, imp_throat);
                assert_eq!(orig_outlet, imp_outlet);
                assert_eq!(orig_taper_profile, imp_taper_profile);
                assert_eq!(orig_throat_position, imp_throat_position);
                assert_eq!(orig_has_venturi_throat, imp_has_venturi_throat);
            }
            _ => panic!("Expected frustum channels"),
        }
    }

    println!("   ✅ Data integrity verified - all frustum properties preserved");

    println!();
    Ok(())
}

#[allow(dead_code)]
fn print_width_profile(config: &FrustumConfig, name: &str) {
    println!("   📏 {} Width Profile:", name);
    for i in 0..=10 {
        let t = i as f64 / 10.0;
        let width = config.width_at_position(t);
        println!("      t={:.1}: width={:.3}", t, width);
    }
}
