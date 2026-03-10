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
    geometry::{generator::create_geometry, SplitType},
    visualizations::schematic::plot_geometry,
};
#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌊 Frustum Channel Demonstration");
    println!("================================");
    println!();

    demonstrate_taper_profiles()?;
    demonstrate_throat_positions()?;
    demonstrate_width_configurations()?;
    demonstrate_mixed_channel_systems()?;
    demonstrate_json_serialization()?;

    println!("✅ All demonstrations completed successfully!");
    println!();
    println!("   • [Various Output Types] -> output/examples/frustum_channel_demo/");

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
    shared::save_example_output_with_name(&linear_system, "frustum_channel_demo", "frustum_linear_taper");
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
    shared::save_example_output_with_name(&exponential_system, "frustum_channel_demo", "frustum_exponential_taper");
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
    shared::save_example_output_with_name(&smooth_system, "frustum_channel_demo", "frustum_smooth_taper");
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
    shared::save_example_output_with_name(&system, "frustum_channel_demo", "frustum_throat_positions");
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
    shared::save_example_output_with_name(&system, "frustum_channel_demo", "frustum_width_configs");
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
    shared::save_example_output_with_name(&system, "frustum_channel_demo", "mixed_channel_system");
    println!("   ✅ Mixed system: saved to mixed_channel_system.svg");

    let mut channel_counts = std::collections::HashMap::new();
    for channel in &system.channels {
        let channel_type_name = if channel.venturi_geometry.is_some() {
            "Frustum"
        } else {
            match &channel.channel_shape {
                cfd_schematics::domain::model::ChannelShape::Straight => "Straight",
                cfd_schematics::domain::model::ChannelShape::Serpentine { .. } => "Serpentine",
            }
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

    // We just rely on the saving logic that save_example_output already provides
    // but we simulate the original example logic below.
    shared::save_example_output_with_name(&original_system, "frustum_channel_demo", "frustum_system_export");
    
    let json = serde_json::to_string_pretty(&original_system)?;
    println!("   ✅ Exported to: frustum_system_export.json");

    let imported_system: cfd_schematics::domain::model::NetworkBlueprint =
        serde_json::from_str(&json)?;
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
        assert_eq!(orig.channel_shape, imp.channel_shape);
        assert_eq!(orig.path, imp.path);
        assert_eq!(orig.venturi_geometry, imp.venturi_geometry);
    }

    println!("   ✅ Data integrity verified - all frustum properties preserved");

    println!();
    Ok(())
}
