use cfd_schematics::geometry::{create_geometry, SplitType};
use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};

fn main() {
    println!("Symmetric Trifurcation Demo");
    println!("---------------------------");

    let splits = [
        SplitType::SymmetricTrifurcation { center_ratio: 0.5 }
    ];
    
    // Default config has channel_width = 1.0
    let config = GeometryConfig::default(); 
    
    println!("Generating geometry with SymmetricTrifurcation (center_ratio=0.5)...");
    let system = create_geometry(
        (100.0, 50.0), 
        &splits, 
        &config, 
        &ChannelTypeConfig::AllStraight
    );

    println!("Generated {} channels", system.channels.len());
    // Total capacity assumption: 3 * parent.
    // Parent = 1.0mm. Total = 3.0mm.
    // Center = 0.5 * 3.0 = 1.5mm.
    // Sides = (1.0 - 0.5) * 3.0 / 2 = 0.75mm.

    let mut passed = true;
    let mut inlet_found = false;
    let mut center_found = false;
    let mut side_count = 0;

    for ch in &system.channels {
        println!("Channel ID {}: Width = {:.4} mm", ch.id, ch.width);
        
        if (ch.width - 1.0).abs() < 0.001 {
            inlet_found = true;
        } else if (ch.width - 1.5).abs() < 0.001 {
            center_found = true;
            println!("  -> Found expected CENTER branch (1.5mm)");
        } else if (ch.width - 0.75).abs() < 0.001 {
            side_count += 1;
            println!("  -> Found expected SIDE branch (0.75mm)");
        } else {
             println!("  -> Found merged/outlet branch ({:.4}mm)", ch.width);
        }
    }

    if inlet_found && center_found && side_count >= 2 {
        println!("\n✅ Verification SUCCESS: Found center (1.5mm) and 2+ symmetric sides (0.75mm).");
    } else {
        println!("\n❌ Verification FAILED: Missing expected channel widths.");
        passed = false;
    }

    if !passed {
        std::process::exit(1);
    }
}
