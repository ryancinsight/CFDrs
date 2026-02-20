use cfd_schematics::geometry::{create_geometry, SplitType};
use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};

fn main() {
    println!("Asymmetric Bifurcation Demo");
    println!("---------------------------");

    let splits = [
        SplitType::AsymmetricBifurcation { ratio: 0.7 }
    ];
    
    // Default config has channel_width = 1.0
    let config = GeometryConfig::default(); 
    
    println!("Generating geometry with AsymmetricBifurcation (ratio=0.7)...");
    let system = create_geometry(
        (100.0, 50.0), 
        &splits, 
        &config, 
        &ChannelTypeConfig::AllStraight
    );

    println!("Generated {} channels", system.channels.len());
    println!("Base Channel Width: {:.3} mm", config.channel_width);

    let mut passed = true;
    let mut inlet_found = false;
    let mut wide_branch_found = false;
    let mut narrow_branch_found = false;

    for ch in &system.channels {
        println!("Channel ID {}: Width = {:.3} mm, Type = {:?}", ch.id, ch.width, ch.channel_type);
        
        // Simple heuristic verification
        if (ch.width - 1.0).abs() < 0.001 {
            inlet_found = true;
        } else if (ch.width - 1.4).abs() < 0.001 {
            wide_branch_found = true;
            println!("  -> Found expected WIDE branch (1.4mm)");
        } else if (ch.width - 0.6).abs() < 0.001 {
            narrow_branch_found = true;
            println!("  -> Found expected NARROW branch (0.6mm)");
        } else {
            // Merge channels might have average width?
            // (1.4 + 0.6) / 2 = 1.0. So outlet/merge channels should be 1.0.
            if (ch.width - 1.0).abs() < 0.001 {
                 println!("  -> Found expected merged/outlet branch (1.0mm)");
            } else {
                 println!("  -> Unexpected width!");
                 passed = false;
            }
        }
    }

    if inlet_found && wide_branch_found && narrow_branch_found {
        println!("\n✅ Verification SUCCESS: All expected channel widths found.");
    } else {
        println!("\n❌ Verification FAILED: Missing expected channel widths.");
        passed = false;
    }

    if !passed {
        std::process::exit(1);
    }
}
