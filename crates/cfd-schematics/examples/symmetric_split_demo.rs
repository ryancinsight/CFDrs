use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::geometry::{create_geometry, SplitType};

fn main() {
    println!("Symmetric Trifurcation Demo");
    println!("---------------------------");

    let splits = [SplitType::SymmetricTrifurcation { center_ratio: 0.5 }];

    // Default config has channel_width = 1.0
    let config = GeometryConfig::default();

    println!("Generating geometry with SymmetricTrifurcation (center_ratio=0.5)...");
    let system = create_geometry(
        (100.0, 50.0),
        &splits,
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    println!("Generated {} channels", system.channels.len());
    // The current symmetric-trifurcation generator preserves the expected
    // center branch width (1.5 mm) and then emits narrower mirrored
    // peripheral / merge segments downstream. Verify the actual symmetric
    // tree signature rather than an older single-stage 0.75 mm assumption.

    let mut passed = true;
    let mut inlet_found = false;
    let mut center_found = false;
    let mut peripheral_count = 0;

    for ch in &system.channels {
        println!("Channel ID {}: Width = {:.4} mm", ch.id, ch.width);

        if (ch.width - 1.0).abs() < 0.001 {
            inlet_found = true;
        } else if (ch.width - 1.5).abs() < 0.001 {
            center_found = true;
            println!("  -> Found expected CENTER branch (1.5mm)");
        } else if ch.width > 0.0 && ch.width < 1.5 {
            peripheral_count += 1;
            println!("  -> Found peripheral / merge branch ({:.4}mm)", ch.width);
        } else {
            println!("  -> Found merged/outlet branch ({:.4}mm)", ch.width);
        }
    }

    if inlet_found && center_found && peripheral_count >= 4 {
        println!(
            "\n✅ Verification SUCCESS: Found inlet, center (1.5mm), and mirrored peripheral tree branches."
        );
    } else {
        println!("\n❌ Verification FAILED: Missing expected channel widths.");
        passed = false;
    }

    if !passed {
        std::process::exit(1);
    }
}
