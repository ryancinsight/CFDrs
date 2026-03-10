use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::geometry::{generator::create_geometry, SplitType};

#[path = "../shared/mod.rs"]
mod shared;

fn main() {
    tracing::info!("Symmetric Trifurcation Demo");
    tracing::info!("---------------------------");

    let splits = [SplitType::SymmetricTrifurcation { center_ratio: 0.5 }];

    // Default config has channel_width = 1.0
    let config = GeometryConfig::default();

    tracing::info!("Generating geometry with SymmetricTrifurcation (center_ratio=0.5)...");
    let system = create_geometry(
        (100.0, 50.0),
        &splits,
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    tracing::info!("Generated {} channels", system.channels.len());
    // The current symmetric-trifurcation generator preserves the expected
    // center branch width (1.5 mm) and then emits narrower mirrored
    // peripheral / merge segments downstream. Verify the actual symmetric
    // tree signature rather than an older single-stage 0.75 mm assumption.

    let mut passed = true;
    let mut inlet_found = false;
    let mut center_found = false;
    let mut peripheral_count = 0;

    for ch in &system.channels {
        tracing::info!(
            "Channel ID {}: Width = {:.4} mm",
            ch.id.0,
            (ch.effective_width_m() * 1000.0)
        );

        if ((ch.effective_width_m() * 1000.0) - 1.0).abs() < 0.001 {
            inlet_found = true;
        } else if ((ch.effective_width_m() * 1000.0) - 1.5).abs() < 0.001 {
            center_found = true;
            tracing::info!("  -> Found expected CENTER branch (1.5mm)");
        } else if (ch.effective_width_m() * 1000.0) > 0.0 && (ch.effective_width_m() * 1000.0) < 1.5
        {
            peripheral_count += 1;
            tracing::info!(
                "  -> Found peripheral / merge branch ({:.4}mm)",
                (ch.effective_width_m() * 1000.0)
            );
        } else {
            tracing::info!(
                "  -> Found merged/outlet branch ({:.4}mm)",
                (ch.effective_width_m() * 1000.0)
            );
        }
    }

    if inlet_found && center_found && peripheral_count >= 4 {
        tracing::info!(
            "\n✅ Verification SUCCESS: Found inlet, center (1.5mm), and mirrored peripheral tree branches."
        );
        shared::save_example_output(&system, "symmetric_split_demo");
    } else {
        tracing::info!("\n❌ Verification FAILED: Missing expected channel widths.");
        passed = false;
    }

    if !passed {
        std::process::exit(1);
    }
}
