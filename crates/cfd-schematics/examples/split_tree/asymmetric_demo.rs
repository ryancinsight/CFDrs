use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::geometry::{create_geometry, SplitType};

#[path = "../shared/mod.rs"]
mod shared;
use shared::save_example_output;

fn main() {
    println!("Asymmetric Bifurcation Demo");
    println!("---------------------------");

    let splits = [SplitType::AsymmetricBifurcation { ratio: 0.7 }];

    // Default config has channel_width = 1.0
    let config = GeometryConfig::default();

    println!("Generating geometry with AsymmetricBifurcation (ratio=0.7)...");
    let system = create_geometry(
        (100.0, 50.0),
        &splits,
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    println!("Generated {} channels", system.channels.len());
    println!("Base Channel Width: {:.3} mm", config.channel_width);

    let mut passed = true;
    let mut base_width_found = false;
    let mut wider_than_base_found = false;
    let mut narrower_than_base_found = false;
    let mut min_width = f64::INFINITY;
    let mut max_width = f64::NEG_INFINITY;

    for ch in &system.channels {
        println!(
            "Channel ID {}: Width = {:.3} mm, Type = {:?}",
            ch.id.0,
            (ch.effective_width_m() * 1000.0),
            ch.channel_shape
        );

        min_width = min_width.min(ch.effective_width_m() * 1000.0);
        max_width = max_width.max(ch.effective_width_m() * 1000.0);

        if ((ch.effective_width_m() * 1000.0) - config.channel_width).abs() < 0.05 {
            base_width_found = true;
            println!("  -> Found base-width trunk or merge channel");
        } else if (ch.effective_width_m() * 1000.0) > config.channel_width * 1.10 {
            wider_than_base_found = true;
            println!("  -> Found widened asymmetric branch");
        } else if (ch.effective_width_m() * 1000.0) < config.channel_width * 0.90 {
            narrower_than_base_found = true;
            println!("  -> Found narrowed asymmetric branch");
        } else {
            println!("  -> Width falls in transitional range after adaptive shaping");
        }
    }

    println!(
        "\nObserved width range: {:.3} mm .. {:.3} mm",
        min_width, max_width
    );

    save_example_output(&system, "asymmetric_demo");

    if base_width_found && wider_than_base_found && narrower_than_base_found {
        println!("\n✅ Verification SUCCESS: Asymmetric widening and narrowing are both present.");
    } else {
        println!("\n❌ Verification FAILED: Missing base, widened, or narrowed channels.");
        passed = false;
    }

    if !passed {
        std::process::exit(1);
    }
}
