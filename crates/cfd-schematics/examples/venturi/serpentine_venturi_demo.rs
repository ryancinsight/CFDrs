//! Serpentine-with-Venturi Demonstration
//!
//! Demonstrates the `serpentine_venturi_rect` preset which places venturi
//! constrictions at each U-turn bend of a serpentine channel — the Dean-flow
//! apex where centripetal acceleration is maximum.
//!
//! Generates several configurations and renders auto-annotated schematics
//! showing venturi throat markers (TH1, TH2, …) at each bend.
//!
//! # Physics
//!
//! At a curved-channel U-turn the centripetal acceleration is:
//!
//! ```text
//! a_c = u² / R_bend
//! ```
//!
//! This drives Dean secondary vortices (De = Re·√(D_h / 2R)) that
//! pre-focus cells toward the centreline before they enter the venturi
//! throat.  Placing the throat at the Dean apex maximises both cell
//! concentration and wall shear — ideal for CTC cavitation in SDT.

use cfd_schematics::interface::presets::serpentine_venturi_rect;
#[path = "../shared/mod.rs"]
mod shared;

#[derive(Clone, Copy)]
struct VenturiGeometry {
    label: &'static str,
    segments: usize,
    segment_length_m: f64,
    channel_width_m: f64,
    throat_width_m: f64,
    channel_height_m: f64,
    throat_length_m: f64,
    bend_radius_m: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ── Configuration matrix ────────────────────────────────────────
    // Millifluidic-scale channels within 96-well-plate footprint
    // (127.76 × 85.47 mm — ANSI/SLAS 1-2004).
    let configurations = [
        VenturiGeometry {
            label: "baseline_4seg",
            segments: 4,
            segment_length_m: 0.015,
            channel_width_m: 0.002,
            throat_width_m: 0.001,
            channel_height_m: 0.001,
            throat_length_m: 0.003,
            bend_radius_m: 0.001,
        },
        VenturiGeometry {
            label: "aggressive_6seg",
            segments: 6,
            segment_length_m: 0.012,
            channel_width_m: 0.003,
            throat_width_m: 0.0012,
            channel_height_m: 0.001,
            throat_length_m: 0.004,
            bend_radius_m: 0.0015,
        },
        VenturiGeometry {
            label: "tight_bend_3seg",
            segments: 3,
            segment_length_m: 0.020,
            channel_width_m: 0.002,
            throat_width_m: 0.0008,
            channel_height_m: 0.001,
            throat_length_m: 0.002,
            bend_radius_m: 0.0005,
        },
        VenturiGeometry {
            label: "mild_8seg",
            segments: 8,
            segment_length_m: 0.010,
            channel_width_m: 0.002,
            throat_width_m: 0.0015,
            channel_height_m: 0.001,
            throat_length_m: 0.003,
            bend_radius_m: 0.002,
        },
    ];

    println!("Serpentine-Venturi Dean-Apex Demonstration");
    println!("==========================================\n");

    for configuration in configurations {
        let bp = serpentine_venturi_rect(
            format!("sv_{}", configuration.label),
            configuration.segments,
            configuration.segment_length_m,
            configuration.channel_width_m,
            configuration.throat_width_m,
            configuration.channel_height_m,
            configuration.throat_length_m,
            configuration.bend_radius_m,
        );

        let n_venturis = bp.venturi_channels().len();
        let n_segments = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        let total_channels = bp.channels.len();

        // Compute Dean number at each bend for this geometry
        let d_h = 2.0 * configuration.channel_width_m * configuration.channel_height_m
            / (configuration.channel_width_m + configuration.channel_height_m);
        // Re at throat (blood: ρ=1060 kg/m³, μ=3.5e-3 Pa·s)
        let rho = 1060.0_f64;
        let mu = 3.5e-3_f64;
        // Approximate mean velocity for Re ~ 100 (typical millifluidic)
        let re_channel = 100.0_f64;
        let de = re_channel * (d_h / (2.0 * configuration.bend_radius_m)).sqrt();
        let contraction_ratio = configuration.channel_width_m / configuration.throat_width_m;

        println!("─── {} ───", configuration.label);
        println!("  Segments:           {n_segments}");
        println!("  Venturi throats:    {n_venturis}");
        println!("  Total channels:     {total_channels}");
        println!(
            "  Channel width:      {:.1} mm",
            configuration.channel_width_m * 1e3
        );
        println!(
            "  Throat width:       {:.2} mm",
            configuration.throat_width_m * 1e3
        );
        println!("  Contraction ratio:  {contraction_ratio:.2}:1");
        println!(
            "  Bend radius:        {:.1} mm",
            configuration.bend_radius_m * 1e3
        );
        println!("  D_h:                {:.3} mm", d_h * 1e3);
        println!("  Dean number (Re=100): {de:.1}");
        println!(
            "  Centripetal accel:  a_c = u²/R = (ν·Re/D_h)² / R = {:.2} m/s²",
            {
                let u = re_channel * mu / (rho * d_h);
                u * u / configuration.bend_radius_m
            }
        );

        // Render schematic
        shared::save_example_output_with_name(&bp, "serpentine_venturi_demo", configuration.label);
    }

    // ── Summary comparison table ────────────────────────────────────
    println!("Summary:");
    println!(
        "{:<22} {:>4} {:>5} {:>6} {:>5} {:>5}",
        "Config", "Segs", "Vents", "De", "CR", "R(mm)"
    );
    println!("{}", "─".repeat(52));
    for configuration in configurations {
        let d_h = 2.0 * configuration.channel_width_m * configuration.channel_height_m
            / (configuration.channel_width_m + configuration.channel_height_m);
        let de = 100.0 * (d_h / (2.0 * configuration.bend_radius_m)).sqrt();
        let cr = configuration.channel_width_m / configuration.throat_width_m;
        println!(
            "{:<22} {:>4} {:>5} {:>6.1} {:>5.2} {:>5.1}",
            configuration.label,
            configuration.segments,
            configuration.segments - 1,
            de,
            cr,
            configuration.bend_radius_m * 1e3,
        );
    }

    Ok(())
}
