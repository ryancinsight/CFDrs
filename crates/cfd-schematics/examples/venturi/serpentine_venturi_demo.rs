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
use cfd_schematics::visualizations::schematic::plot_blueprint_auto_annotated;
use cfd_schematics::visualizations::traits::RenderConfig;
#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ── Configuration matrix ────────────────────────────────────────
    // Millifluidic-scale channels within 96-well-plate footprint
    // (127.76 × 85.47 mm — ANSI/SLAS 1-2004).
    let configs: Vec<(&str, usize, f64, f64, f64, f64, f64, f64)> = vec![
        // (label, segments, seg_len_m, width_m, throat_m, height_m, throat_len_m, bend_r_m)
        (
            "baseline_4seg",
            4,
            0.015, // 15 mm segments
            0.002, // 2 mm channel width
            0.001, // 1 mm throat
            0.001, // 1 mm height
            0.003, // 3 mm throat length
            0.001, // 1 mm bend radius
        ),
        (
            "aggressive_6seg",
            6,
            0.012,  // 12 mm segments
            0.003,  // 3 mm channel width
            0.0012, // 1.2 mm throat (contraction ratio 2.5:1)
            0.001,  // 1 mm height
            0.004,  // 4 mm throat length
            0.0015, // 1.5 mm bend radius
        ),
        (
            "tight_bend_3seg",
            3,
            0.020,  // 20 mm segments
            0.002,  // 2 mm channel width
            0.0008, // 0.8 mm throat (contraction ratio 2.5:1)
            0.001,  // 1 mm height
            0.002,  // 2 mm throat length
            0.0005, // 0.5 mm bend radius — high Dean number
        ),
        (
            "mild_8seg",
            8,
            0.010,  // 10 mm segments
            0.002,  // 2 mm channel width
            0.0015, // 1.5 mm throat (contraction ratio 1.33:1)
            0.001,  // 1 mm height
            0.003,  // 3 mm throat length
            0.002,  // 2 mm bend radius — lower Dean number
        ),
    ];

    println!("Serpentine-Venturi Dean-Apex Demonstration");
    println!("==========================================\n");

    for (label, segments, seg_len, width, throat, height, throat_len, bend_r) in &configs {
        let bp = serpentine_venturi_rect(
            format!("sv_{label}"),
            *segments,
            *seg_len,
            *width,
            *throat,
            *height,
            *throat_len,
            *bend_r,
        );

        let n_venturis = bp.venturi_channels().len();
        let n_segments = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        let total_channels = bp.channels.len();

        // Compute Dean number at each bend for this geometry
        let d_h = 2.0 * width * height / (width + height);
        // Re at throat (blood: ρ=1060 kg/m³, μ=3.5e-3 Pa·s)
        let rho = 1060.0_f64;
        let mu = 3.5e-3_f64;
        // Approximate mean velocity for Re ~ 100 (typical millifluidic)
        let re_channel = 100.0_f64;
        let de = re_channel * (d_h / (2.0 * bend_r)).sqrt();
        let contraction_ratio = width / throat;

        println!("─── {label} ───");
        println!("  Segments:           {n_segments}");
        println!("  Venturi throats:    {n_venturis}");
        println!("  Total channels:     {total_channels}");
        println!("  Channel width:      {:.1} mm", width * 1e3);
        println!("  Throat width:       {:.2} mm", throat * 1e3);
        println!("  Contraction ratio:  {contraction_ratio:.2}:1");
        println!("  Bend radius:        {:.1} mm", bend_r * 1e3);
        println!("  D_h:                {:.3} mm", d_h * 1e3);
        println!("  Dean number (Re=100): {de:.1}");
        println!(
            "  Centripetal accel:  a_c = u²/R = (ν·Re/D_h)² / R = {:.2} m/s²",
            {
                let u = re_channel * mu / (rho * d_h);
                u * u / bend_r
            }
        );

        // Render schematic
        shared::save_example_output_with_name(&bp, "serpentine_venturi_demo", label);
    }

    // ── Summary comparison table ────────────────────────────────────
    println!("Summary:");
    println!(
        "{:<22} {:>4} {:>5} {:>6} {:>5} {:>5}",
        "Config", "Segs", "Vents", "De", "CR", "R(mm)"
    );
    println!("{}", "─".repeat(52));
    for (label, segments, _, width, throat, height, _, bend_r) in &configs {
        let d_h = 2.0 * width * height / (width + height);
        let de = 100.0 * (d_h / (2.0 * bend_r)).sqrt();
        let cr = width / throat;
        println!(
            "{:<22} {:>4} {:>5} {:>6.1} {:>5.2} {:>5.1}",
            label,
            segments,
            segments - 1,
            de,
            cr,
            bend_r * 1e3,
        );
    }

    Ok(())
}
