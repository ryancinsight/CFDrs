//! Regenerate the three schematic SVG figures used in the Milestone 12 report.
//!
//! Outputs are written to `report/figures/` at the workspace root.
//!
//! The three figures are:
//! - **Figure 4** – `selected_ga_schematic.svg`  
//!   Option 1 selected design (ultrasound, no venturi):
//!   `151254-TS-q100ml-g100kPa-dt0um-tl2-w8000um-h2000-n18`  
//!   (TrifurcationSerpentine, 100 mL/min, 100 kPa, 8000 µm width, 2000 µm height)
//!
//! - **Figure 5** – `selected_cifx_combined_schematic.svg`  
//!   Combined SDT+leukapheresis top design:
//!   `405279-CIFX-pt3-pcf540-tcf530-btf850-q100ml-g200kPa-dt50um-tl10-w6000um`  
//!   (IncrementalFiltration 3-pre-tri, 100 mL/min, 200 kPa, 50 µm throat, tl×10)
//!
//! - **Figure 6** – `selected_cif_schematic.svg`  
//!   RBC-protected SDT top design:
//!   `389191-CIFX-pt2-pcf580-tcf560-btf850-q100ml-g400kPa-dt45um-tl5-w6000um`  
//!   (IncrementalFiltration 2-pre-tri, 100 mL/min, 400 kPa, 45 µm throat, tl×5)
//!
//! All schematics are rendered with `RenderConfig::well_plate_96()` (1278×855 px,
//! 127.76 × 85.47 mm axis range) via `map_to_plate_coords()` in `blueprint.rs`.
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example regen_report_schematics --no-default-features
//! ```

use cfd_optim::{export::save_schematic_svg, CrossSectionShape, DesignCandidate, DesignTopology};
use std::path::PathBuf;

fn main() {
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("cfd-optim crate has a parent")
        .parent()
        .expect("crates/ has a workspace root")
        .to_path_buf();
    let figures_dir = workspace_root.join("report").join("figures");
    std::fs::create_dir_all(&figures_dir).expect("report/figures/ must be creatable");

    // ── Figure 4 — Option 1 selected no-venturi branch network ─────────────
    // Design ID: 151254-TS-q100ml-g100kPa-dt0um-tl2-w8000um-h2000-n18
    // Topology:  TrifurcationSerpentine
    let option1_ts = DesignCandidate {
        id: "151254-TS-q100ml-g100kPa-dt0um-tl2-w8000um-h2000-n18".to_string(),
        topology: DesignTopology::TrifurcationSerpentine,
        flow_rate_m3_s: 1.667e-6,  // 100 mL/min
        inlet_gauge_pa: 100_000.0, // 100 kPa
        throat_diameter_m: 0.0,    // no venturi in Option 1
        inlet_diameter_m: 4e-3,
        throat_length_m: 1e-6,  // placeholder (unused for no-venturi topology)
        channel_width_m: 8e-3,  // 8000 µm
        channel_height_m: 2e-3, // 2000 µm
        serpentine_segments: 18,
        segment_length_m: 0.045, // 45 mm = treatment zone width
        bend_radius_m: 4.5e-3,   // half well pitch
        feed_hematocrit: 0.45,
        trifurcation_center_frac: 1.0 / 3.0,
        cif_pretri_center_frac: 1.0 / 3.0,
        cif_terminal_tri_center_frac: 1.0 / 3.0,
        cif_terminal_bi_treat_frac: 0.68,
        asymmetric_narrow_frac: 0.5,
        trifurcation_left_frac: 1.0 / 3.0,
        cross_section_shape: CrossSectionShape::Rectangular,
    };

    // ── Figure 5 — combined SDT+leukapheresis CIFX-pt3 ─────────────────────
    // Design ID: 405279-CIFX-pt3-pcf540-tcf530-btf850-q100ml-g200kPa-dt50um-tl10-w6000um
    // Topology:  IncrementalFiltrationTriBiSeparator { n_pretri: 3 }
    let cifx_pt3_combined = DesignCandidate {
        id: "405279-CIFX-pt3-pcf540-tcf530-btf850-q100ml-g200kPa-dt50um-tl10-w6000um".to_string(),
        topology: DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri: 3 },
        flow_rate_m3_s: 1.667e-6,  // 100 mL/min
        inlet_gauge_pa: 200_000.0, // 200 kPa
        throat_diameter_m: 50e-6,  // 50 µm
        inlet_diameter_m: 4e-3,
        throat_length_m: 500e-6, // tl10 × 50 µm
        channel_width_m: 6e-3,
        channel_height_m: 1e-3,
        serpentine_segments: 1,
        segment_length_m: 0.045,
        bend_radius_m: 4.5e-3,
        feed_hematocrit: 0.45,
        trifurcation_center_frac: 0.54,     // same as pcf540
        cif_pretri_center_frac: 0.54,       // pcf540
        cif_terminal_tri_center_frac: 0.53, // tcf530
        cif_terminal_bi_treat_frac: 0.85,   // btf850
        asymmetric_narrow_frac: 0.5,
        trifurcation_left_frac: 1.0 / 3.0,
        cross_section_shape: CrossSectionShape::Rectangular,
    };

    // ── Figure 6 — RBC-protected SDT CIFX-pt2 ──────────────────────────────
    // Design ID: 389191-CIFX-pt2-pcf580-tcf560-btf850-q100ml-g400kPa-dt45um-tl5-w6000um
    // Topology:  IncrementalFiltrationTriBiSeparator { n_pretri: 2 }
    let cifx_pt2_rbc = DesignCandidate {
        id: "389191-CIFX-pt2-pcf580-tcf560-btf850-q100ml-g400kPa-dt45um-tl5-w6000um".to_string(),
        topology: DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri: 2 },
        flow_rate_m3_s: 1.667e-6,  // 100 mL/min
        inlet_gauge_pa: 400_000.0, // 400 kPa
        throat_diameter_m: 45e-6,  // 45 µm
        inlet_diameter_m: 4e-3,
        throat_length_m: 225e-6, // tl5 × 45 µm
        channel_width_m: 6e-3,
        channel_height_m: 1e-3,
        serpentine_segments: 1,
        segment_length_m: 0.045,
        bend_radius_m: 4.5e-3,
        feed_hematocrit: 0.45,
        trifurcation_center_frac: 0.58,     // same as pcf580
        cif_pretri_center_frac: 0.58,       // pcf580
        cif_terminal_tri_center_frac: 0.56, // tcf560
        cif_terminal_bi_treat_frac: 0.85,   // btf850
        asymmetric_narrow_frac: 0.5,
        trifurcation_left_frac: 1.0 / 3.0,
        cross_section_shape: CrossSectionShape::Rectangular,
    };

    // ── Render all three figures ────────────────────────────────────────────
    let figures = [
        (
            &option1_ts,
            "selected_ga_schematic.svg",
            "Figure 4 (Option 1 no-venturi)",
        ),
        (
            &cifx_pt3_combined,
            "selected_cifx_combined_schematic.svg",
            "Figure 5 (CIFX-pt3 combined)",
        ),
        (
            &cifx_pt2_rbc,
            "selected_cif_schematic.svg",
            "Figure 6 (CIFX-pt2 RBC-protected)",
        ),
    ];

    let mut success = 0usize;
    let mut failures = 0usize;

    for (candidate, filename, label) in &figures {
        let path = figures_dir.join(filename);
        match save_schematic_svg(candidate, &path) {
            Ok(()) => {
                println!("  ✓ {label}  →  {}", path.display());
                success += 1;
            }
            Err(e) => {
                eprintln!("  ✗ {label}  FAILED: {e}");
                failures += 1;
            }
        }
    }

    println!();
    println!("  {success} figure(s) generated, {failures} failure(s).");
    if failures > 0 {
        std::process::exit(1);
    }
}
