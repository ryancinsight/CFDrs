//! Diverse top-5 final designs — one best-in-class design per optimization mode.
//!
//! Selects the rank-1 design from each of 5 distinct [`OptimMode`] variants,
//! validates hard constraints (FDA shear, SBS-96 footprint, hemolysis safety),
//! generates full mesh artifacts via [`DesignPipeline`], and prints a
//! side-by-side comparison table.
//!
//! # Output
//!
//! ```text
//! crates/cfd-optim/outputs/final_designs/
//!   sdt_cavitation/       — fluid.stl, chip.stl, schematic.svg, schematic.json, constant/polyMesh/
//!   sdt_therapy/          — ...
//!   hydro_cavitation_sdt/ — ...
//!   rbc_protected_sdt/    — ...
//!   combined_sdt_leuka/   — ...
//! ```
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_top5_final --features mesh-export
//! ```

use std::path::Path;

use cfd_optim::{
    constraints::FDA_MAX_WALL_SHEAR_PA,
    pipeline::{DesignArtifacts, DesignPipeline},
    score_description, OptimMode, RankedDesign, SdtOptimizer, SdtWeights,
};

/// The 5 optimization modes selected to produce maximally diverse final designs.
///
/// Each mode focuses on a different clinical / engineering objective, ensuring
/// the final set covers cavitation intensity, selective therapy, RBC safety,
/// combined leukapheresis, and overall SDT therapy quality.
fn diverse_modes() -> Vec<(&'static str, OptimMode)> {
    vec![
        ("sdt_cavitation", OptimMode::SdtCavitation),
        ("sdt_therapy", OptimMode::SdtTherapy),
        (
            "hydro_cavitation_sdt",
            OptimMode::HydrodynamicCavitationSDT,
        ),
        ("rbc_protected_sdt", OptimMode::RbcProtectedSdt),
        (
            "combined_sdt_leuka",
            OptimMode::CombinedSdtLeukapheresis {
                leuka_weight: 0.5,
                sdt_weight: 0.5,
                patient_weight_kg: 3.0,
            },
        ),
    ]
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("final_designs");
    std::fs::create_dir_all(&out_root)?;

    let weights = SdtWeights::default();
    let modes = diverse_modes();
    let pipeline = DesignPipeline::new(&out_root);

    println!();
    println!("{}", "=".repeat(120));
    println!(
        "  sdt_top5_final  |  Best-in-class design from each of {} optimization modes",
        modes.len()
    );
    println!("  Constraints: FDA main-channel shear <= {FDA_MAX_WALL_SHEAR_PA} Pa  |  SBS-96 plate fit  |  hemolysis safe");
    println!("{}", "=".repeat(120));

    // ── Phase 1: Optimize, validate, and export ─────────────────────────────

    let mut results: Vec<(&str, RankedDesign, DesignArtifacts)> = Vec::new();

    for (slug, mode) in &modes {
        println!();
        println!("{}", "-".repeat(120));
        println!(
            "  Mode: {slug:<25}  ({})",
            score_description(*mode)
        );
        println!("{}", "-".repeat(120));

        // Get top-1 for this mode
        let top1 = match SdtOptimizer::new(*mode, weights).top_k(1) {
            Ok(designs) => designs,
            Err(e) => {
                eprintln!("  ERROR: optimization failed for {slug}: {e}");
                continue;
            }
        };

        let design = &top1[0];
        let m = &design.metrics;

        // ── Hard constraint validation ──────────────────────────────────────

        // FDA main-channel shear limit
        assert!(
            m.fda_main_compliant,
            "FDA shear violation for {slug}: main-channel shear {:.1} Pa exceeds {FDA_MAX_WALL_SHEAR_PA} Pa",
            m.max_main_channel_shear_pa,
        );

        // SBS-96 footprint
        assert!(
            m.plate_fits,
            "SBS-96 footprint violation for {slug}: design exceeds 127.76 x 85.47 mm plate boundary",
        );

        // Pressure feasibility
        assert!(
            m.pressure_feasible,
            "Pressure feasibility violation for {slug}: total dP {:.0} Pa exceeds available gauge",
            m.total_pressure_drop_pa,
        );

        // Hemolysis safety for therapy modes that enforce it
        if matches!(
            mode,
            OptimMode::SdtTherapy
                | OptimMode::HydrodynamicCavitationSDT
                | OptimMode::CombinedSdtLeukapheresis { .. }
        ) {
            assert!(
                m.hemolysis_index_per_pass < cfd_optim::constraints::HI_PASS_LIMIT,
                "Hemolysis limit exceeded for {slug}: HI={:.2e} > limit {:.2e}",
                m.hemolysis_index_per_pass,
                cfd_optim::constraints::HI_PASS_LIMIT,
            );
        }

        println!(
            "  Validated: FDA={} plate={} dP_ok={} HI={:.2e}",
            if m.fda_main_compliant { "OK" } else { "!!" },
            if m.plate_fits { "OK" } else { "!!" },
            if m.pressure_feasible { "OK" } else { "!!" },
            m.hemolysis_index_per_pass,
        );

        // ── Mesh artifact export ────────────────────────────────────────────

        let artifacts = match pipeline.export_design(design) {
            Ok(a) => a,
            Err(e) => {
                eprintln!("  ERROR: mesh export failed for {slug}: {e}");
                continue;
            }
        };

        println!(
            "  Mesh: {} vertices, {} faces, watertight={}",
            artifacts.vertex_count, artifacts.face_count, artifacts.watertight,
        );
        println!("  STL:  {}", artifacts.fluid_stl.display());
        if let Some(ref chip) = artifacts.chip_stl {
            println!("  Chip: {}", chip.display());
        }
        println!("  OF:   {}", artifacts.openfoam_dir.display());
        println!("  SVG:  {}", artifacts.schematic_svg.display());
        println!("  JSON: {}", artifacts.schematic_json.display());

        results.push((slug, design.clone(), artifacts));
    }

    // ── Phase 2: Comparison table ───────────────────────────────────────────

    println!();
    println!("{}", "=".repeat(120));
    println!("  COMPARISON TABLE — Best-in-class from {} modes", results.len());
    println!("{}", "=".repeat(120));
    println!(
        "  {:<25} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "Mode", "Score", "Sigma", "Shear", "HI/pass", "Sep3", "Cov%", "dP kPa", "Verts", "WT"
    );
    println!("  {}", "-".repeat(115));

    for (slug, design, artifacts) in &results {
        let m = &design.metrics;
        let sigma = if m.cavitation_number.is_finite() {
            format!("{:.3}", m.cavitation_number)
        } else {
            "inf".to_string()
        };
        println!(
            "  {:<25} {:>7.4} {:>8} {:>7.1} {:>8.2e} {:>8.4} {:>7.0}% {:>8.1} {:>8} {:>8}",
            slug,
            design.score,
            sigma,
            m.max_main_channel_shear_pa,
            m.hemolysis_index_per_pass,
            m.three_pop_sep_efficiency,
            m.well_coverage_fraction * 100.0,
            m.total_pressure_drop_pa * 1e-3,
            artifacts.vertex_count,
            if artifacts.watertight { "yes" } else { "no" },
        );
    }
    println!("  {}", "-".repeat(115));

    // ── Phase 3: Extended metrics table ──────────────────────────────────────

    println!();
    println!(
        "  {:<25} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "Mode", "CancCav", "RBCprot", "Sono", "OncSel", "ThWindow", "LysRisk", "FDA_all"
    );
    println!("  {}", "-".repeat(95));

    for (slug, design, _) in &results {
        let m = &design.metrics;
        println!(
            "  {:<25} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.2e} {:>8}",
            slug,
            m.cancer_targeted_cavitation,
            m.rbc_venturi_protection,
            m.sonoluminescence_proxy,
            m.oncology_selectivity_index,
            m.therapeutic_window_score,
            m.lysis_risk_index,
            if m.fda_overall_compliant { "yes" } else { "no" },
        );
    }
    println!("  {}", "-".repeat(95));

    println!();
    println!("  Output root: {}", out_root.display());
    println!(
        "  {} designs exported with full mesh artifacts.",
        results.len()
    );
    println!();

    Ok(())
}
