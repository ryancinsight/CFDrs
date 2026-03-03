//! Master pipeline: optimise -> mesh -> compare.
//!
//! 1. Runs [`SdtOptimizer`] in [`OptimMode::SdtTherapy`] mode to find the top 5
//!    designs from the full parametric sweep.
//! 2. For each ranked design, generates mesh artifacts (fluid STL, chip STL,
//!    OpenFOAM polyMesh, schematic SVG, interchange JSON) via [`DesignPipeline`].
//! 3. Prints a comparison report summarising scores, key metrics, and mesh
//!    statistics for all five designs.
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_pipeline --features mesh-export
//! ```

use cfd_optim::{DesignArtifacts, DesignPipeline, OptimMode, SdtOptimizer, SdtWeights};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("pipeline");
    std::fs::create_dir_all(&out_dir)?;

    // ── 1. Optimise ─────────────────────────────────────────────────────────
    println!("=== SDT Pipeline: SdtTherapy mode ===\n");
    println!("Running parametric sweep ...");

    let optimizer = SdtOptimizer::new(OptimMode::SdtTherapy, SdtWeights::default());
    let top5 = optimizer.top_k(5)?;

    // ── 2. Print design summaries ───────────────────────────────────────────
    println!("\n--- Top 5 Designs ---\n");
    println!(
        "{:>4}  {:<40}  {:>7}  {:>8}  {:>8}  {:>8}  {:>7}  {:>6}  {:>9}",
        "#", "Candidate ID", "Score", "Sigma", "Tau(Pa)", "HI/pass", "Dose%", "Cov%", "dP(Pa)"
    );
    println!("{}", "-".repeat(120));

    for design in &top5 {
        let c = &design.candidate;
        let m = &design.metrics;

        let id_display = if c.id.len() > 40 {
            &c.id[..40]
        } else {
            &c.id
        };
        let sigma_str = if m.cavitation_number.is_finite() {
            format!("{:>8.3}", m.cavitation_number)
        } else {
            format!("{:>8}", "inf")
        };

        println!(
            "{:>4}  {:<40}  {:>7.4}  {}  {:>8.1}  {:>8.2e}  {:>6.1}%  {:>5.0}%  {:>9.0}",
            design.rank,
            id_display,
            design.score,
            sigma_str,
            m.max_main_channel_shear_pa,
            m.hemolysis_index_per_pass,
            m.cancer_dose_fraction * 100.0,
            m.well_coverage_fraction * 100.0,
            m.total_pressure_drop_pa,
        );
    }

    // ── 3. Export mesh artifacts ─────────────────────────────────────────────
    println!("\n--- Mesh Export ---\n");

    let pipeline = DesignPipeline::new(&out_dir);
    let all_artifacts: Vec<DesignArtifacts> = pipeline.export_designs(&top5)?;

    println!(
        "{:>4}  {:<40}  {:>8}  {:>8}  {:>10}  {}",
        "#", "Design ID", "Vertices", "Faces", "Watertight", "Fluid STL"
    );
    println!("{}", "-".repeat(110));

    for (i, artifacts) in all_artifacts.iter().enumerate() {
        let id_display = if artifacts.design_id.len() > 40 {
            &artifacts.design_id[..40]
        } else {
            &artifacts.design_id
        };
        println!(
            "{:>4}  {:<40}  {:>8}  {:>8}  {:>10}  {}",
            i + 1,
            id_display,
            artifacts.vertex_count,
            artifacts.face_count,
            if artifacts.watertight { "yes" } else { "NO" },
            artifacts.fluid_stl.display(),
        );
    }

    // ── 4. Detailed per-design summary ──────────────────────────────────────
    println!("\n--- Per-Design Details ---\n");

    for (design, artifacts) in top5.iter().zip(all_artifacts.iter()) {
        let m = &design.metrics;
        println!("Design #{}: {}", design.rank, design.candidate.id);
        println!("  Score:          {:.4}", design.score);
        println!("  Topology:       {:?}", design.candidate.topology);
        println!("  Cavitation no.: {}", if m.cavitation_number.is_finite() {
            format!("{:.3}", m.cavitation_number)
        } else {
            "inf (no venturi)".to_string()
        });
        println!("  Cav potential:  {:.3}", m.cavitation_potential);
        println!("  Main shear:     {:.1} Pa (FDA: {})", m.max_main_channel_shear_pa, if m.fda_main_compliant { "OK" } else { "FAIL" });
        println!("  HI/pass:        {:.2e}", m.hemolysis_index_per_pass);
        println!("  Cancer dose:    {:.1}%", m.cancer_dose_fraction * 100.0);
        println!("  Well coverage:  {:.0}%", m.well_coverage_fraction * 100.0);
        println!("  Pressure drop:  {:.0} Pa (feasible: {})", m.total_pressure_drop_pa, m.pressure_feasible);
        println!("  Residence time: {:.2} s", m.mean_residence_time_s);
        println!("  Sep3 eff:       {:.4}", m.three_pop_sep_efficiency);
        println!("  Sono proxy:     {:.3}", m.sonoluminescence_proxy);
        println!("  Mesh:           {} verts, {} faces, watertight: {}", artifacts.vertex_count, artifacts.face_count, artifacts.watertight);
        println!("  Fluid STL:      {}", artifacts.fluid_stl.display());
        if let Some(ref chip) = artifacts.chip_stl {
            println!("  Chip STL:       {}", chip.display());
        }
        println!("  OpenFOAM:       {}", artifacts.openfoam_dir.display());
        println!("  Schematic SVG:  {}", artifacts.schematic_svg.display());
        println!("  Schematic JSON: {}", artifacts.schematic_json.display());
        println!();
    }

    println!("=== Pipeline complete. Outputs in: {} ===", out_dir.display());
    Ok(())
}
