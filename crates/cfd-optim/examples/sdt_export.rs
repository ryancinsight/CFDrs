//! Export top-5 designs to JSON and SVG for all optimisation modes.
//!
//! Runs the parametric sweep for each [`OptimMode`], exports:
//!
//! - `cfd-optim/outputs/<mode>.json` — full ranked design list with all
//!   metrics, suitable for downstream processing or archival.
//! - `cfd-optim/outputs/<mode>.svg`  — bar-chart comparison of the top 5
//!   designs with annotated key metrics.
//! - `cfd-optim/outputs/<mode>/rank<N>_<topology>.svg` — 2D channel schematic
//!   for each of the top 5 designs, showing the spatial channel layout.
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_export
//! ```

use cfd_optim::{
    evo::GeneticOptimizer, save_all_modes_json, save_comparison_svg, save_schematic_svg,
    save_top5_json, OptimMode, RankedDesign, SdtOptimizer, SdtWeights,
};

fn main() {
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("outputs");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let weights = SdtWeights::default();

    let modes: &[(&str, OptimMode, usize)] = &[
        ("cavitation", OptimMode::SdtCavitation, 5),
        ("uniform_exposure", OptimMode::UniformExposure, 5),
        ("cell_separation", OptimMode::CellSeparation, 5),
        ("three_pop", OptimMode::ThreePopSeparation, 5),
        ("sdt_therapy", OptimMode::SdtTherapy, 5),
        // Current parametric pediatric space yields 4 feasible designs under hard constraints.
        (
            "pediatric_3kg",
            OptimMode::PediatricLeukapheresis {
                patient_weight_kg: 3.0,
            },
            4,
        ),
        (
            "combined",
            OptimMode::Combined {
                cavitation_weight: 0.6,
                exposure_weight: 0.4,
            },
            5,
        ),
    ];

    // Accumulate all parametric results for combined JSON.
    let mut parametric_all: Vec<(&str, Vec<RankedDesign>)> = Vec::new();

    for (slug, mode, top_k) in modes {
        println!("\n── {slug} ──────────────────────────────────────────");
        match SdtOptimizer::new(*mode, weights).top_k(*top_k) {
            Ok(designs) => {
                let json_path = out_dir.join(format!("{slug}.json"));
                let svg_path = out_dir.join(format!("{slug}.svg"));

                match save_top5_json(&designs, &json_path) {
                    Ok(()) => println!("  JSON → {}", json_path.display()),
                    Err(e) => eprintln!("  JSON export failed: {e}"),
                }
                match save_comparison_svg(&designs, &svg_path, *mode) {
                    Ok(()) => println!("  SVG  → {}", svg_path.display()),
                    Err(e) => eprintln!("  SVG export failed: {e}"),
                }

                // Per-design 2D channel schematics
                let sch_dir = out_dir.join(slug);
                std::fs::create_dir_all(&sch_dir).expect("create schematic subdir");
                for d in &designs {
                    let sch_path = sch_dir.join(format!(
                        "rank{:02}_{}.svg",
                        d.rank,
                        d.candidate.topology.short(),
                    ));
                    match save_schematic_svg(&d.candidate, &sch_path) {
                        Ok(()) => println!("  schematic → {}", sch_path.display()),
                        Err(e) => eprintln!("  schematic error rank {}: {e}", d.rank),
                    }
                }

                // Print a brief table
                println!(
                    "  {:>4}  {:<38}  {:>7}  {:>9}  {:>8}  {:>6}",
                    "#", "Candidate ID", "Score", "σ / Sep", "HI/pass", "Cov%"
                );
                for d in &designs {
                    let m = &d.metrics;
                    let col3 = match mode {
                        OptimMode::SdtCavitation | OptimMode::Combined { .. } => {
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            }
                        }
                        OptimMode::CellSeparation => {
                            format!("{:>9.4}", m.cell_separation_efficiency)
                        }
                        OptimMode::ThreePopSeparation | OptimMode::SdtTherapy => {
                            format!("{:>9.4}", m.three_pop_sep_efficiency)
                        }
                        OptimMode::HydrodynamicCavitationSDT => {
                            format!("{:>8.1}%", m.cancer_dose_fraction * 100.0)
                        }
                        OptimMode::UniformExposure => format!("{:>9.4}", m.flow_uniformity),
                        OptimMode::PediatricLeukapheresis { .. } => {
                            format!("{:>8.1}%", m.wbc_recovery * 100.0)
                        }
                        OptimMode::CombinedSdtLeukapheresis { .. } => {
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            }
                        }
                        OptimMode::RbcProtectedSdt => {
                            format!("{:>9.4}", m.therapeutic_window_score)
                        }
                    };
                    let id_trunc = if d.candidate.id.len() > 38 {
                        &d.candidate.id[..38]
                    } else {
                        &d.candidate.id
                    };
                    println!(
                        "  {:>4}  {:<38}  {:>7.4}  {}  {:>8.2e}  {:>5.0}%",
                        d.rank,
                        id_trunc,
                        d.score,
                        col3,
                        m.hemolysis_index_per_pass,
                        m.well_coverage_fraction * 100.0,
                    );
                }

                parametric_all.push((*slug, designs));
            }
            Err(e) => eprintln!("  Optimization failed: {e}"),
        }
    }

    // Combined parametric JSON (all 5 modes in one file).
    if !parametric_all.is_empty() {
        let combined_path = out_dir.join("parametric_all.json");
        let mode_slices: Vec<(&str, &[RankedDesign])> = parametric_all
            .iter()
            .map(|(slug, designs)| (*slug, designs.as_slice()))
            .collect();
        match save_all_modes_json(&mode_slices, &combined_path) {
            Ok(()) => println!("\n  parametric_all.json → {}", combined_path.display()),
            Err(e) => eprintln!("  Combined parametric JSON failed: {e}"),
        }
    }

    // ── Evolutionary GA search (all 14 topologies, all modes) ────────────────
    println!("\n{}", "=".repeat(60));
    println!("  Evolutionary GA  |  pop=40  gen=80  mode-specific top_k");
    println!("{}", "=".repeat(60));

    // Accumulate all GA results for combined JSON.
    let mut evo_all: Vec<(String, Vec<RankedDesign>)> = Vec::new();

    for (slug, mode, top_k) in modes {
        println!("\n── evo/{slug} ──────────────────────────────────────────");
        let optimizer = GeneticOptimizer::new(*mode, weights)
            .with_population(40)
            .with_max_generations(80)
            .with_top_k(*top_k);

        match optimizer.run() {
            Ok(result) => {
                let designs = &result.top_designs;
                let evo_slug = format!("evo_{slug}");
                let json_path = out_dir.join(format!("{evo_slug}.json"));
                let svg_path = out_dir.join(format!("{evo_slug}.svg"));

                match save_top5_json(designs, &json_path) {
                    Ok(()) => println!("  JSON → {}", json_path.display()),
                    Err(e) => eprintln!("  JSON export failed: {e}"),
                }
                match save_comparison_svg(designs, &svg_path, *mode) {
                    Ok(()) => println!("  SVG  → {}", svg_path.display()),
                    Err(e) => eprintln!("  SVG export failed: {e}"),
                }

                // Per-design 2D schematics
                let sch_dir = out_dir.join(&evo_slug);
                std::fs::create_dir_all(&sch_dir).expect("create evo schematic subdir");
                for d in designs {
                    let sch_path = sch_dir.join(format!(
                        "rank{:02}_{}.svg",
                        d.rank,
                        d.candidate.topology.short(),
                    ));
                    match save_schematic_svg(&d.candidate, &sch_path) {
                        Ok(()) => println!("  schematic → {}", sch_path.display()),
                        Err(e) => eprintln!("  schematic error rank {}: {e}", d.rank),
                    }
                }

                // Brief results table (reuses the same format as the parametric sweep)
                println!(
                    "  {:>4}  {:<38}  {:>7}  {:>9}  {:>8}  {:>6}",
                    "#", "Candidate ID", "Score", "σ / Sep", "HI/pass", "Cov%"
                );
                for d in designs {
                    let m = &d.metrics;
                    let col3 = match mode {
                        OptimMode::SdtCavitation | OptimMode::Combined { .. } => {
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            }
                        }
                        OptimMode::CellSeparation => {
                            format!("{:>9.4}", m.cell_separation_efficiency)
                        }
                        OptimMode::ThreePopSeparation | OptimMode::SdtTherapy => {
                            format!("{:>9.4}", m.three_pop_sep_efficiency)
                        }
                        OptimMode::HydrodynamicCavitationSDT => {
                            format!("{:>8.1}%", m.cancer_dose_fraction * 100.0)
                        }
                        OptimMode::UniformExposure => format!("{:>9.4}", m.flow_uniformity),
                        OptimMode::PediatricLeukapheresis { .. } => {
                            format!("{:>8.1}%", m.wbc_recovery * 100.0)
                        }
                        OptimMode::CombinedSdtLeukapheresis { .. } => {
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            }
                        }
                        OptimMode::RbcProtectedSdt => {
                            format!("{:>9.4}", m.therapeutic_window_score)
                        }
                    };
                    let id_trunc = if d.candidate.id.len() > 38 {
                        &d.candidate.id[..38]
                    } else {
                        &d.candidate.id
                    };
                    println!(
                        "  {:>4}  {:<38}  {:>7.4}  {}  {:>8.2e}  {:>5.0}%",
                        d.rank,
                        id_trunc,
                        d.score,
                        col3,
                        m.hemolysis_index_per_pass,
                        m.well_coverage_fraction * 100.0,
                    );
                }

                // Collect for combined GA JSON.
                evo_all.push((evo_slug, result.top_designs.clone()));
            }
            Err(e) => eprintln!("  GA failed for {slug}: {e}"),
        }
    }

    // Combined GA JSON (all 5 modes in one file).
    if !evo_all.is_empty() {
        let evo_combined_path = out_dir.join("evo_all.json");
        let mode_slices: Vec<(&str, &[RankedDesign])> = evo_all
            .iter()
            .map(|(slug, designs)| (slug.as_str(), designs.as_slice()))
            .collect();
        match save_all_modes_json(&mode_slices, &evo_combined_path) {
            Ok(()) => println!("\n  evo_all.json → {}", evo_combined_path.display()),
            Err(e) => eprintln!("  Combined GA JSON failed: {e}"),
        }
    }

    println!("\n── Done ─────────────────────────────────────────────");
    println!("  Outputs written to: {}", out_dir.display());
}
