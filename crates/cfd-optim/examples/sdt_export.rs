//! Export top-5 designs to JSON and SVG for all optimisation modes.
//!
//! Runs the parametric sweep for each [`OptimMode`], exports:
//!
//! - `outputs/cfd-optim/<mode>.json` — full ranked design list with all
//!   metrics, suitable for downstream processing or archival.
//! - `outputs/cfd-optim/<mode>.svg`  — bar-chart comparison of the top 5
//!   designs with annotated key metrics.
//! - `outputs/cfd-optim/<mode>/rank<N>_<topology>.svg` — 2D channel schematic
//!   for each of the top 5 designs, showing the spatial channel layout.
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_export
//! ```

use std::path::PathBuf;

use cfd_optim::{
    evo::GeneticOptimizer,
    save_comparison_svg, save_schematic_svg, save_top5_json,
    OptimMode, SdtOptimizer, SdtWeights,
};

fn main() {
    let out_dir = PathBuf::from("outputs/cfd-optim");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let weights = SdtWeights::default();

    let modes: &[(&str, OptimMode)] = &[
        ("cavitation",       OptimMode::SdtCavitation),
        ("uniform_exposure", OptimMode::UniformExposure),
        ("cell_separation",  OptimMode::CellSeparation),
        ("three_pop",        OptimMode::ThreePopSeparation),
        ("combined",         OptimMode::Combined {
            cavitation_weight: 0.6,
            exposure_weight: 0.4,
        }),
    ];

    for (slug, mode) in modes {
        println!("\n── {slug} ──────────────────────────────────────────");
        match SdtOptimizer::new(*mode, weights).top_k(5) {
            Ok(designs) => {
                let json_path = out_dir.join(format!("{slug}.json"));
                let svg_path  = out_dir.join(format!("{slug}.svg"));

                match save_top5_json(&designs, &json_path) {
                    Ok(()) => println!("  JSON → {}", json_path.display()),
                    Err(e)  => eprintln!("  JSON export failed: {e}"),
                }
                match save_comparison_svg(&designs, &svg_path, *mode) {
                    Ok(()) => println!("  SVG  → {}", svg_path.display()),
                    Err(e)  => eprintln!("  SVG export failed: {e}"),
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
                        OptimMode::SdtCavitation | OptimMode::Combined { .. } =>
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            },
                        OptimMode::CellSeparation =>
                            format!("{:>9.4}", m.cell_separation_efficiency),
                        OptimMode::ThreePopSeparation =>
                            format!("{:>9.4}", m.three_pop_sep_efficiency),
                        OptimMode::UniformExposure =>
                            format!("{:>9.4}", m.flow_uniformity),
                    };
                    let id_trunc = if d.candidate.id.len() > 38 {
                        &d.candidate.id[..38]
                    } else {
                        &d.candidate.id
                    };
                    println!(
                        "  {:>4}  {:<38}  {:>7.4}  {}  {:>8.2e}  {:>5.0}%",
                        d.rank, id_trunc, d.score, col3,
                        m.hemolysis_index_per_pass,
                        m.well_coverage_fraction * 100.0,
                    );
                }
            }
            Err(e) => eprintln!("  Optimization failed: {e}"),
        }
    }

    // ── Evolutionary GA search (all 14 topologies, all modes) ────────────────
    println!("\n{}", "=".repeat(60));
    println!("  Evolutionary GA  |  pop=40  gen=80  top_k=5");
    println!("{}", "=".repeat(60));

    for (slug, mode) in modes {
        println!("\n── evo/{slug} ──────────────────────────────────────────");
        let optimizer = GeneticOptimizer::new(*mode, weights)
            .with_population(40)
            .with_max_generations(80)
            .with_top_k(5);

        match optimizer.run() {
            Ok(designs) => {
                let evo_slug  = format!("evo_{slug}");
                let json_path = out_dir.join(format!("{evo_slug}.json"));
                let svg_path  = out_dir.join(format!("{evo_slug}.svg"));

                match save_top5_json(&designs, &json_path) {
                    Ok(()) => println!("  JSON → {}", json_path.display()),
                    Err(e)  => eprintln!("  JSON export failed: {e}"),
                }
                match save_comparison_svg(&designs, &svg_path, *mode) {
                    Ok(()) => println!("  SVG  → {}", svg_path.display()),
                    Err(e)  => eprintln!("  SVG export failed: {e}"),
                }

                // Per-design 2D schematics
                let sch_dir = out_dir.join(&evo_slug);
                std::fs::create_dir_all(&sch_dir).expect("create evo schematic subdir");
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

                // Brief results table (reuses the same format as the parametric sweep)
                println!(
                    "  {:>4}  {:<38}  {:>7}  {:>9}  {:>8}  {:>6}",
                    "#", "Candidate ID", "Score", "σ / Sep", "HI/pass", "Cov%"
                );
                for d in &designs {
                    let m = &d.metrics;
                    let col3 = match mode {
                        OptimMode::SdtCavitation | OptimMode::Combined { .. } =>
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            },
                        OptimMode::CellSeparation =>
                            format!("{:>9.4}", m.cell_separation_efficiency),
                        OptimMode::ThreePopSeparation =>
                            format!("{:>9.4}", m.three_pop_sep_efficiency),
                        OptimMode::UniformExposure =>
                            format!("{:>9.4}", m.flow_uniformity),
                    };
                    let id_trunc = if d.candidate.id.len() > 38 {
                        &d.candidate.id[..38]
                    } else {
                        &d.candidate.id
                    };
                    println!(
                        "  {:>4}  {:<38}  {:>7.4}  {}  {:>8.2e}  {:>5.0}%",
                        d.rank, id_trunc, d.score, col3,
                        m.hemolysis_index_per_pass,
                        m.well_coverage_fraction * 100.0,
                    );
                }
            }
            Err(e) => eprintln!("  GA failed for {slug}: {e}"),
        }
    }

    println!("\n── Done ─────────────────────────────────────────────");
    println!("  Outputs written to: {}", out_dir.display());
}
