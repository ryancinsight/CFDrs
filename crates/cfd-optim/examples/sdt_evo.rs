//! Evolutionary genetic algorithm search across all 14 millifluidic topology
//! families optimised for SDT cavitation and uniform exposure.
//!
//! Demonstrates the [`GeneticOptimizer`] which searches the continuous design
//! space (bounded by the same physical ranges as the parametric sweep) using:
//!
//! - Tournament selection (k = 3)
//! - Simulated Binary Crossover (SBX, η = 2)
//! - Polynomial mutation (η_m = 20)
//!
//! The best 5 designs found across all generations are printed and exported.
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_evo
//! ```

use std::path::PathBuf;

use cfd_optim::{
    save_comparison_svg, save_schematic_svg, save_top5_json,
    evo::GeneticOptimizer,
    OptimMode, SdtWeights,
};

fn main() {
    let out_dir = PathBuf::from("outputs/cfd-optim");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let weights = SdtWeights::default();

    println!("{}", "=".repeat(100));
    println!("  cfd-optim  |  Evolutionary Optimizer  |  Single-inlet / single-outlet designs");
    println!("{}", "=".repeat(100));
    println!("  Topologies searched: all 14 families (SingleVenturi … TrifurcationSerpentine)");
    println!("  Algorithm: real-coded GA with SBX crossover + polynomial mutation");
    println!("  Genome: 8 continuous genes (topology, Q, P_gauge, d_throat, w_ch, n_segs, L_seg, R_bend)");
    println!("{}", "-".repeat(100));

    let modes: &[(&str, OptimMode)] = &[
        ("evo_cavitation",  OptimMode::SdtCavitation),
        ("evo_exposure",    OptimMode::UniformExposure),
    ];

    for (slug, mode) in modes {
        println!("\n── Evolutionary search: {} ──", slug);

        let optimizer = GeneticOptimizer::new(*mode, weights)
            .with_population(60)
            .with_max_generations(100)
            .with_top_k(5);

        match optimizer.run() {
            Ok(designs) => {
                println!(
                    "  {:>4}  {:<42}  {:>7}  {:>9}  {:>8}  {:>6}",
                    "#", "Candidate ID", "Score", "σ / Unif", "HI/pass", "Cov%"
                );
                println!("{}", "-".repeat(100));

                for d in &designs {
                    let m = &d.metrics;
                    let metric_col = match mode {
                        OptimMode::SdtCavitation =>
                            if m.cavitation_number.is_finite() {
                                format!("{:>9.3}", m.cavitation_number)
                            } else {
                                format!("{:>9}", "∞")
                            },
                        OptimMode::UniformExposure =>
                            format!("{:>9.4}", m.flow_uniformity),
                        _ => format!("{:>9}", "—"),
                    };
                    let id_trunc = if d.candidate.id.len() > 42 {
                        &d.candidate.id[..42]
                    } else {
                        &d.candidate.id
                    };
                    println!(
                        "  {:>4}  {:<42}  {:>7.4}  {}  {:>8.2e}  {:>5.0}%",
                        d.rank, id_trunc, d.score, metric_col,
                        m.hemolysis_index_per_pass,
                        m.well_coverage_fraction * 100.0,
                    );
                }

                // Detailed printout of the winner
                if let Some(winner) = designs.first() {
                    let c = &winner.candidate;
                    let m = &winner.metrics;
                    println!("\n  ── Winner detail ──");
                    println!("    ID         : {}", c.id);
                    println!("    Topology   : {}", c.topology.name());
                    println!("    Flow rate  : {:.2} mL/min  ({:.3e} m³/s)",
                        c.flow_rate_m3_s * 6e7, c.flow_rate_m3_s);
                    println!("    Inlet gauge: {:.0} kPa", c.inlet_gauge_pa * 1e-3);
                    if c.topology.has_venturi() {
                        println!("    Throat Ø   : {:.0} µm  (L = {:.0} µm)",
                            c.throat_diameter_m * 1e6, c.throat_length_m * 1e6);
                        println!("    σ          : {:.4}  (cavitation if < 1)",
                            m.cavitation_number);
                    }
                    println!("    Channel    : {:.0} × {:.0} µm",
                        c.channel_width_m * 1e6, c.channel_height_m * 1e6);
                    println!("    Segments   : {}  ×  {:.1} mm",
                        c.serpentine_segments, c.segment_length_m * 1e3);
                    println!("    Main shear : {:.1} Pa  (FDA ≤ 150 Pa: {})",
                        m.max_main_channel_shear_pa,
                        if m.fda_main_compliant { "PASS" } else { "FAIL" });
                    println!("    Score      : {:.4}", winner.score);
                }

                // Export
                let json_path = out_dir.join(format!("{slug}.json"));
                let svg_path  = out_dir.join(format!("{slug}.svg"));
                if let Err(e) = save_top5_json(&designs, &json_path) {
                    eprintln!("  JSON export failed: {e}");
                } else {
                    println!("\n  JSON → {}", json_path.display());
                }
                if let Err(e) = save_comparison_svg(&designs, &svg_path, *mode) {
                    eprintln!("  SVG export failed: {e}");
                } else {
                    println!("  SVG  → {}", svg_path.display());
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
            }
            Err(e) => eprintln!("  Evolutionary search failed: {e}"),
        }
    }

    println!("\n{}", "=".repeat(100));
    println!("  Done.");
    println!("{}", "=".repeat(100));
}
