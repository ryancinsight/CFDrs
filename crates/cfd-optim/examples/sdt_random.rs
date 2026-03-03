//! Random-search design space sampler for millifluidic SDT devices.
//!
//! Generates [`N_RANDOM`] uniformly-sampled [`DesignCandidate`] objects drawn
//! from the **continuous** millifluidic design space (all 24 topology families),
//! evaluates their physics metrics, and prints:
//!
//! 1. A per-topology feasibility table — how many random designs passed hard
//!    constraints for each of the 24 topology families.
//! 2. The top-[`TOP_K`] diverse designs, using the same topology + flow + width
//!    bucket diversity filter as [`SdtOptimizer::top_k`].
//!
//! # Purpose
//!
//! - Broad initial exploration of the design space before GA seeding.
//! - Benchmarking: identifies which topologies consistently produce feasible
//!   candidates vs. which are sensitive to parameter tuning.
//! - Dataset generation: [`N_RANDOM`] evaluations in ~ a few seconds (pure Rust).
//!
//! Unlike the parametric sweep ([`SdtOptimizer`]), candidates here are sampled
//! continuously — each is an independent draw, not a point on a fixed grid.
//! This makes it far less likely to miss good regions between grid points.
//!
//! # Run with
//! ```bash
//! cargo run -p cfd-optim --example sdt_random --release
//! ```

use std::collections::HashMap;

use cfd_optim::{
    compute_metrics, sample_random_candidates, score_candidate, OptimMode, RankedDesign, SdtWeights,
};

const N_RANDOM: usize = 5_000;
const TOP_K: usize = 10;

fn main() {
    let weights = SdtWeights::default();
    let mode = OptimMode::SdtCavitation;

    println!("{}", "=".repeat(100));
    println!(
        "  cfd-optim  |  Random Design Search  |  {} candidates  |  top-{}  |  SdtCavitation",
        N_RANDOM, TOP_K
    );
    println!("{}", "=".repeat(100));
    println!("  Sampling uniformly from the continuous design space (all 24 topology families).");
    println!("  Unlike the parametric grid, every candidate is an independent random draw.");
    println!("{}", "-".repeat(100));

    // ── Generate random candidates ────────────────────────────────────────────
    let mut rng = rand::thread_rng();
    let candidates = sample_random_candidates(N_RANDOM, &mut rng);
    eprintln!(
        "[sdt_random] generated {} random candidates",
        candidates.len()
    );

    // ── Evaluate all candidates ───────────────────────────────────────────────
    // per_topology: short_code → (feasible_count, total_count, name)
    let mut per_topology: HashMap<String, (usize, usize, String)> = HashMap::new();
    let mut evaluated: Vec<RankedDesign> = Vec::with_capacity(N_RANDOM / 2);

    for candidate in candidates {
        let short = candidate.topology.short().to_string();
        let name = candidate.topology.name().to_string();
        let entry = per_topology.entry(short).or_insert((0, 0, name));
        entry.1 += 1; // total

        match compute_metrics(&candidate) {
            Ok(metrics) => {
                let score = score_candidate(&metrics, mode, &weights);
                if score > 0.0 {
                    entry.0 += 1; // feasible
                    evaluated.push(RankedDesign {
                        rank: 0,
                        candidate,
                        metrics,
                        score,
                    });
                }
            }
            Err(_) => {} // skip invalid parameter combinations
        }
    }

    let total_feasible = evaluated.len();
    eprintln!(
        "[sdt_random] {} / {} candidates passed hard constraints ({:.1}%)",
        total_feasible,
        N_RANDOM,
        total_feasible as f64 / N_RANDOM as f64 * 100.0
    );

    // ── Per-topology feasibility report ──────────────────────────────────────
    println!("\n  Per-Topology Feasibility Report");
    println!(
        "  {:<8}  {:<42}  {:>8}  {:>8}  {:>8}",
        "Short", "Topology Name", "Feasible", "Total", "Rate%"
    );
    println!("  {}", "-".repeat(79));

    let mut topo_rows: Vec<(String, String, usize, usize)> = per_topology
        .into_iter()
        .map(|(short, (feas, tot, name))| (short, name, feas, tot))
        .collect();
    // Sort by feasible count descending, then short-code alphabetically
    topo_rows.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));

    for (short, name, feas, tot) in &topo_rows {
        let rate = if *tot > 0 {
            *feas as f64 / *tot as f64 * 100.0
        } else {
            0.0
        };
        let name_trunc = if name.len() > 42 { &name[..42] } else { name };
        println!(
            "  {:<8}  {:<42}  {:>8}  {:>8}  {:>7.1}%",
            short, name_trunc, feas, tot, rate
        );
    }

    let topologies_with_feasible = topo_rows.iter().filter(|(_, _, f, _)| *f > 0).count();
    println!(
        "\n  {}/{} topology families produced at least 1 feasible candidate.",
        topologies_with_feasible,
        topo_rows.len()
    );

    // Early exit if nothing feasible
    if evaluated.is_empty() {
        println!("\n  No feasible candidates found — try increasing N_RANDOM.");
        return;
    }

    // ── Sort and diversity-aware top-K selection ──────────────────────────────
    evaluated.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            // Tiebreak: lower cavitation number (closer to inception) wins for cav mode
            .then_with(|| {
                a.metrics
                    .cavitation_number
                    .partial_cmp(&b.metrics.cavitation_number)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    let mut diverse: Vec<RankedDesign> = Vec::with_capacity(TOP_K);
    'outer: for d in evaluated {
        for prev in &diverse {
            if d.candidate.topology == prev.candidate.topology
                && (d.candidate.flow_rate_m3_s - prev.candidate.flow_rate_m3_s).abs() < 0.5e-8
                && (d.candidate.channel_width_m - prev.candidate.channel_width_m).abs() < 0.5e-3
            {
                continue 'outer; // same cluster — skip
            }
        }
        diverse.push(d);
        if diverse.len() == TOP_K {
            break;
        }
    }
    for (i, d) in diverse.iter_mut().enumerate() {
        d.rank = i + 1;
    }

    if diverse.is_empty() {
        println!("\n  Not enough diverse candidates found — try increasing N_RANDOM.");
        return;
    }

    // ── Print top-K ───────────────────────────────────────────────────────────
    println!("\n  ── Top-{TOP_K} Diverse Designs ──");
    println!(
        "\n  {:>4}  {:<42}  {:>7}  {:>9}  {:>8}  {:>6}",
        "#", "Candidate ID", "Score", "σ", "HI/pass", "Cov%"
    );
    println!("{}", "-".repeat(100));

    for d in &diverse {
        let m = &d.metrics;
        let sigma_col = if m.cavitation_number.is_finite() {
            format!("{:>9.3}", m.cavitation_number)
        } else {
            format!("{:>9}", "∞")
        };
        let id_trunc = if d.candidate.id.len() > 42 {
            &d.candidate.id[..42]
        } else {
            &d.candidate.id
        };
        println!(
            "  {:>4}  {:<42}  {:>7.4}  {}  {:>8.2e}  {:>5.0}%",
            d.rank,
            id_trunc,
            d.score,
            sigma_col,
            m.hemolysis_index_per_pass,
            m.well_coverage_fraction * 100.0,
        );
    }

    // ── Detailed winner ───────────────────────────────────────────────────────
    if let Some(winner) = diverse.first() {
        let c = &winner.candidate;
        let m = &winner.metrics;
        println!("\n  ── Winner detail ──");
        println!("    ID         : {}", c.id);
        println!("    Topology   : {}", c.topology.name());
        println!(
            "    Flow rate  : {:.2} mL/min  ({:.3e} m³/s)",
            c.flow_rate_m3_s * 6e7,
            c.flow_rate_m3_s
        );
        println!("    Inlet gauge: {:.0} kPa", c.inlet_gauge_pa * 1e-3);
        if c.topology.has_venturi() {
            println!(
                "    Throat Ø   : {:.0} µm  (L = {:.0} µm)",
                c.throat_diameter_m * 1e6,
                c.throat_length_m * 1e6
            );
            println!(
                "    σ          : {:.4}  (cavitation if < 1)",
                m.cavitation_number
            );
        }
        println!(
            "    Channel    : {:.0} × {:.0} µm",
            c.channel_width_m * 1e6,
            c.channel_height_m * 1e6
        );
        println!(
            "    Segments   : {}  ×  {:.1} mm",
            c.serpentine_segments,
            c.segment_length_m * 1e3
        );
        println!(
            "    Main shear : {:.1} Pa  (FDA ≤ 150 Pa: {})",
            m.max_main_channel_shear_pa,
            if m.fda_main_compliant { "PASS" } else { "FAIL" }
        );
        println!(
            "    PAI/pass   : {:.2e}  (limit 5e-4: {})",
            m.platelet_activation_index,
            if m.platelet_activation_index <= 5e-4 {
                "PASS"
            } else {
                "WARN"
            }
        );
        println!(
            "    ΔP total   : {:.0} Pa  ({:.2} kPa)",
            m.total_pressure_drop_pa,
            m.total_pressure_drop_pa * 1e-3
        );
        println!("    Score      : {:.4}", winner.score);
    }

    println!("\n{}", "=".repeat(100));
    println!("  Done.  Evaluated {N_RANDOM} random designs.");
    println!("{}", "=".repeat(100));
}
