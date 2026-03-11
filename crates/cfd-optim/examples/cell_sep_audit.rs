//! Cell-separation physics audit using the memory-efficient EvaluatedPool.
//!
//! Evaluates the Milestone 12 candidate space **once**, then scores under
//! all three OptimizationGoals without re-evaluating physics.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example cell_sep_audit --no-default-features --release
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, BlueprintObjectiveEvaluation, EvaluatedPool,
    OptimizationGoal,
};
use std::time::Instant;

fn print_header() {
    println!(
        "{:<50} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Candidate", "Score", "CancCtr%", "WbcCtr%", "RbcPer%", "SepEff", "ResTime"
    );
    println!("{}", "-".repeat(118));
}

fn print_row(e: &BlueprintObjectiveEvaluation) {
    println!(
        "{:<50} {:>8.5} {:>10.1} {:>10.1} {:>10.1} {:>10.4} {:>10.3}",
        e.candidate_id,
        e.score_or_zero(),
        e.separation.cancer_center_fraction * 100.0,
        e.separation.wbc_center_fraction * 100.0,
        e.separation.rbc_peripheral_fraction * 100.0,
        e.separation.separation_efficiency,
        e.residence.treatment_residence_time_s,
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t0 = Instant::now();

    // Build blueprint candidates.
    println!("\n  Building Milestone 12 candidate space...");
    let t_build = Instant::now();
    let candidates = build_milestone12_blueprint_candidate_space()?;
    let n = candidates.len();
    let d_build = t_build.elapsed();
    println!("  {} candidates built in {:.2}s", n, d_build.as_secs_f64());

    // Evaluate all candidates once.
    println!("  Evaluating physics (1-D network solves)...");
    let t_eval = Instant::now();
    let pool = EvaluatedPool::from_candidates(&candidates);
    let d_eval = t_eval.elapsed();
    let heap_mb = pool.heap_bytes() as f64 / (1024.0 * 1024.0);
    println!(
        "  Pool ready: {} / {} evaluated in {:.1}s  ({:.1} MB heap)\n",
        pool.len(),
        n,
        d_eval.as_secs_f64(),
        heap_mb,
    );

    // ── Option 1: AsymmetricSplitResidenceSeparation ────────────────────
    println!("{}", "=".repeat(118));
    println!("  OPTION 1 — AsymmetricSplitResidenceSeparation");
    println!("  Updated physics: SE_CANCER=1.85, Fahraeus margination correction");
    println!("{}", "=".repeat(118));

    let t1 = Instant::now();
    let opt1 = pool.top_k(20, OptimizationGoal::AsymmetricSplitResidenceSeparation)?;
    println!("  [{:.3}s] Top 20 designs:\n", t1.elapsed().as_secs_f64());
    print_header();
    for e in &opt1 {
        print_row(e);
    }

    // ── Option 2: AsymmetricSplitVenturiCavitationSelectivity ──────────────────────────────
    println!("\n\n{}", "=".repeat(118));
    println!("  OPTION 2 — AsymmetricSplitVenturiCavitationSelectivity");
    println!("  Updated physics: SE_CANCER=1.85, Fahraeus margination correction");
    println!("{}", "=".repeat(118));

    let t2 = Instant::now();
    let opt2 = pool.top_k(
        20,
        OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
    )?;
    println!("  [{:.3}s] Top 20 designs:\n", t2.elapsed().as_secs_f64());
    print_header();
    for e in &opt2 {
        print_row(e);
    }

    // ── Option 3: InPlaceDeanSerpentineRefinement ──────────────────────────────
    println!("\n\n{}", "=".repeat(118));
    println!("  OPTION 3 — InPlaceDeanSerpentineRefinement");
    println!("{}", "=".repeat(118));

    let t3 = Instant::now();
    match pool.top_k(20, OptimizationGoal::InPlaceDeanSerpentineRefinement) {
        Ok(opt3) => {
            println!("  [{:.3}s] Top 20 designs:\n", t3.elapsed().as_secs_f64());
            print_header();
            for e in &opt3 {
                print_row(e);
            }
        }
        Err(e) => {
            println!("  [{:.3}s] {e}", t3.elapsed().as_secs_f64());
        }
    }

    let total = t0.elapsed();
    println!(
        "\n=== Audit complete: {:.1}s total (build: {:.1}s, eval: {:.1}s, heap: {:.1} MB) ===",
        total.as_secs_f64(),
        d_build.as_secs_f64(),
        d_eval.as_secs_f64(),
        heap_mb,
    );
    Ok(())
}
