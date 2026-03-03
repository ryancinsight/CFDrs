//! Paediatric leukapheresis: end-to-end pipeline example.
//!
//! Demonstrates three progressive analyses on millifluidic leukapheresis chips:
//!
//! ## Phase 1 — Literature benchmark replication
//! Runs two reference designs through `compute_metrics()` and compares model
//! predictions against reported WBC recovery and RBC removal:
//!
//! - **Nivedita et al. (2017)**: Inertial spiral, 400 × 80 µm, 14 turns, 1.8 mL/min, HCT 2 %
//! - **Wu Z. et al. (2019)**: Constriction-expansion array, 300 × 60 µm, 20 cycles, 150 µL/min
//!
//! ## Phase 2 — Clinical scaling
//! Determines how many parallel Nivedita chips achieve ≥ 10 mL/min while
//! keeping ECV ≤ 25.5 mL (10 % of a 3 kg neonate's 255 mL blood volume).
//!
//! ## Phase 3 — Evolutionary optimisation
//! Genetic algorithm (population 80, 150 generations) over 19 topology families.
//! Objective: `PediatricLeukapheresis` (40 % WBC recovery + 30 % RBC removal +
//! 20 % purity + 10 % throughput feasibility).
//!
//! # Run
//! ```bash
//! cargo run --example pediatric_leukapheresis --release
//! ```
//!
//! # Output
//! Files are written to `examples/pediatric_leukapheresis/outputs/`.

mod chip_model;
mod metrics_report;
mod paper_benchmark;
mod pediatric_blood;

use cfd_optim::{compute_metrics, evo::GeneticOptimizer, save_comparison_svg, OptimMode, SdtWeights};

use pediatric_blood::{
    CLINICAL_FLOW_ML_MIN,
    DILUTED_BLOOD_VISCOSITY_PAS,
    DILUTED_HCT,
    NEONATAL_RBC_DIAMETER_M,
    NEONATE_3KG_ECV_BUDGET_ML,
    NEONATE_3KG_WEIGHT_KG,
    NIVEDITA_RBC_REMOVAL,
    NIVEDITA_WBC_RECOVERY,
    TBV_PER_KG_ML,
    MAX_ECV_FRACTION,
    WHOLE_BLOOD_HCT,
    WBC_DIAMETER_M,
    WU_WBC_PURITY,
    WU_WBC_RECOVERY,
};

fn main() {
    // Output directory: examples/pediatric_leukapheresis/outputs/
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("pediatric_leukapheresis")
        .join("outputs");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let patient_weight_kg = NEONATE_3KG_WEIGHT_KG;
    let mode              = OptimMode::PediatricLeukapheresis { patient_weight_kg };
    let weights           = SdtWeights::default();

    // ═════════════════════════════════════════════════════════════════════════
    println!("\n{}", "=".repeat(100));
    println!("  cfd-optim  |  Paediatric Leukapheresis Pipeline");
    println!(
        "  Patient: {patient_weight_kg:.0} kg neonate  |  \
         ECV budget: {:.1} mL  |  Target: {CLINICAL_FLOW_ML_MIN:.0} mL/min",
        NEONATE_3KG_ECV_BUDGET_ML,
    );
    println!("  Objective: 40% WBC recovery + 30% RBC removal + 20% WBC purity + 10% throughput");
    println!("{}", "=".repeat(100));

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 1 — Literature benchmark replication
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n{}", "-".repeat(100));
    println!("  PHASE 1 — Literature benchmark replication");
    println!("  Compare model predictions vs. reported WBC recovery / RBC removal");
    println!("{}", "-".repeat(100));

    let nivedita = paper_benchmark::nivedita_spiral();
    let wu       = paper_benchmark::wu_constriction();

    let nivedita_metrics = match compute_metrics(&nivedita) {
        Ok(m)  => m,
        Err(e) => { eprintln!("  ERROR (Nivedita metrics): {e}"); return; }
    };
    let wu_metrics = match compute_metrics(&wu) {
        Ok(m)  => m,
        Err(e) => { eprintln!("  ERROR (Wu metrics): {e}"); return; }
    };

    println!();
    metrics_report::print_benchmark_header();
    metrics_report::print_benchmark_comparison(
        &nivedita.id,
        nivedita_metrics.wbc_recovery,
        1.0 - nivedita_metrics.rbc_pass_fraction,
        NIVEDITA_WBC_RECOVERY,
        NIVEDITA_RBC_REMOVAL,
    );
    // Wu reports WBC purity (not RBC removal); use it as proxy for column alignment.
    metrics_report::print_benchmark_comparison(
        &wu.id,
        wu_metrics.wbc_recovery,
        1.0 - wu_metrics.rbc_pass_fraction,
        WU_WBC_RECOVERY,
        WU_WBC_PURITY,
    );

    println!();
    println!("  Notes:");
    let dh_nivedita = 2.0 * paper_benchmark::NIVEDITA_WIDTH_M * paper_benchmark::NIVEDITA_HEIGHT_M
        / (paper_benchmark::NIVEDITA_WIDTH_M + paper_benchmark::NIVEDITA_HEIGHT_M);
    println!(
        "   - Nivedita spiral: D_h = {:.0} um; \
         kappa_WBC = a_WBC/D_h ~{:.3} (marginal if kappa < 0.15)",
        dh_nivedita * 1e6,
        WBC_DIAMETER_M / dh_nivedita,
    );
    println!(
        "   - Model blood assumptions: diluted HCT {:.0}% (vs whole-blood {:.0}%), \
         mu_diluted ≈ {:.2} mPa.s",
        DILUTED_HCT * 100.0,
        WHOLE_BLOOD_HCT * 100.0,
        DILUTED_BLOOD_VISCOSITY_PAS * 1e3,
    );
    println!(
        "   - Cell sizes: WBC {:.1} um, neonatal RBC {:.1} um",
        WBC_DIAMETER_M * 1e6,
        NEONATAL_RBC_DIAMETER_M * 1e6,
    );
    println!(
        "   - CFL correction is weak at HCT = {:.0}% — inertial term dominates",
        paper_benchmark::NIVEDITA_HCT * 100.0,
    );
    println!("   - 1-D lumped model is conservative vs. device-specific 3-D optimisation");
    println!(
        "   - Wu design: HCT = {:.1}% → CFL enhancement Gamma ~1.0 (near-zero packing)",
        paper_benchmark::WU_HCT * 100.0,
    );

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2 — Clinical scaling
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n{}", "-".repeat(100));
    println!("  PHASE 2 — Clinical scaling to {CLINICAL_FLOW_ML_MIN:.0} mL/min");
    println!("  Scale the Nivedita spiral by parallel chip replication");
    println!("{}", "-".repeat(100));
    println!();

    let scaling = chip_model::scale_to_clinical(
        paper_benchmark::NIVEDITA_FLOW_M3S,
        nivedita_metrics.total_ecv_ml,   // single-chip ECV from Phase 1
        NEONATE_3KG_ECV_BUDGET_ML,
        CLINICAL_FLOW_ML_MIN,
    );
    scaling.print(CLINICAL_FLOW_ML_MIN);
    let parallel_nivedita = chip_model::make_parallel_candidate(&nivedita, scaling.n_parallel);
    println!("  Parallel assembly candidate ID: {}", parallel_nivedita.id);

    // Per-chip performance is unchanged by parallelisation — report from Phase 1.
    println!();
    println!("  Per-chip separation performance (each of {} replicas):", scaling.n_parallel);
    println!("    WBC recovery : {:.1}%", nivedita_metrics.wbc_recovery * 100.0);
    println!("    RBC removal  : {:.1}%", (1.0 - nivedita_metrics.rbc_pass_fraction) * 100.0);
    println!("    WBC purity   : {:.1}%", nivedita_metrics.wbc_purity * 100.0);
    println!("    HI/pass      : {:.2e}", nivedita_metrics.hemolysis_index_per_pass);

    let assembly_feasible = nivedita_metrics.wbc_recovery >= 0.70
        && scaling.ecv_ok
        && scaling.total_flow_ml_min >= CLINICAL_FLOW_ML_MIN;
    println!();
    println!(
        "  Clinical feasibility: {}",
        if assembly_feasible { "FEASIBLE" } else { "NEEDS FURTHER OPTIMISATION" }
    );

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 3 — Evolutionary optimisation
    // ─────────────────────────────────────────────────────────────────────────
    println!("\n{}", "-".repeat(100));
    println!("  PHASE 3 — Evolutionary optimisation  |  Mode: PediatricLeukapheresis");
    println!("  Patient: {patient_weight_kg:.0} kg  |  Objective: 40% WBC + 30% RBC-rem + 20% purity + 10% Q");
    println!("  Algorithm: real-coded GA  |  Population: 80  |  Generations: 150  |  Top-k: 5");
    println!("{}", "-".repeat(100));
    println!();

    let ga_optimizer = GeneticOptimizer::new(mode, weights)
        .with_population(80)
        .with_max_generations(150)
        .with_top_k(5);

    let ga_result = match ga_optimizer.run() {
        Ok(d)  => d,
        Err(e) => { eprintln!("  ERROR (GA): {e}"); return; }
    };
    let ga_top5 = &ga_result.top_designs;

    metrics_report::print_leuka_table_header();
    for d in ga_top5 {
        metrics_report::print_leuka_row(d);
    }

    println!("\n  -- Detailed breakdown (GA top-5) --");
    for d in ga_top5 {
        println!();
        print_leuka_detailed(d, patient_weight_kg);
    }

    // Save GA outputs.
    metrics_report::save_leuka_results(ga_top5, "optimized", &out_dir, patient_weight_kg);

    // Head-to-head top-2 comparison SVG.
    if ga_top5.len() >= 2 {
        let top2 = &ga_top5[..2];
        let cmp_path = out_dir.join("top2_comparison.svg");
        match save_comparison_svg(top2, &cmp_path, mode) {
            Ok(_)  => println!("  Saved: {}", cmp_path.display()),
            Err(e) => eprintln!("  WARN (top2 SVG): {e}"),
        }
    }

    // ═════════════════════════════════════════════════════════════════════════
    println!("\n{}", "=".repeat(100));
    println!("  Done.  Results in: {}", out_dir.display());
    println!("  Plate: ANSI/SLAS 1-2004 96-well  |  Blood: Casson model  |  FDA <= 150 Pa");
    println!("{}", "=".repeat(100));
}

// ── Detailed design printer ────────────────────────────────────────────────────

fn print_leuka_detailed(d: &cfd_optim::RankedDesign, patient_weight_kg: f64) {
    let c             = &d.candidate;
    let m             = &d.metrics;
    let bv_ml         = patient_weight_kg * TBV_PER_KG_ML;
    let ecv_budget_ml = bv_ml * MAX_ECV_FRACTION;

    println!("  Rank #{}: {}", d.rank, c.id);
    println!("    Topology   : {}", c.topology.name());
    println!(
        "    Flow rate  : {:.2} mL/min  ({:.3e} m3/s)",
        c.flow_rate_m3_s * 6.0e7,
        c.flow_rate_m3_s,
    );
    println!(
        "    Channel    : {:.0} um wide x {:.0} um tall",
        c.channel_width_m * 1e6,
        c.channel_height_m * 1e6,
    );
    println!("    HCT feed   : {:.1}%", c.feed_hematocrit * 100.0);
    println!(
        "    WBC recov  : {:.1}%  |  RBC removal: {:.1}%  |  Purity: {:.1}%",
        m.wbc_recovery * 100.0,
        (1.0 - m.rbc_pass_fraction) * 100.0,
        m.wbc_purity * 100.0,
    );
    println!(
        "    ECV        : {:.2} mL  (budget {:.1} mL  ->  {})",
        m.total_ecv_ml,
        ecv_budget_ml,
        if m.total_ecv_ml <= ecv_budget_ml { "PASS" } else { "FAIL" },
    );
    println!(
        "    HI/pass    : {:.3e}  |  FDA shear: {:.1} Pa ({})",
        m.hemolysis_index_per_pass,
        m.max_main_channel_shear_pa,
        if m.fda_main_compliant { "PASS" } else { "FAIL" },
    );
    println!("    Score      : {:.4}", d.score);
}
