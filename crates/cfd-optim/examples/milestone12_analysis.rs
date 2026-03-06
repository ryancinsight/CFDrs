//! Milestone 12 deep-analysis companion (7 parts).
//!
//! Performs computationally intensive analyses not required for every report
//! regeneration.  All outputs go to `outputs/milestone12_analysis/`.
//!
//! | Part | Content |
//! |------|---------|
//! | 1 | Multi-mode parametric sweep (6 modes) |
//! | 2 | Cavitation envelope: σ vs (throat, flow, pressure) — SingleVenturi grid |
//! | 3 | Joint feasibility: primitive selective venturi with σ < 1 AND cancer enrichment > 20% |
//! | 4 | NSGA-II Pareto front: max cancer_cav · min lysis_risk · max sep3 |
//! | 5 | Design robustness: ±10%/±20% sensitivity sweep (flow, pressure, throat) |
//! | 6 | Wall shear percentiles & diffuser pressure recovery (ASTM F1841-20) |
//! | 7 | Low-flow band (30–60 mL/min) primitive selective analysis with gauge compensation |
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_analysis --no-default-features
//! ```

use cfd_optim::{
    analysis::{
        compute_sdt_pareto_front, robustness_sweep, RobustnessReport, SdtParetoFront,
        STANDARD_PERTURBATIONS,
    },
    build_candidate_space, compute_metrics, save_top5_json, score_candidate, score_description,
    DesignTopology, OptimMode, RankedDesign, SdtMetrics, SdtOptimizer, SdtWeights,
};
use cfd_schematics::{write_well_plate_diagram_svg, CandidateZoneData};
use serde::Serialize;
use std::collections::HashSet;
use std::path::Path;

// ── LowFlowBandRow struct ─────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize)]
struct LowFlowBandRow {
    seed_label: String,
    seed_candidate_id: String,
    flow_ml_min: f64,
    best_gauge_kpa: Option<f64>,
    feasible: bool,
    cavitation_ready_any: bool,
    best_combined_score: f64,
    cavitation_number: Option<f64>,
    cancer_targeted_cavitation: Option<f64>,
    wbc_recovery: Option<f64>,
    rbc_venturi_exposure_fraction: Option<f64>,
    cancer_rbc_cavitation_bias_index: Option<f64>,
    cif_remerge_proximity_score: Option<f64>,
    selective_cavitation_delivery_index: Option<f64>,
    clotting_risk_index: Option<f64>,
    projected_hemolysis_15min_pediatric_pct: Option<f64>,
    total_ecv_ml: Option<f64>,
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("milestone12_analysis");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let weights = SdtWeights::default();
    let combined_mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };

    println!("=== Milestone 12 Deep Analysis (7 parts) ===\n");
    println!("Outputs: {}\n", out_dir.display());

    // ── Part 1: Multi-mode parametric sweep ───────────────────────────────────
    let modes: &[(&str, OptimMode)] = &[
        ("CombinedLeuka", combined_mode),
        ("SdtCavitation", OptimMode::SdtCavitation),
        ("ThreePopSep", OptimMode::ThreePopSeparation),
        ("SdtTherapy", OptimMode::SdtTherapy),
        ("HydroSDT", OptimMode::HydrodynamicCavitationSDT),
        ("RbcProtected", OptimMode::RbcProtectedSdt),
    ];

    let mut mode_results: Vec<(&str, OptimMode, Vec<RankedDesign>)> = Vec::new();

    for (label, mode) in modes {
        println!("{}", "=".repeat(110));
        println!(
            "  PART 1 — Parametric Sweep  |  Mode: {label}  ({})",
            score_description(*mode)
        );
        println!("{}", "=".repeat(110));
        println!("  Constraints: ΔP ≤ gauge  |  main-channel shear ≤ 150 Pa  |  plate fits");
        println!("{}", "-".repeat(110));

        match SdtOptimizer::new(*mode, weights).top_5() {
            Ok(designs) => {
                print_table_header();
                for d in &designs {
                    print_row(d);
                }
                let json = out_dir.join(format!("p1_{label}.json"));
                if let Err(e) = save_top5_json(&designs, &json) {
                    eprintln!("  WARN: JSON save failed: {e}");
                } else {
                    println!("  Saved: {}", json.display());
                }
                mode_results.push((label, *mode, designs));
            }
            Err(e) => {
                eprintln!("  ERROR ({label}): {e}");
                mode_results.push((label, *mode, Vec::new()));
            }
        }
    }

    let sdt_therapy_top5 = mode_results
        .iter()
        .find(|(l, _, _)| *l == "SdtTherapy")
        .map(|(_, _, d)| d.clone())
        .unwrap_or_default();
    let rbc_protected_top5 = mode_results
        .iter()
        .find(|(l, _, _)| *l == "RbcProtected")
        .map(|(_, _, d)| d.clone())
        .unwrap_or_default();
    let combined_top5 = mode_results
        .iter()
        .find(|(l, _, _)| *l == "CombinedLeuka")
        .map(|(_, _, d)| d.clone())
        .unwrap_or_default();

    // ── Part 2: Cavitation envelope ───────────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 2 — Cavitation Envelope  |  SingleVenturi: σ vs (throat, flow, pressure)");
    println!("{}", "=".repeat(110));
    println!("  σ < 1 = active hydrodynamic cavitation at venturi throat.");
    println!("{}", "-".repeat(110));

    let all_candidates = build_candidate_space();
    println!(
        "  Parametric space: {} candidates total",
        all_candidates.len()
    );

    let sv_candidates: Vec<_> = all_candidates
        .iter()
        .filter(|c| matches!(c.topology, DesignTopology::SingleVenturi))
        .collect();

    let mut cav_envelope: Vec<(f64, f64, f64, f64, bool)> = Vec::new();
    for c in &sv_candidates {
        if let Ok(m) = compute_metrics(c) {
            let sigma = m.cavitation_number;
            let active = sigma.is_finite() && sigma < 1.0 && m.pressure_feasible;
            cav_envelope.push((
                c.throat_diameter_m,
                c.flow_rate_m3_s,
                c.inlet_gauge_pa,
                sigma,
                active,
            ));
        }
    }

    let active_rows: Vec<_> = cav_envelope.iter().filter(|r| r.4).collect();
    println!(
        "  {:>8}  {:>10}  {:>10}  {:>8}",
        "throat µm", "flow mL/min", "gauge kPa", "σ"
    );
    println!("  {}", "-".repeat(50));
    if active_rows.is_empty() {
        println!("  (no candidates achieved σ < 1 within the current pressure envelope)");
        let min_sigma = cav_envelope
            .iter()
            .filter_map(|(_, _, _, s, _)| if s.is_finite() { Some(*s) } else { None })
            .fold(f64::INFINITY, f64::min);
        println!("  Minimum σ recorded: {min_sigma:.4}");
    } else {
        for (throat, flow, gauge, sigma, _) in &active_rows {
            println!(
                "  {:>8.0}  {:>10.2}  {:>10.0}  {:>8.4}  ← CAVITATING",
                throat * 1e6,
                flow * 6e7,
                gauge * 1e-3,
                sigma,
            );
        }
        println!(
            "\n  {} operating point(s) with σ < 1 confirmed.",
            active_rows.len()
        );
    }

    let cav_csv_path = out_dir.join("cavitation_envelope.csv");
    let mut cav_csv = String::from("throat_um,flow_ml_min,gauge_kpa,sigma,active_cavitation\n");
    for (throat, flow, gauge, sigma, active) in &cav_envelope {
        let s = if sigma.is_finite() {
            format!("{sigma:.6}")
        } else {
            "inf".to_string()
        };
        cav_csv.push_str(&format!(
            "{:.1},{:.3},{:.1},{},{}\n",
            throat * 1e6,
            flow * 6e7,
            gauge * 1e-3,
            s,
            active
        ));
    }
    if let Err(e) = std::fs::write(&cav_csv_path, &cav_csv) {
        eprintln!("  WARN: cav envelope CSV: {e}");
    } else {
        println!("  Saved: {}", cav_csv_path.display());
    }

    // ── Part 3: Joint feasibility scan ───────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 3 — Joint Feasibility: primitive selective venturi with σ < 1 AND cancer enrichment > 20%");
    println!("{}", "=".repeat(110));
    println!("  Sort by: cancer_targeted_cavitation + 0.5 × three_pop_sep_efficiency");
    println!("{}", "-".repeat(110));

    let primitive_selective_candidates: Vec<_> = all_candidates
        .iter()
        .filter(|c| matches!(c.topology, DesignTopology::PrimitiveSelectiveTree { .. }))
        .collect();

    let mut joint_feasible: Vec<(f64, &cfd_optim::DesignCandidate, cfd_optim::SdtMetrics)> =
        Vec::new();
    for c in &primitive_selective_candidates {
        if let Ok(m) = compute_metrics(c) {
            let sigma_ok = m.cavitation_number.is_finite() && m.cavitation_number < 1.0;
            if m.pressure_feasible
                && m.fda_main_compliant
                && sigma_ok
                && m.cancer_center_fraction > 0.20
            {
                let key = m.cancer_targeted_cavitation + 0.5 * m.three_pop_sep_efficiency;
                joint_feasible.push((key, c, m));
            }
        }
    }
    joint_feasible.sort_by(|a, b| b.0.total_cmp(&a.0));

    if joint_feasible.is_empty() {
        println!(
            "  (no primitive selective candidates met joint feasibility criteria in current sweep)"
        );
        let mut best_primitive: Vec<_> = primitive_selective_candidates
            .iter()
            .filter_map(|c| compute_metrics(c).ok().map(|m| (c, m)))
            .collect();
        best_primitive.sort_by(|a, b| {
            b.1.cancer_targeted_cavitation
                .total_cmp(&a.1.cancer_targeted_cavitation)
        });
        if let Some((bc, bm)) = best_primitive.first() {
            println!(
                "  Nearest: {}  σ={:.4}  cancer_ctr={:.1}%  cancer_cav={:.4}  feasible={}",
                bc.id,
                bm.cavitation_number,
                bm.cancer_center_fraction * 100.0,
                bm.cancer_targeted_cavitation,
                bm.pressure_feasible,
            );
        }
    } else {
        println!(
            "  {:>4}  {:<36}  {:>8}  {:>8}  {:>8}  {:>8}  {:>9}",
            "#", "Candidate ID", "CancCav", "Sep3", "CancCtr%", "σ", "ΔP kPa"
        );
        println!("  {}", "-".repeat(90));
        for (rank, (_, c, m)) in joint_feasible.iter().take(10).enumerate() {
            let sigma_str = if m.cavitation_number.is_finite() {
                format!("{:.4}", m.cavitation_number)
            } else {
                "∞".to_string()
            };
            println!(
                "  {:>4}  {:<36}  {:>8.4}  {:>8.4}  {:>7.1}%  {:>8}  {:>8.1}",
                rank + 1,
                truncate(&c.id, 36),
                m.cancer_targeted_cavitation,
                m.three_pop_sep_efficiency,
                m.cancer_center_fraction * 100.0,
                sigma_str,
                m.total_pressure_drop_pa * 1e-3,
            );
        }
        println!(
            "\n  {} primitive selective candidate(s) meet joint criteria.",
            joint_feasible.len()
        );
    }

    let jf_csv_path = out_dir.join("joint_feasibility.csv");
    let mut jf_csv = String::from(
        "rank,candidate_id,topology,cancer_targeted_cav,sep3,cancer_ctr_pct,sigma,delta_p_kpa\n",
    );
    for (rank, (_, c, m)) in joint_feasible.iter().take(20).enumerate() {
        let sigma_str = if m.cavitation_number.is_finite() {
            format!("{:.6}", m.cavitation_number)
        } else {
            "inf".to_string()
        };
        jf_csv.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.3},{},{:.3}\n",
            rank + 1,
            c.id,
            c.topology.short(),
            m.cancer_targeted_cavitation,
            m.three_pop_sep_efficiency,
            m.cancer_center_fraction * 100.0,
            sigma_str,
            m.total_pressure_drop_pa * 1e-3,
        ));
    }
    if let Err(e) = std::fs::write(&jf_csv_path, &jf_csv) {
        eprintln!("  WARN: joint feasibility CSV: {e}");
    } else {
        println!("  Saved: {}", jf_csv_path.display());
    }

    // ── Pre-compute RBC feasible set (Parts 4–6) ─────────────────────────────
    let rbc_all_ranked = match SdtOptimizer::new(OptimMode::RbcProtectedSdt, weights).all_ranked() {
        Ok(r) => r,
        Err(e) => {
            eprintln!("  WARN: RbcProtectedSdt all-ranked: {e}");
            Vec::new()
        }
    };
    let rbc_feasible_owned: Vec<RankedDesign> = rbc_all_ranked
        .iter()
        .filter(|d| d.metrics.pressure_feasible && d.metrics.fda_main_compliant)
        .cloned()
        .collect();
    let rbc_top10: Vec<RankedDesign> = rbc_all_ranked.iter().take(10).cloned().collect();

    // ── Part 4: Multi-Objective Pareto Front ─────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 4 — Multi-Objective Pareto Front  (3 Objectives · NSGA-II)");
    println!("{}", "=".repeat(110));
    println!("  Objectives: max(cancer_cav) · min(lysis_risk) · max(sep3_efficiency)");
    println!("  Input: pressure_feasible && fda_main_compliant candidates");
    println!("{}", "-".repeat(110));

    let pareto_front: SdtParetoFront = compute_sdt_pareto_front(&rbc_feasible_owned);
    println!("  Pareto front size: {} members", pareto_front.len());

    if !pareto_front.is_empty() {
        if let Some(d) = pareto_front.best_cancer_cav() {
            let m = &d.metrics;
            println!(
                "    Max cancer_cav  → {} | cav={:.4} lysis={:.4e} sep3={:.4}",
                truncate(&d.candidate.id, 32),
                m.cancer_targeted_cavitation,
                m.lysis_risk_index,
                m.three_pop_sep_efficiency
            );
        }
        if let Some(d) = pareto_front.safest_rbc() {
            let m = &d.metrics;
            println!(
                "    Min lysis_risk  → {} | cav={:.4} lysis={:.4e} sep3={:.4}",
                truncate(&d.candidate.id, 32),
                m.cancer_targeted_cavitation,
                m.lysis_risk_index,
                m.three_pop_sep_efficiency
            );
        }
        if let Some(d) = pareto_front.best_sep3() {
            let m = &d.metrics;
            println!(
                "    Max sep3_eff    → {} | cav={:.4} lysis={:.4e} sep3={:.4}",
                truncate(&d.candidate.id, 32),
                m.cancer_targeted_cavitation,
                m.lysis_risk_index,
                m.three_pop_sep_efficiency
            );
        }

        let top_crowded = pareto_front.top_k_by_crowding(5);
        println!(
            "\n  {:>4}  {:<36}  {:>12}  {:>10}  {:>10}",
            "#", "Candidate ID", "cancer_cav", "lysis_risk", "sep3"
        );
        println!("  {}", "-".repeat(80));
        for (rank, d) in top_crowded.iter().enumerate() {
            let m = &d.metrics;
            println!(
                "  {:>4}  {:<36}  {:>12.4}  {:>10.4e}  {:>10.4}",
                rank + 1,
                truncate(&d.candidate.id, 36),
                m.cancer_targeted_cavitation,
                m.lysis_risk_index,
                m.three_pop_sep_efficiency
            );
        }

        let pareto_json_path = out_dir.join("pareto_front.json");
        let pareto_members: Vec<RankedDesign> = pareto_front.members.clone();
        if let Err(e) = save_top5_json(&pareto_members, &pareto_json_path) {
            eprintln!("  WARN: {e}");
        } else {
            println!("\n  Saved: {}", pareto_json_path.display());
        }

        let pareto_cands: Vec<CandidateZoneData> = top_crowded
            .iter()
            .take(5)
            .map(|d| CandidateZoneData {
                label: truncate(&d.candidate.id, 8).to_string(),
                cancer_cav: d.metrics.cancer_targeted_cavitation,
                lysis_risk: d.metrics.lysis_risk_index,
                therapy_frac: d.metrics.therapy_channel_fraction,
            })
            .collect();
        let plate_svg_path = out_dir.join("pareto_treatment_zone_plate.svg");
        match write_well_plate_diagram_svg(&pareto_cands, &plate_svg_path) {
            Ok(()) => println!("  Saved: {}", plate_svg_path.display()),
            Err(e) => eprintln!("  WARN: well-plate SVG: {e}"),
        }
    }

    // ── Part 5: Design Robustness ────────────────────────────────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 5 — Design Robustness  (±10%/±20% Parameter Sensitivity)");
    println!("{}", "=".repeat(110));
    println!(
        "  Mode: RbcProtectedSdt  ·  Top-{} candidates",
        rbc_protected_top5.len().min(5)
    );
    println!("  Parameters perturbed: flow_rate, inlet_pressure, throat_diameter (if venturi)");
    println!("{}", "-".repeat(110));

    let robustness_reports: Vec<RobustnessReport> = rbc_protected_top5
        .iter()
        .map(|d| {
            robustness_sweep(
                &d.candidate,
                OptimMode::RbcProtectedSdt,
                &weights,
                &STANDARD_PERTURBATIONS,
            )
        })
        .collect();

    println!(
        "  {:>4}  {:<36}  {:>8}  {:>9}  {:>9}  {:>7}  {}",
        "#", "Candidate ID", "IsRobust", "Scr(nom)", "Scr(min)", "CV%", "WorstParam"
    );
    println!("  {}", "-".repeat(100));
    for (rank, rep) in robustness_reports.iter().enumerate() {
        println!(
            "  {:>4}  {:<36}  {:>8}  {:>9.4}  {:>9.4}  {:>6.1}%  {}",
            rank + 1,
            truncate(&rep.candidate_id, 36),
            if rep.is_robust { "YES" } else { "NO" },
            rep.score_nominal,
            rep.score_min,
            rep.score_cv * 100.0,
            rep.worst_case_param
        );
    }
    let robust_count = robustness_reports.iter().filter(|r| r.is_robust).count();
    println!(
        "\n  Robust designs (CV < 10%): {}/{}",
        robust_count,
        robustness_reports.len()
    );

    let robustness_json_path = out_dir.join("robustness_top5.json");
    let rob_json = serde_json::to_string_pretty(&robustness_reports).unwrap_or_default();
    if let Err(e) = std::fs::write(&robustness_json_path, rob_json) {
        eprintln!("  WARN: {e}");
    } else {
        println!("  Saved: {}", robustness_json_path.display());
    }

    // ── Part 6: Wall Shear Percentiles & Diffuser Recovery ───────────────────
    println!("\n{}", "=".repeat(110));
    println!(
        "  PART 6 — Wall Shear Percentile Statistics & Diffuser Pressure Recovery  |  v2.1 Physics"
    );
    println!("{}", "=".repeat(110));
    println!(
        "  ASTM F1841-20: P95 ≤ 150 Pa (FDA sustained), P99 ≤ 300 Pa (FDA transient threshold)"
    );
    println!("  Diffuser recovery: Idelchik Diagram 6-21, C_D = 0.80");
    println!("{}", "-".repeat(110));

    println!(
        "  {:>4}  {:<30}  {:>10}  {:>10}  {:>10}  {:>8}  {:>10}  {:>10}",
        "#", "Topology", "P95 [Pa]", "P99 [Pa]", "Mean [Pa]", "CV", "FDA_P%", "ΔP_diff[kPa]"
    );
    println!("  {}", "-".repeat(102));
    for (rank, d) in rbc_top10.iter().take(10).enumerate() {
        let m = &d.metrics;
        println!(
            "  {:>4}  {:<30}  {:>10.2}  {:>10.2}  {:>10.2}  {:>8.3}  {:>10}  {:>10.3}",
            rank + 1,
            truncate(d.candidate.topology.name(), 30),
            m.wall_shear_p95_pa,
            m.wall_shear_p99_pa,
            m.wall_shear_mean_pa,
            m.wall_shear_cv,
            if m.fda_shear_percentile_compliant {
                "PASS"
            } else {
                "FAIL"
            },
            m.diffuser_recovery_pa * 1e-3,
        );
    }

    let percentile_pass_count = rbc_top10
        .iter()
        .filter(|d| d.metrics.fda_shear_percentile_compliant)
        .count();
    let total_with_diffuser = rbc_top10
        .iter()
        .filter(|d| d.metrics.diffuser_recovery_pa > 0.0)
        .count();
    println!(
        "\n  FDA shear percentile compliance: {percentile_pass_count}/{} candidates PASS",
        rbc_top10.len()
    );
    println!(
        "  Candidates with diffuser recovery > 0: {total_with_diffuser}/{}",
        rbc_top10.len()
    );

    let wall_shear_json = out_dir.join("wall_shear_diffuser_top10.json");
    if let Err(e) = save_top5_json(&rbc_top10, &wall_shear_json) {
        eprintln!("  WARN: {e}");
    } else {
        println!("  Saved: {}", wall_shear_json.display());
    }
    let wall_shear_md_path = out_dir.join("wall_shear_diffuser_analysis.md");
    if let Err(e) = write_wall_shear_diffuser_markdown(
        &rbc_top10,
        percentile_pass_count,
        total_with_diffuser,
        &wall_shear_md_path,
    ) {
        eprintln!("  WARN: {e}");
    } else {
        println!("  Saved: {}", wall_shear_md_path.display());
    }

    // ── Part 7: Low-flow primitive selective band analysis ───────────────────
    println!("\n{}", "=".repeat(110));
    println!("  PART 7 — Low-Flow Band Analysis  (30–60 mL/min, gauge-compensated)");
    println!("{}", "=".repeat(110));
    println!("  Seeds: CombinedTop, SdtTherapyTop, RbcProtectedTop");
    println!("  Objective: preserve σ<1 + hard constraints + combined selectivity at lower flow");
    println!("{}", "-".repeat(110));

    let low_flow_grid_ml_min: [f64; 9] = [30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 150.0, 200.0];
    let gauge_grid_kpa: [f64; 6] = [100.0, 150.0, 200.0, 300.0, 400.0, 500.0];
    let mut low_flow_seed_ids = HashSet::new();
    let mut low_flow_rows: Vec<LowFlowBandRow> = Vec::new();

    let mut push_seed = |label: &str, ranked: Option<&RankedDesign>| {
        if let Some(d) = ranked {
            if low_flow_seed_ids.insert(d.candidate.id.clone()) {
                low_flow_rows.extend(scan_low_flow_band(
                    label,
                    d,
                    &low_flow_grid_ml_min,
                    &gauge_grid_kpa,
                    combined_mode,
                    &weights,
                ));
            }
        }
    };
    push_seed("CombinedTop", combined_top5.first());
    push_seed("SdtTherapyTop", sdt_therapy_top5.first());
    push_seed("RbcProtectedTop", rbc_protected_top5.first());

    println!(
        "  {:<16} {:>9} {:>12} {:>12} {:>12} {:>10}",
        "Seed", "Flow", "Gauge kPa", "CombScore", "SelDeliv", "ClotRisk"
    );
    println!("  {}", "-".repeat(76));
    for row in low_flow_rows.iter().filter(|r| r.flow_ml_min <= 60.0) {
        let gauge = row
            .best_gauge_kpa
            .map_or_else(|| "NA".to_string(), |v| format!("{v:.0}"));
        let scdi = row
            .selective_cavitation_delivery_index
            .map_or_else(|| "NA".to_string(), |v| format!("{v:.4}"));
        let clot = row
            .clotting_risk_index
            .map_or_else(|| "NA".to_string(), |v| format!("{v:.4}"));
        println!(
            "  {:<16} {:>7.0}mL {:>10}  {:>10.4}  {:>10}  {:>10}",
            row.seed_label, row.flow_ml_min, gauge, row.best_combined_score, scdi, clot
        );
    }

    let low_flow_json_path = out_dir.join("low_flow_band_analysis.json");
    match serde_json::to_string_pretty(&low_flow_rows) {
        Ok(payload) => {
            if let Err(e) = std::fs::write(&low_flow_json_path, payload) {
                eprintln!("  WARN: {e}");
            } else {
                println!("  Saved: {}", low_flow_json_path.display());
            }
        }
        Err(e) => eprintln!("  WARN: low-flow json: {e}"),
    }
    let low_flow_md_path = out_dir.join("low_flow_band_analysis.md");
    if let Err(e) = write_low_flow_band_markdown(&low_flow_rows, &low_flow_md_path) {
        eprintln!("  WARN: {e}");
    } else {
        println!("  Saved: {}", low_flow_md_path.display());
    }

    println!("\n{}", "=".repeat(110));
    println!("  Done.  Results in: {}", out_dir.display());
    println!("  Units: Pa, m³/s, m, s  |  Plate: ANSI/SLAS 1-2004 96-well (127.76 × 85.47 mm)");
    println!("{}", "=".repeat(110));
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn ml_min_to_m3_s(flow_ml_min: f64) -> f64 {
    flow_ml_min / 6e7
}

fn truncate(s: &str, n: usize) -> &str {
    if s.len() <= n {
        s
    } else {
        &s[..n]
    }
}

fn print_table_header() {
    println!(
        "  {:>4}  {:<38}  {:>7}  {:>8}  {:>8}  {:>7}  {:>7}  {:>7}  {:>7}",
        "#", "Candidate ID", "Score", "Sep3Eff", "HI/pass", "σ", "Unif", "Cov%", "ΔP kPa"
    );
    println!("  {}", "-".repeat(100));
}

fn print_row(d: &RankedDesign) {
    let m = &d.metrics;
    let sigma = if m.cavitation_number.is_finite() {
        format!("{:>7.3}", m.cavitation_number)
    } else {
        format!("{:>7}", "∞")
    };
    println!(
        "  {:>4}  {:<38}  {:>7.4}  {:>8.4}  {:>8.2e}  {}  {:>7.4}  {:>6.0}%  {:>7.1}",
        d.rank,
        truncate(&d.candidate.id, 38),
        d.score,
        m.three_pop_sep_efficiency,
        m.hemolysis_index_per_pass,
        sigma,
        m.flow_uniformity,
        m.well_coverage_fraction * 100.0,
        m.total_pressure_drop_pa * 1e-3,
    );
}

fn scan_low_flow_band(
    seed_label: &str,
    seed: &RankedDesign,
    flow_grid_ml_min: &[f64],
    gauge_grid_kpa: &[f64],
    score_mode: OptimMode,
    weights: &SdtWeights,
) -> Vec<LowFlowBandRow> {
    let mut rows = Vec::with_capacity(flow_grid_ml_min.len());
    for &flow_ml_min in flow_grid_ml_min {
        let mut best: Option<(f64, f64, SdtMetrics)> = None;
        let mut cavitation_ready_any = false;

        for &gauge_kpa in gauge_grid_kpa {
            let mut cand = seed.candidate.clone();
            cand.flow_rate_m3_s = ml_min_to_m3_s(flow_ml_min);
            cand.inlet_gauge_pa = gauge_kpa * 1000.0;

            let Ok(metrics) = compute_metrics(&cand) else {
                continue;
            };
            let cav_ready =
                metrics.cavitation_number.is_finite() && metrics.cavitation_number < 1.0;
            cavitation_ready_any |= cav_ready;

            let score = score_candidate(&metrics, score_mode, weights);
            if score <= 0.0 || !cav_ready {
                continue;
            }

            match &best {
                None => best = Some((score, gauge_kpa, metrics)),
                Some((best_score, best_gauge, _)) => {
                    if score > *best_score + 1e-12
                        || ((score - *best_score).abs() <= 1e-12 && gauge_kpa < *best_gauge)
                    {
                        best = Some((score, gauge_kpa, metrics));
                    }
                }
            }
        }

        if let Some((best_score, best_gauge_kpa, m)) = best {
            rows.push(LowFlowBandRow {
                seed_label: seed_label.to_string(),
                seed_candidate_id: seed.candidate.id.clone(),
                flow_ml_min,
                best_gauge_kpa: Some(best_gauge_kpa),
                feasible: true,
                cavitation_ready_any,
                best_combined_score: best_score,
                cavitation_number: Some(m.cavitation_number),
                cancer_targeted_cavitation: Some(m.cancer_targeted_cavitation),
                wbc_recovery: Some(m.wbc_recovery),
                rbc_venturi_exposure_fraction: Some(m.rbc_venturi_exposure_fraction),
                cancer_rbc_cavitation_bias_index: Some(m.cancer_rbc_cavitation_bias_index),
                cif_remerge_proximity_score: Some(m.cif_remerge_proximity_score),
                selective_cavitation_delivery_index: Some(m.selective_cavitation_delivery_index),
                clotting_risk_index: Some(m.clotting_risk_index),
                projected_hemolysis_15min_pediatric_pct: Some(
                    m.projected_hemolysis_15min_pediatric_3kg * 100.0,
                ),
                total_ecv_ml: Some(m.total_ecv_ml),
            });
        } else {
            rows.push(LowFlowBandRow {
                seed_label: seed_label.to_string(),
                seed_candidate_id: seed.candidate.id.clone(),
                flow_ml_min,
                best_gauge_kpa: None,
                feasible: false,
                cavitation_ready_any,
                best_combined_score: 0.0,
                cavitation_number: None,
                cancer_targeted_cavitation: None,
                wbc_recovery: None,
                rbc_venturi_exposure_fraction: None,
                cancer_rbc_cavitation_bias_index: None,
                cif_remerge_proximity_score: None,
                selective_cavitation_delivery_index: None,
                clotting_risk_index: None,
                projected_hemolysis_15min_pediatric_pct: None,
                total_ecv_ml: None,
            });
        }
    }
    rows
}

fn write_low_flow_band_markdown(
    rows: &[LowFlowBandRow],
    out_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut md = String::new();
    md.push_str("# Low-Flow Band Analysis (30–60 mL/min)\n\n");
    md.push_str(
        "Generated by `cargo run -p cfd-optim --example milestone12_analysis` — Part 7.\n\n",
    );
    md.push_str("| Seed | Seed Candidate | Flow mL/min | Best Gauge kPa | Feasible | σ | Combined Score | Sel. Delivery | CancCav | ClotRisk | ECV mL |\n");
    md.push_str("|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|\n");

    for row in rows {
        let gauge = row
            .best_gauge_kpa
            .map_or_else(|| "-".to_string(), |v| format!("{v:.0}"));
        let sigma = row
            .cavitation_number
            .map_or_else(|| "-".to_string(), |v| format!("{v:.4}"));
        let combined = if row.feasible {
            format!("{:.4}", row.best_combined_score)
        } else {
            "-".to_string()
        };
        let sel_delivery = row
            .selective_cavitation_delivery_index
            .map_or_else(|| "-".to_string(), |v| format!("{v:.4}"));
        let cancav = row
            .cancer_targeted_cavitation
            .map_or_else(|| "-".to_string(), |v| format!("{v:.4}"));
        let clot = row
            .clotting_risk_index
            .map_or_else(|| "-".to_string(), |v| format!("{v:.4}"));
        let ecv = row
            .total_ecv_ml
            .map_or_else(|| "-".to_string(), |v| format!("{v:.3}"));
        md.push_str(&format!(
            "| {} | `{}` | {:.0} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            row.seed_label,
            row.seed_candidate_id,
            row.flow_ml_min,
            gauge,
            if row.feasible { "PASS" } else { "FAIL" },
            sigma,
            combined,
            sel_delivery,
            cancav,
            clot,
            ecv,
        ));
    }
    std::fs::write(out_path, md)?;
    Ok(())
}

fn write_wall_shear_diffuser_markdown(
    designs: &[RankedDesign],
    percentile_pass_count: usize,
    total_with_diffuser: usize,
    out_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut md = String::new();
    md.push_str("# Wall Shear Percentile & Diffuser Recovery Analysis\n\n");
    md.push_str(
        "Generated by `cargo run -p cfd-optim --example milestone12_analysis` — Part 6.\n\n",
    );
    md.push_str("ASTM F1841-20 percentile criteria: `P95 <= 150 Pa` and `P99 <= 300 Pa`.\n\n");
    md.push_str("| Rank | Candidate ID | Topology | P95 Pa | P99 Pa | Mean Pa | CV | FDA_Pct | Diffuser kPa |\n");
    md.push_str("|---:|---|---|---:|---:|---:|---:|---|---:|\n");
    for d in designs {
        let m = &d.metrics;
        md.push_str(&format!(
            "| {} | `{}` | {} | {:.2} | {:.2} | {:.2} | {:.3} | {} | {:.3} |\n",
            d.rank,
            d.candidate.id,
            d.candidate.topology.name(),
            m.wall_shear_p95_pa,
            m.wall_shear_p99_pa,
            m.wall_shear_mean_pa,
            m.wall_shear_cv,
            if m.fda_shear_percentile_compliant {
                "PASS"
            } else {
                "FAIL"
            },
            m.diffuser_recovery_pa * 1e-3,
        ));
    }
    md.push_str(&format!(
        "\n- FDA shear percentile compliance: **{percentile_pass_count}/{}** PASS.\n",
        designs.len()
    ));
    md.push_str(&format!(
        "- Candidates with diffuser recovery > 0: **{total_with_diffuser}/{}**.\n",
        designs.len()
    ));
    std::fs::write(out_path, md)?;
    Ok(())
}
