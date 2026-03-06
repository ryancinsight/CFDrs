//! Milestone 12 unified report generator.
//!
//! Consolidates concept selection, GA optimisation, multi-fidelity validation,
//! and figure generation in one runnable example.  All outputs go to
//! `report/milestone12/` and `report/figures/`.
//!
//! # Parts
//! 1. **Two-concept selection** — parametric sweep over shared selective
//!    center-enrichment layouts, ranked by SelectiveAcousticTherapy
//!    (Option 1) and CombinedSdtLeukapheresis (Option 2), top-5 each.
//! 2. **GA warm-start search** — top-100 parametric seeds → HydrodynamicCavitationSDT GA
//!    (pop=80, gen=120) for deeper optimisation of venturi-cavitation designs.
//! 3. **Multi-fidelity validation** — 2D FVM + 3D FEM pressure-drop confirmation on
//!    the selected combined-score Option 2 design plus the HydroSDT GA design.
//! 4. **Figure generation** — schematics selected dynamically from optimizer output
//!    (Figure 4 = best selective acoustic design; Figure 5 = best combined-score
//!    venturi design; Figure 6 = best HydroSDT GA design).
//! 5. **Canonical report markdown** — `report/milestone12_results.md`.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_report --no-default-features
//! # optional fast iteration mode
//! $env:M12_FAST=1; cargo run -p cfd-optim --example milestone12_report --no-default-features
//! ```

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::venturi_flow::{VenturiGeometry, VenturiSolver2D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::fluid::CassonBlood;
use cfd_mesh::VenturiMeshBuilder;
use cfd_optim::{
    analysis::{robustness_sweep, RobustnessReport, STANDARD_PERTURBATIONS},
    build_candidate_space, compute_metrics,
    constraints::{
        BLOOD_DENSITY_KG_M3 as RHO, BLOOD_VAPOR_PRESSURE_PA as P_VAPOR_PA, M12_GA_HYDRO_SEED,
        P_ATM_PA,
    },
    metrics_cache::MetricsCache,
    reporting::{
        pct_diff, shortlist_report, write_milestone12_narrative_report, write_milestone12_results,
        Milestone12NarrativeInput, ValidationRow,
    },
    save_schematic_svg, save_top5_json, score_candidate, DesignCandidate, DesignTopology,
    GeneticOptimizer, OptimMode, PrimitiveSplitSequence, RankedDesign, SdtMetrics, SdtWeights,
};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

// ── Output structures ─────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
struct ConceptTopline {
    concept: String,
    mode: String,
    candidate_id: String,
    score: f64,
    topology: String,
    treatment_zone_mode: String,
    active_venturi_throat_count: usize,
    cavitation_number: f64,
    cancer_center_fraction: f64,
    therapeutic_window_score: f64,
    hemolysis_index_per_pass: f64,
    wall_shear_p95_pa: f64,
    mean_residence_time_s: f64,
    total_ecv_ml: f64,
    therapy_channel_fraction: f64,
}

#[derive(Debug, Serialize)]
struct Venturi2DResult {
    u_inlet_m_s: f64,
    u_throat_m_s: f64,
    dp_throat_pa: f64,
    dp_recovery_pa: f64,
    sigma_2d: f64,
    dp_bernoulli_1d_pa: f64,
}

#[derive(Debug, Serialize)]
struct Venturi3DResult {
    u_inlet_m_s: f64,
    u_throat_m_s: f64,
    dp_throat_pa: f64,
    dp_recovery_pa: f64,
    mass_error: f64,
    resolution: (usize, usize),
}

struct EvaluatedCandidate {
    candidate: DesignCandidate,
    metrics: SdtMetrics,
    ga_score: f64,
    option1_score: Option<f64>,
    combined_score: Option<f64>,
}

struct RunConfig {
    fast_mode: bool,
    fast_stride: usize,
    fast_max_candidates: usize,
    fast_nonventuri_reserve: usize,
    fast_eval_min: usize,
    fast_eval_max: usize,
    ga_seed_take: usize,
    ga_population: usize,
    ga_generations: usize,
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn report_eligible_nonventuri(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.therapy_channel_fraction > 0.0
}

fn is_selective_report_topology(topology: DesignTopology) -> bool {
    matches!(topology, DesignTopology::PrimitiveSelectiveTree { .. })
}

fn report_eligible_venturi_oncology(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.therapy_channel_fraction > 0.0
        && metrics.cavitation_number < 1.0
}

fn option2_mode() -> OptimMode {
    // Adult oncology patient (70 kg) — CTC treatment via hydrodynamic cavitation.
    // Projected 15-min hemolysis scales with blood volume (≈ 5 600 mL for an
    // adult vs 240 mL for a 3 kg neonate), producing a physically meaningful
    // hi15_gate that doesn't collapse the score.
    OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 70.0,
    }
}

fn selective_serpentine_variant(candidate: &DesignCandidate) -> Option<DesignCandidate> {
    let _ = candidate;
    None
}

fn selective_extra_throat_variant(candidate: &DesignCandidate) -> Option<DesignCandidate> {
    let DesignTopology::PrimitiveSelectiveTree { .. } = candidate.topology else {
        return None;
    };
    if !candidate.uses_venturi_treatment() || candidate.centerline_venturi_throat_count >= 4 {
        return None;
    }
    let next_throat_count = candidate.centerline_venturi_throat_count + 1;
    let mut extra = candidate.clone();
    extra.id = format!("{}-VT{}", candidate.id, next_throat_count);
    extra.centerline_venturi_throat_count = next_throat_count;
    Some(extra)
}

fn expand_m12_candidates(raw_candidates: Vec<DesignCandidate>) -> Vec<DesignCandidate> {
    let mut expanded = Vec::with_capacity(raw_candidates.len() * 5 / 4);
    for candidate in raw_candidates {
        if let Some(serpentine) = selective_serpentine_variant(&candidate) {
            expanded.push(serpentine);
        }
        if let Some(extra_throat) = selective_extra_throat_variant(&candidate) {
            expanded.push(extra_throat);
        }
        expanded.push(candidate);
    }
    expanded
}

fn run_config_from_env() -> RunConfig {
    let fast_mode = std::env::var("M12_FAST")
        .ok()
        .map(|v| {
            let s = v.trim().to_ascii_lowercase();
            s == "1" || s == "true" || s == "yes"
        })
        .unwrap_or(false);
    if fast_mode {
        let fast_stride = std::env::var("M12_FAST_STRIDE")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(3)
            .max(1);
        let fast_max_candidates = std::env::var("M12_FAST_MAX_CANDIDATES")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(8_000)
            .max(200);
        let fast_nonventuri_reserve = std::env::var("M12_FAST_NONVENTURI_RESERVE")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(1_500)
            .max(100);
        let fast_eval_min = std::env::var("M12_FAST_EVAL_MIN")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(120)
            .max(40);
        let fast_eval_max = std::env::var("M12_FAST_EVAL_MAX")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(600)
            .max(fast_eval_min);
        RunConfig {
            fast_mode: true,
            fast_stride,
            fast_max_candidates,
            fast_nonventuri_reserve,
            fast_eval_min,
            fast_eval_max,
            ga_seed_take: 50,
            ga_population: 16,
            ga_generations: 12,
        }
    } else {
        RunConfig {
            fast_mode: false,
            fast_stride: 1,
            fast_max_candidates: usize::MAX,
            fast_nonventuri_reserve: 0,
            fast_eval_min: usize::MAX,
            fast_eval_max: usize::MAX,
            ga_seed_take: 150,
            ga_population: 80,
            ga_generations: 120,
        }
    }
}

fn fast_mode_candidate_filter(c: &DesignCandidate) -> bool {
    let flow_ml_min = c.flow_rate_m3_s * 6.0e7;
    let gauge_kpa = c.inlet_gauge_pa * 1.0e-3;
    if !c.topology.has_venturi() {
        return flow_ml_min <= 250.0
            && gauge_kpa <= 300.0
            && c.channel_width_m >= 1.0e-3
            && c.channel_width_m <= 6.0e-3
            && c.channel_height_m >= 0.5e-3
            && c.channel_height_m <= 2.0e-3;
    }

    let tl_factor = if c.throat_diameter_m > 0.0 {
        c.throat_length_m / c.throat_diameter_m
    } else {
        0.0
    };

    match c.topology {
        DesignTopology::PrimitiveSelectiveTree { sequence } => {
            matches!(
                sequence,
                PrimitiveSplitSequence::Tri
                    | PrimitiveSplitSequence::TriBi
                    | PrimitiveSplitSequence::TriTri
                    | PrimitiveSplitSequence::TriBiBi
                    | PrimitiveSplitSequence::TriBiTri
                    | PrimitiveSplitSequence::TriTriBi
                    | PrimitiveSplitSequence::TriTriTri
            ) && (1..=4).contains(&c.centerline_venturi_throat_count)
                && flow_ml_min >= 60.0
                && flow_ml_min <= 550.0
                && gauge_kpa >= 25.0
                && gauge_kpa <= 550.0
                && c.throat_diameter_m >= 25.0e-6
                && c.throat_diameter_m <= 120.0e-6
                && tl_factor <= 5.0
                && c.channel_width_m >= 3.0e-3
                && c.channel_width_m <= 8.0e-3
                && c.channel_height_m >= 0.5e-3
                && c.channel_height_m <= 2.5e-3
                && c.cif_pretri_center_frac >= 0.45
                && c.cif_pretri_center_frac <= 0.62
                && c.cif_terminal_tri_center_frac >= 0.33
                && c.cif_terminal_tri_center_frac <= 0.58
        }
        _ => false,
    }
}

fn fast_eval_priority(c: &DesignCandidate) -> (u8, i64, i64, i64) {
    let family = match c.topology {
        DesignTopology::PrimitiveSelectiveTree {
            sequence: PrimitiveSplitSequence::TriTri,
        } => 0,
        DesignTopology::PrimitiveSelectiveTree { .. } => 1,
        _ => 1,
    };
    // Dual-centre priority: candidates close to EITHER the acoustic sweet spot
    // (q=120, g=200, d=55) or the cavitation sweet spot (q=300, g=400, d=35)
    // rank highly, ensuring both Option 1 and Option 2 are well-sampled.
    let q_ml = c.flow_rate_m3_s * 6.0e7;
    let g_kpa = c.inlet_gauge_pa * 1.0e-3;
    let d_um = c.throat_diameter_m * 1.0e6;
    let acoustic_pen = (q_ml - 120.0).abs() + (g_kpa - 200.0).abs() + (d_um - 55.0).abs();
    let cav_pen = (q_ml - 300.0).abs() + (g_kpa - 400.0).abs() + (d_um - 35.0).abs();
    let best_pen = acoustic_pen.min(cav_pen).round() as i64;
    (family, best_pen, 0, 0)
}

/// Priority for acoustic-only candidates: centre on q = 120 mL/min, g = 200 kPa,
/// d = 55 μm — the sweet spot for selective routing without cavitation.
fn fast_eval_priority_acoustic(c: &DesignCandidate) -> (u8, i64, i64, i64) {
    let family = match c.topology {
        DesignTopology::PrimitiveSelectiveTree {
            sequence: PrimitiveSplitSequence::TriTri,
        } => 0,
        DesignTopology::PrimitiveSelectiveTree { .. } => 1,
        _ => 1,
    };
    let q_pen = (c.flow_rate_m3_s * 6.0e7 - 120.0).abs().round() as i64;
    let g_pen = (c.inlet_gauge_pa * 1.0e-3 - 200.0).abs().round() as i64;
    let d_pen = (c.throat_diameter_m * 1.0e6 - 55.0).abs().round() as i64;
    (family, q_pen, g_pen, d_pen)
}

fn validate_venturi_candidate(
    track: &str,
    d: &RankedDesign,
    out_dir: &Path,
) -> Result<ValidationRow, Box<dyn std::error::Error>> {
    let c = &d.candidate;
    let m = &d.metrics;
    let id = &c.id;

    if !c.topology.has_venturi() {
        return Err(format!("validation requires venturi topology: {id}").into());
    }

    let q = c.per_venturi_flow();
    let a_in = c.inlet_area_m2();
    let a_th = c.throat_area_m2();
    let v_in_1d = q / a_in.max(1e-30);
    let v_th_1d = q / a_th.max(1e-30);
    let dp_bernoulli = 0.5 * RHO * (v_th_1d * v_th_1d - v_in_1d * v_in_1d);

    // 2D FVM
    let geom2d = VenturiGeometry::<f64>::new(
        c.inlet_diameter_m,
        c.throat_diameter_m,
        3e-3,
        2e-3,
        c.throat_length_m,
        4e-3,
        c.channel_height_m,
    );
    let blood_2d = BloodModel::Casson(CassonBlood::<f64>::normal_blood());
    let cr = c.inlet_diameter_m / c.throat_diameter_m.max(1e-12);
    let ny_2d = (4.0 * cr).round().clamp(40.0, 200.0) as usize;
    let beta_2d = (1.0 - 4.0 * c.throat_diameter_m / c.inlet_diameter_m.max(1e-12)).clamp(0.0, 0.9);
    let mut solver2d = VenturiSolver2D::new_stretched(geom2d, blood_2d, RHO, 60, ny_2d, beta_2d);
    let area_2d = c.inlet_diameter_m * c.channel_height_m;
    let u_inlet_2d = q / area_2d.max(1e-30);

    let sol2d = solver2d
        .solve(u_inlet_2d)
        .map_err(|e| format!("2D FVM failed for {id}: {e}"))?;
    let p_abs_inlet = P_ATM_PA + c.inlet_gauge_pa;
    let p_abs_throat = p_abs_inlet + sol2d.dp_throat;
    let dyn_p = 0.5 * RHO * sol2d.u_throat * sol2d.u_throat;
    let sigma_2d = if dyn_p > 1e-12 {
        (p_abs_throat - P_VAPOR_PA) / dyn_p
    } else {
        f64::INFINITY
    };
    let dp_2d = -sol2d.dp_throat;
    if !dp_2d.is_finite() || !sigma_2d.is_finite() {
        return Err(format!("non-finite 2D validation values for {id}").into());
    }
    let r2d = Venturi2DResult {
        u_inlet_m_s: u_inlet_2d,
        u_throat_m_s: sol2d.u_throat,
        dp_throat_pa: dp_2d,
        dp_recovery_pa: -sol2d.dp_recovery,
        sigma_2d,
        dp_bernoulli_1d_pa: dp_bernoulli,
    };
    if let Ok(json) = serde_json::to_string_pretty(&r2d) {
        std::fs::write(out_dir.join(format!("{id}_2d_venturi.json")), json)?;
    }

    // 3D FEM (deterministic resolution).
    // Resolution (60,10) balances fidelity against solve time for high-CR throats.
    let res3d = (60_usize, 10_usize);
    let builder3d = VenturiMeshBuilder::<f64>::new(
        c.inlet_diameter_m,
        c.throat_diameter_m,
        5.0 * c.inlet_diameter_m,
        3.0 * c.inlet_diameter_m,
        c.throat_length_m,
        7.0 * c.inlet_diameter_m,
        5.0 * c.inlet_diameter_m,
    )
    .with_resolution(res3d.0, res3d.1)
    .with_circular(false);

    let config3d = VenturiConfig3D::<f64> {
        inlet_flow_rate: q,
        resolution: res3d,
        circular: false,
        rect_height: Some(c.channel_height_m),
        ..Default::default()
    };

    let sol3d = VenturiSolver3D::new(builder3d, config3d)
        .solve(CarreauYasudaBlood::<f64>::normal_blood())
        .map_err(|e| format!("3D FEM failed for {id}: {e}"))?;
    let dp_3d = sol3d.dp_throat.abs();
    let mass_err_3d = sol3d.mass_error.abs();
    if !dp_3d.is_finite() || !mass_err_3d.is_finite() {
        return Err(format!("non-finite 3D validation values for {id}").into());
    }
    let r3d = Venturi3DResult {
        u_inlet_m_s: sol3d.u_inlet,
        u_throat_m_s: sol3d.u_throat,
        dp_throat_pa: dp_3d,
        dp_recovery_pa: sol3d.dp_recovery.abs(),
        mass_error: mass_err_3d,
        resolution: res3d,
    };
    if let Ok(json) = serde_json::to_string_pretty(&r3d) {
        std::fs::write(out_dir.join(format!("{id}_3d_result.json")), json)?;
    }

    let row = ValidationRow {
        track: track.to_string(),
        id: id.clone(),
        topology: c.topology.short().to_string(),
        dp_1d_bernoulli_pa: dp_bernoulli,
        dp_2d_fvm_pa: dp_2d,
        dp_3d_fem_pa: dp_3d,
        agreement_1d_2d_pct: pct_diff(dp_bernoulli, dp_2d),
        agreement_2d_3d_pct: pct_diff(dp_2d, dp_3d),
        mass_error_3d_pct: mass_err_3d * 100.0,
        sigma_1d: m.cavitation_number,
        sigma_2d,
        score: d.score,
    };
    if !row.agreement_1d_2d_pct.is_finite()
        || !row.agreement_2d_3d_pct.is_finite()
        || !row.mass_error_3d_pct.is_finite()
    {
        return Err(format!("non-finite validation agreement values for {id}").into());
    }
    if let Ok(json) = serde_json::to_string_pretty(&row) {
        std::fs::write(out_dir.join(format!("{id}_validation.json")), json)?;
    }

    Ok(row)
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("cfd-optim crate has a parent")
        .parent()
        .expect("crates/ has a workspace root")
        .to_path_buf();

    let out_dir = workspace_root.join("report").join("milestone12");
    let figures_dir = workspace_root.join("report").join("figures");
    std::fs::create_dir_all(&out_dir)?;
    std::fs::create_dir_all(&figures_dir)?;

    println!("=== Milestone 12 Unified Report Generator ===\n");
    println!("Outputs: {}", out_dir.display());

    let weights = SdtWeights::default();
    let option1_mode = OptimMode::SelectiveAcousticTherapy;
    let oncology_mode = option2_mode();
    let run_cfg = run_config_from_env();
    let raw_candidates = expand_m12_candidates(build_candidate_space());
    let all_candidates: Vec<DesignCandidate> = if run_cfg.fast_mode {
        let mut acoustic_reserve: Vec<DesignCandidate> = Vec::new();
        let mut venturi_filtered: Vec<DesignCandidate> = Vec::new();
        for candidate in raw_candidates {
            if !candidate.uses_venturi_treatment()
                && is_selective_report_topology(candidate.topology)
            {
                acoustic_reserve.push(candidate);
            } else if candidate.uses_venturi_treatment() && fast_mode_candidate_filter(&candidate) {
                venturi_filtered.push(candidate);
            }
        }
        acoustic_reserve.sort_by(|a, b| {
            fast_eval_priority_acoustic(a)
                .cmp(&fast_eval_priority_acoustic(b))
                .then_with(|| a.id.cmp(&b.id))
        });
        acoustic_reserve.truncate(run_cfg.fast_nonventuri_reserve);
        let mut filtered: Vec<DesignCandidate> = venturi_filtered
            .into_iter()
            .enumerate()
            .filter(|(idx, _)| idx % run_cfg.fast_stride == 0)
            .map(|(_, c)| c)
            .collect();
        filtered.sort_by(|a, b| {
            fast_eval_priority(a)
                .cmp(&fast_eval_priority(b))
                .then_with(|| a.id.cmp(&b.id))
        });
        let venturi = filtered;
        let keep_nonventuri = acoustic_reserve
            .len()
            .min(run_cfg.fast_nonventuri_reserve)
            .min(run_cfg.fast_max_candidates);
        let keep_venturi = venturi
            .len()
            .min(run_cfg.fast_max_candidates.saturating_sub(keep_nonventuri));
        let mut selected = Vec::with_capacity(keep_nonventuri + keep_venturi);
        let mut vent_iter = venturi.into_iter().take(keep_venturi);
        let mut non_iter = acoustic_reserve.into_iter().take(keep_nonventuri);
        loop {
            let mut pushed = false;
            if let Some(c) = vent_iter.next() {
                selected.push(c);
                pushed = true;
            }
            if let Some(c) = non_iter.next() {
                selected.push(c);
                pushed = true;
            }
            if !pushed {
                break;
            }
        }
        selected
    } else {
        raw_candidates
    };
    let total_candidates = all_candidates.len();
    println!(
        "Mode: {} (M12_FAST={})",
        if run_cfg.fast_mode {
            "FAST regeneration"
        } else {
            "full regeneration"
        },
        if run_cfg.fast_mode { "1" } else { "0" }
    );
    println!(
        "Candidate pool before metrics: {} (stride={}, cap={}, non-venturi reserve={})",
        total_candidates,
        run_cfg.fast_stride,
        run_cfg.fast_max_candidates,
        run_cfg.fast_nonventuri_reserve
    );

    // ── Part 1: Two-concept selection ─────────────────────────────────────────
    println!("\n[1/6] Two-concept parametric selection …");
    let ga_mode = OptimMode::HydrodynamicCavitationSDT;
    let mut option1_pool: Vec<RankedDesign> = Vec::new();
    let mut option2_pool: Vec<RankedDesign> = Vec::new();
    let mut scored_seeds: Vec<(f64, DesignCandidate)> = Vec::new();
    let mut seen_ids: HashSet<String> = HashSet::with_capacity(total_candidates);
    let mut duplicate_ids_skipped: usize = 0;
    let mut unique_candidates: Vec<DesignCandidate> = Vec::with_capacity(total_candidates);
    let metrics_cache = MetricsCache::new(workspace_root.join("target").join("m12_metrics_cache"))?;

    for c in all_candidates {
        if !seen_ids.insert(c.id.clone()) {
            duplicate_ids_skipped += 1;
            continue;
        }
        unique_candidates.push(c);
    }

    let evaluated: Vec<EvaluatedCandidate> = if run_cfg.fast_mode {
        let mut out = Vec::new();
        let mut option1_hits = 0usize;
        let mut option2_hits = 0usize;
        let mut eval_count = 0usize;

        for candidate in unique_candidates {
            if eval_count >= run_cfg.fast_eval_max {
                break;
            }
            let metrics =
                match metrics_cache.get_or_compute(&candidate, || compute_metrics(&candidate)) {
                    Ok(m) => m,
                    Err(_) => continue,
                };
            eval_count += 1;

            if candidate.uses_venturi_treatment()
                && is_selective_report_topology(candidate.topology)
            {
                let mut acoustic_candidate = candidate.clone();
                acoustic_candidate.id = format!("{}-ACS", candidate.id);
                acoustic_candidate.treatment_zone_mode =
                    cfd_optim::TreatmentZoneMode::UltrasoundOnly;
                if let Ok(acoustic_metrics) = metrics_cache
                    .get_or_compute(&acoustic_candidate, || compute_metrics(&acoustic_candidate))
                {
                    if report_eligible_nonventuri(&acoustic_metrics) {
                        option1_hits += 1;
                        out.push(EvaluatedCandidate {
                            candidate: acoustic_candidate,
                            metrics: acoustic_metrics.clone(),
                            ga_score: 0.0,
                            option1_score: Some(score_candidate(
                                &acoustic_metrics,
                                option1_mode,
                                &weights,
                            )),
                            combined_score: None,
                        });
                    }
                }
            }

            let ga_score = score_candidate(&metrics, ga_mode, &weights);
            let option1_score = if !candidate.uses_venturi_treatment()
                && is_selective_report_topology(candidate.topology)
                && report_eligible_nonventuri(&metrics)
            {
                Some(score_candidate(&metrics, option1_mode, &weights))
            } else {
                None
            };

            let mut oncology_score = None;
            if candidate.uses_venturi_treatment()
                && is_selective_report_topology(candidate.topology)
                && report_eligible_venturi_oncology(&metrics)
            {
                oncology_score = Some(score_candidate(&metrics, oncology_mode, &weights));
            }

            if option1_score.is_some() {
                option1_hits += 1;
            }
            if oncology_score.is_some() {
                option2_hits += 1;
            }

            out.push(EvaluatedCandidate {
                candidate,
                metrics,
                ga_score,
                option1_score,
                combined_score: oncology_score,
            });

            if eval_count >= run_cfg.fast_eval_min && option1_hits >= 10 && option2_hits >= 10 {
                break;
            }
        }

        println!(
            "      Fast evaluation budget: evaluated {} candidates (min={}, max={})",
            eval_count, run_cfg.fast_eval_min, run_cfg.fast_eval_max
        );
        out
    } else {
        let cache = metrics_cache.clone();
        unique_candidates
            .into_par_iter()
            .flat_map_iter(|candidate| {
                let mut entries = Vec::new();
                let metrics = match cache.get_or_compute(&candidate, || compute_metrics(&candidate))
                {
                    Ok(m) => m,
                    Err(_) => return entries.into_iter(),
                };

                if candidate.uses_venturi_treatment()
                    && is_selective_report_topology(candidate.topology)
                {
                    let mut acoustic_candidate = candidate.clone();
                    acoustic_candidate.id = format!("{}-ACS", candidate.id);
                    acoustic_candidate.treatment_zone_mode =
                        cfd_optim::TreatmentZoneMode::UltrasoundOnly;
                    if let Ok(acoustic_metrics) = cache.get_or_compute(&acoustic_candidate, || {
                        compute_metrics(&acoustic_candidate)
                    }) {
                        if report_eligible_nonventuri(&acoustic_metrics) {
                            entries.push(EvaluatedCandidate {
                                candidate: acoustic_candidate,
                                metrics: acoustic_metrics.clone(),
                                ga_score: 0.0,
                                option1_score: Some(score_candidate(
                                    &acoustic_metrics,
                                    option1_mode,
                                    &weights,
                                )),
                                combined_score: None,
                            });
                        }
                    }
                }

                let ga_score = score_candidate(&metrics, ga_mode, &weights);
                let option1_score = if !candidate.uses_venturi_treatment()
                    && is_selective_report_topology(candidate.topology)
                    && report_eligible_nonventuri(&metrics)
                {
                    Some(score_candidate(&metrics, option1_mode, &weights))
                } else {
                    None
                };

                let mut oncology_score = None;
                if candidate.uses_venturi_treatment()
                    && is_selective_report_topology(candidate.topology)
                    && report_eligible_venturi_oncology(&metrics)
                {
                    oncology_score = Some(score_candidate(&metrics, oncology_mode, &weights));
                }

                entries.push(EvaluatedCandidate {
                    candidate,
                    metrics,
                    ga_score,
                    option1_score,
                    combined_score: oncology_score,
                });

                entries.into_iter()
            })
            .collect()
    };

    let mut nonventuri_total = 0usize;
    let mut nonventuri_pressure_plate = 0usize;
    let mut nonventuri_fda_main = 0usize;
    let mut nonventuri_eligible = 0usize;
    let mut venturi_total = 0usize;
    let mut venturi_sigma_lt1_any = 0usize;
    let mut sigma_lt1_pressure_ok = 0usize;
    let mut sigma_lt1_plate_ok = 0usize;
    let mut sigma_lt1_fda_main_ok = 0usize;
    let mut sigma_lt1_fda_overall_ok = 0usize;
    let mut sigma_lt1_ppfda_overall_ok = 0usize;
    let mut venturi_pressure_plate_fda = 0usize;
    let mut venturi_sigma_lt1 = 0usize;
    let mut venturi_oncology_eligible = 0usize;
    for ev in evaluated {
        if ev.candidate.uses_venturi_treatment() {
            venturi_total += 1;
            if ev.metrics.cavitation_number.is_finite() && ev.metrics.cavitation_number < 1.0 {
                venturi_sigma_lt1_any += 1;
                if ev.metrics.pressure_feasible {
                    sigma_lt1_pressure_ok += 1;
                }
                if ev.metrics.plate_fits {
                    sigma_lt1_plate_ok += 1;
                }
                if ev.metrics.fda_main_compliant {
                    sigma_lt1_fda_main_ok += 1;
                }
                if ev.metrics.fda_overall_compliant {
                    sigma_lt1_fda_overall_ok += 1;
                }
                if ev.metrics.pressure_feasible
                    && ev.metrics.plate_fits
                    && ev.metrics.fda_overall_compliant
                {
                    sigma_lt1_ppfda_overall_ok += 1;
                }
            }
            if ev.metrics.pressure_feasible
                && ev.metrics.plate_fits
                && ev.metrics.fda_main_compliant
            {
                venturi_pressure_plate_fda += 1;
                if ev.metrics.cavitation_number.is_finite() && ev.metrics.cavitation_number < 1.0 {
                    venturi_sigma_lt1 += 1;
                }
            }
            if report_eligible_venturi_oncology(&ev.metrics) {
                venturi_oncology_eligible += 1;
            }
        } else {
            nonventuri_total += 1;
            if ev.metrics.pressure_feasible && ev.metrics.plate_fits {
                nonventuri_pressure_plate += 1;
            }
            if ev.metrics.fda_main_compliant {
                nonventuri_fda_main += 1;
            }
            if report_eligible_nonventuri(&ev.metrics) {
                nonventuri_eligible += 1;
            }
        }
        if ev.candidate.uses_venturi_treatment() && ev.ga_score > 0.0 {
            scored_seeds.push((ev.ga_score, ev.candidate.clone()));
        }
        if let Some(score) = ev.option1_score {
            option1_pool.push(RankedDesign {
                rank: 0,
                candidate: ev.candidate.clone(),
                metrics: ev.metrics.clone(),
                score,
            });
        }
        if let Some(score) = ev.combined_score {
            option2_pool.push(RankedDesign {
                rank: 0,
                candidate: ev.candidate,
                metrics: ev.metrics,
                score,
            });
        }
    }

    let option1_pool_len = option1_pool.len();
    let option2_pool_len = option2_pool.len();
    let cache_stats = metrics_cache.stats();
    println!(
        "      Candidate deduplication: {} unique, {} duplicate IDs skipped",
        seen_ids.len(),
        duplicate_ids_skipped
    );
    println!(
        "      Gate diagnostics (Option1 selective acoustic): total={} pressure+plate={} fda_main={} eligible={}",
        nonventuri_total, nonventuri_pressure_plate, nonventuri_fda_main, nonventuri_eligible
    );
    println!(
        "      Gate diagnostics (Option2 combined selective venturi): total={} sigma<1(any)={} pressure+plate+fda_main={} sigma<1(pressure+plate+fda_main)={} oncology_eligible={}",
        venturi_total,
        venturi_sigma_lt1_any,
        venturi_pressure_plate_fda,
        venturi_sigma_lt1,
        venturi_oncology_eligible
    );
    println!(
        "      Metric cache: {} hits, {} misses",
        cache_stats.hits, cache_stats.misses
    );
    println!(
        "      Sigma<1 breakdown: pressure={} plate={} fda_main={} fda_overall={} pressure+plate+fda_overall={}",
        sigma_lt1_pressure_ok,
        sigma_lt1_plate_ok,
        sigma_lt1_fda_main_ok,
        sigma_lt1_fda_overall_ok,
        sigma_lt1_ppfda_overall_ok
    );
    let option1_ranked = shortlist_report(option1_pool, 5, "option1_ultrasound")?;
    let option2_ranked = shortlist_report(option2_pool, 5, "option2_venturi_dtcv")?;

    save_top5_json(
        &option1_ranked,
        &out_dir.join("two_concept_option1_ultrasound_top5.json"),
    )?;
    save_top5_json(
        &option2_ranked,
        &out_dir.join("two_concept_option2_venturi_top5.json"),
    )?;

    println!(
        "      Option 1 (selective acoustic): {} (score={:.4})",
        option1_ranked[0].candidate.topology.short(),
        option1_ranked[0].score
    );
    println!(
        "      Option 2 (combined selective venturi): {} (score={:.4}  σ={:.3})",
        option2_ranked[0].candidate.topology.short(),
        option2_ranked[0].score,
        option2_ranked[0].metrics.cavitation_number
    );

    // ── Part 2: GA warm-start search ─────────────────────────────────────────
    println!("\n[2/6] Building parametric seed pool for GA warm-start …");
    println!("      Parametric space: {} candidates", total_candidates);

    scored_seeds.sort_by(|a, b| b.0.total_cmp(&a.0));
    let seeds: Vec<_> = scored_seeds
        .into_iter()
        .take(run_cfg.ga_seed_take)
        .map(|(_, c)| c)
        .collect();

    println!(
        "      Warm-start seeds: {} feasible (top {} selected)",
        seeds.len(),
        run_cfg.ga_seed_take
    );
    println!(
        "      Running HydroSDT GA (pop={}, gen={}) …",
        run_cfg.ga_population, run_cfg.ga_generations
    );

    let ga_result = GeneticOptimizer::new(ga_mode, weights)
        .with_seeds(seeds.clone())
        .with_population(run_cfg.ga_population)
        .with_max_generations(run_cfg.ga_generations)
        .with_rng_seed(M12_GA_HYDRO_SEED)
        .run()?;

    let mut ga_top: Vec<RankedDesign> = ga_result
        .top_designs
        .iter()
        .filter(|design| is_selective_report_topology(design.candidate.topology))
        .cloned()
        .collect();
    if ga_top.is_empty() {
        ga_top = seeds
            .iter()
            .filter(|candidate| is_selective_report_topology(candidate.topology))
            .filter_map(|candidate| {
                compute_metrics(candidate).ok().map(|metrics| RankedDesign {
                    rank: 0,
                    candidate: candidate.clone(),
                    metrics: metrics.clone(),
                    score: score_candidate(&metrics, ga_mode, &weights),
                })
            })
            .collect();
        ga_top.sort_by(|a, b| b.score.total_cmp(&a.score));
        ga_top.truncate(5);
        for (idx, design) in ga_top.iter_mut().enumerate() {
            design.rank = idx + 1;
        }
    }
    if ga_top.is_empty() {
        ga_top = option2_ranked
            .iter()
            .filter(|design| is_selective_report_topology(design.candidate.topology))
            .take(5)
            .map(|design| RankedDesign {
                rank: 0,
                candidate: design.candidate.clone(),
                metrics: design.metrics.clone(),
                score: score_candidate(&design.metrics, ga_mode, &weights),
            })
            .collect();
        ga_top.sort_by(|a, b| b.score.total_cmp(&a.score));
        for (idx, design) in ga_top.iter_mut().enumerate() {
            design.rank = idx + 1;
        }
    }
    if ga_top.is_empty() {
        return Err(
            "HydroSDT comparison track produced no primitive split-sequence designs".into(),
        );
    }
    println!(
        "      GA rank-1: {} (score={:.4}  σ={:.3}  cancer_cav={:.3})",
        ga_top[0].candidate.id,
        ga_top[0].score,
        ga_top[0].metrics.cavitation_number,
        ga_top[0].metrics.cancer_targeted_cavitation
    );
    save_top5_json(&ga_top, &out_dir.join("ga_hydrosdt_top5.json"))?;
    println!("      Saved: ga_hydrosdt_top5.json");

    // ── Part 3: Multi-fidelity validation ────────────────────────────────────
    let mut validation_rows: Vec<ValidationRow> = Vec::new();
    let option2_robustness: Vec<RobustnessReport>;
    if run_cfg.fast_mode {
        println!("\n[3/6] Fast mode: skipping multi-fidelity validation");
        println!("[4/6] Fast mode: skipping robustness sweep");
        option2_robustness = Vec::new();
    } else {
        println!("\n[3/6] Running 2D FVM + 3D FEM validation on selected venturi designs …");

        let combined_validation = validate_venturi_candidate(
            "Option 2 Combined Selective Venturi",
            &option2_ranked[0],
            &out_dir,
        )?;
        validation_rows.push(combined_validation.clone());
        validation_rows.push(validate_venturi_candidate(
            "GA HydroSDT Venturi",
            &ga_top[0],
            &out_dir,
        )?);

        // ── Part 4: Robustness screening ─────────────────────────────────────
        println!("\n[4/6] Running robustness screening on Option 2 shortlist …");
        option2_robustness = option2_ranked
            .iter()
            .map(|d| {
                robustness_sweep(
                    &d.candidate,
                    oncology_mode,
                    &weights,
                    &STANDARD_PERTURBATIONS,
                )
            })
            .collect();
        let robust_count = option2_robustness.iter().filter(|r| r.is_robust).count();
        println!(
            "      Robust candidates: {}/{}",
            robust_count,
            option2_robustness.len()
        );
        let option2_robustness_json = serde_json::to_string_pretty(&option2_robustness)?;
        std::fs::write(
            out_dir.join("option2_combined_robustness_top5.json"),
            option2_robustness_json,
        )?;
    }

    // ── Part 5: Figure generation ─────────────────────────────────────────────
    println!("\n[5/6] Generating report figures …");

    // Figure 4 — best selective acoustic design
    save_figure(
        &option1_ranked[0].candidate,
        &figures_dir.join("selected_option1_schematic.svg"),
        "Figure 4 (Option 1 selective acoustic)",
    );

    // Figure 5 — best combined selective venturi design
    save_figure(
        &option2_ranked[0].candidate,
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Figure 5 (Option 2 combined selective venturi)",
    );

    // Figure 6 — best HydroSDT GA design
    save_figure(
        &ga_top[0].candidate,
        &figures_dir.join("top_hydrosdt_schematic.svg"),
        "Figure 6 (HydroSDT GA rank-1)",
    );

    // ── Part 6: Canonical markdown + JSON ────────────────────────────────────
    println!("\n[6/6] Writing canonical report and summary JSON …");

    // Two-concept selection summary
    let toplines = vec![
        concept_topline(
            "Option 1: Selective acoustic center treatment",
            "SelectiveAcousticTherapy",
            &option1_ranked[0],
        ),
        concept_topline(
            "Option 2: Selective venturi treatment ranked by combined score",
            "CombinedSdtLeukapheresis",
            &option2_ranked[0],
        ),
    ];
    let topline_json = serde_json::to_string_pretty(&toplines)?;
    std::fs::write(
        out_dir.join("two_concept_selection_summary.json"),
        topline_json,
    )?;

    let canonical_results = workspace_root.join("report").join("milestone12_results.md");
    write_milestone12_results(
        total_candidates,
        option1_pool_len,
        option2_pool_len,
        &option1_ranked,
        &option2_ranked,
        &ga_top[0],
        &validation_rows,
        &option2_robustness,
        &canonical_results,
    )?;
    println!("Canonical report: {}", canonical_results.display());

    let narrative_artifacts = write_milestone12_narrative_report(
        &workspace_root,
        &canonical_results,
        &Milestone12NarrativeInput {
            total_candidates,
            option1_pool_len,
            option2_pool_len,
            option1_ranked: &option1_ranked,
            option2_ranked: &option2_ranked,
            ga_top: &ga_top,
            validation_rows: &validation_rows,
            option2_robustness: &option2_robustness,
            ga_best_per_gen: &ga_result.best_per_gen,
            fast_mode: run_cfg.fast_mode,
        },
    )?;
    println!(
        "Narrative report: {}",
        narrative_artifacts.narrative_path.display()
    );
    println!(
        "Figure manifest: {} ({} figures)",
        narrative_artifacts.figure_manifest_path.display(),
        narrative_artifacts.figure_count
    );

    println!("\n=== Outputs written to {} ===", out_dir.display());
    println!("Figures written to {}", figures_dir.display());

    Ok(())
}

// ── Figure saving helper ──────────────────────────────────────────────────────

fn save_figure(candidate: &DesignCandidate, path: &Path, label: &str) {
    match save_schematic_svg(candidate, path) {
        Ok(()) => println!("      ✓ {label}  →  {}", path.display()),
        Err(e) => eprintln!("      ✗ {label}  FAILED: {e}"),
    }
}

// ── Concept topline builder ───────────────────────────────────────────────────

fn concept_topline(concept: &str, mode: &str, d: &RankedDesign) -> ConceptTopline {
    ConceptTopline {
        concept: concept.to_string(),
        mode: mode.to_string(),
        candidate_id: d.candidate.id.clone(),
        score: d.score,
        topology: d.candidate.topology.short().to_string(),
        treatment_zone_mode: d.metrics.treatment_zone_mode.clone(),
        active_venturi_throat_count: d.metrics.active_venturi_throat_count,
        cavitation_number: d.metrics.cavitation_number,
        cancer_center_fraction: d.metrics.cancer_center_fraction,
        therapeutic_window_score: d.metrics.therapeutic_window_score,
        hemolysis_index_per_pass: d.metrics.hemolysis_index_per_pass,
        wall_shear_p95_pa: d.metrics.wall_shear_p95_pa,
        mean_residence_time_s: d.metrics.mean_residence_time_s,
        total_ecv_ml: d.metrics.total_ecv_ml,
        therapy_channel_fraction: d.metrics.therapy_channel_fraction,
    }
}
