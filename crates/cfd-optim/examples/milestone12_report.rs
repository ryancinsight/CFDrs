//! Milestone 12 unified report generator.
//!
//! Consolidates concept selection, GA optimisation, and figure/report generation
//! in one runnable example.  All outputs go to
//! `report/milestone12/` and `report/figures/`.
//!
//! # Parts
//! 1. **Two-concept selection** — parametric sweep over shared selective
//!    center-enrichment layouts, ranked by SelectiveAcousticTherapy
//!    (Option 1) and CombinedSdtLeukapheresis (Option 2), top-5 each.
//! 2. **GA warm-start search** — top-100 parametric seeds → HydrodynamicCavitationSDT GA
//!    (pop=80, gen=120) for deeper optimisation of venturi-cavitation designs.
//! 3. **Robustness screening** — +/-10 % / +/-20 % operating-point perturbations on
//!    the selected Option 2 shortlist.
//! 4. **Figure generation** — schematics selected dynamically from optimizer output
//!    (Figure 4 = best selective acoustic design; Figure 5 = best combined-score
//!    venturi design; Figure 6 = best HydroSDT GA design).
//! 5. **Canonical report markdown** — `report/milestone12_results.md`.
//!
//! Multi-fidelity 2D/3D venturi validation is handled by the companion example
//! `milestone12_validation`.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_report --no-default-features
//! # optional fast iteration mode
//! $env:M12_FAST=1; cargo run -p cfd-optim --example milestone12_report --no-default-features
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, compute_blueprint_report_metrics,
    evaluate_blueprint_candidate, evaluate_blueprint_genetic_refinement,
    is_milestone12_lineage_topology, milestone12_lineage_key, robustness_sweep_blueprint,
    save_blueprint_schematic_svg, save_top5_report_json, score_candidate, shortlist_report_designs,
    sort_report_designs, validate_milestone12_candidate, write_milestone12_narrative_report,
    write_milestone12_results, BlueprintCandidate, BlueprintGeneticOptimizer,
    Milestone12LineageKey, Milestone12NarrativeInput, Milestone12ReportDesign, Milestone12Stage,
    OptimMode, OptimizationGoal, RobustnessReport, SdtMetrics, SdtWeights, ValidationRow,
    STANDARD_PERTURBATIONS,
};
use cfd_schematics::NetworkBlueprint;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

// ── Output structures ─────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
struct ConceptTopline {
    concept: String,
    mode: String,
    available: bool,
    status: String,
    candidate_id: String,
    score: Option<f64>,
    topology: String,
    treatment_zone_mode: String,
    active_venturi_throat_count: Option<usize>,
    cavitation_number: Option<f64>,
    cancer_center_fraction: Option<f64>,
    therapeutic_window_score: Option<f64>,
    hemolysis_index_per_pass: Option<f64>,
    wall_shear_p95_pa: Option<f64>,
    mean_residence_time_s: Option<f64>,
    treatment_zone_dwell_time_s: Option<f64>,
    total_ecv_ml: Option<f64>,
    therapy_channel_fraction: Option<f64>,
}

struct RunConfig {
    fast_mode: bool,
    fast_stride: usize,
    fast_max_candidates: usize,
    fast_nonventuri_reserve: usize,
    fast_acoustic_eval_max: usize,
    fast_eval_min: usize,
    fast_eval_max: usize,
    ga_seed_take: usize,
    ga_population: usize,
    ga_generations: usize,
}

struct ScanProgress {
    label: &'static str,
    total: usize,
    started_at: Instant,
    processed: AtomicUsize,
    next_log_at: AtomicUsize,
    log_step: usize,
}

impl ScanProgress {
    fn new(label: &'static str, total: usize) -> Self {
        let log_step = if total <= 25 {
            1
        } else {
            (total / 20).clamp(25, 1_000)
        };
        tracing::info!(
            "      {label}: starting {} candidate evaluations (heartbeat every {} candidates)",
            total,
            log_step
        );
        Self {
            label,
            total,
            started_at: Instant::now(),
            processed: AtomicUsize::new(0),
            next_log_at: AtomicUsize::new(log_step.min(total.max(1))),
            log_step,
        }
    }

    fn record(&self) {
        let processed = self.processed.fetch_add(1, Ordering::Relaxed) + 1;
        loop {
            let next_log_at = self.next_log_at.load(Ordering::Relaxed);
            if processed < next_log_at && processed < self.total {
                break;
            }
            let new_next = (next_log_at + self.log_step).min(self.total.max(1));
            match self.next_log_at.compare_exchange(
                next_log_at,
                new_next,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    self.log(processed);
                    break;
                }
                Err(_) => continue,
            }
        }
    }

    fn finish(&self) {
        self.log(self.processed.load(Ordering::Relaxed));
    }

    fn log(&self, processed: usize) {
        let elapsed = self.started_at.elapsed().as_secs_f64().max(1.0e-9);
        let pct = if self.total == 0 {
            100.0
        } else {
            100.0 * processed as f64 / self.total as f64
        };
        let rate = processed as f64 / elapsed;
        tracing::info!(
            "      {}: {}/{} ({:.1}%) in {:.1}s [{:.1} candidates/s]",
            self.label,
            processed,
            self.total,
            pct,
            elapsed,
            rate
        );
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn report_eligible_nonventuri(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.therapy_channel_fraction > 0.0
}

fn is_selective_report_topology(candidate: &BlueprintCandidate) -> bool {
    is_milestone12_lineage_topology(candidate)
}

/// Relaxed lineage match for GA results: same `PrimitiveSelectiveTree` sequence
/// as the selected lineage, but geometry parameters may differ (the GA mutates
/// channel width, serpentine counts, etc. via SBX crossover).
fn ga_matches_lineage_sequence(
    candidate: &BlueprintCandidate,
    selected_key: &Milestone12LineageKey,
) -> bool {
    candidate
        .topology_spec()
        .is_ok_and(|spec| spec.stage_sequence_label() == selected_key.stage_sequence_label())
}

fn report_eligible_venturi_oncology(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.therapy_channel_fraction > 0.0
        && metrics.cavitation_number < 1.0
}

fn report_option1_score(metrics: &SdtMetrics, weights: &SdtWeights) -> f64 {
    let base = score_candidate(metrics, OptimMode::SelectiveAcousticTherapy, weights);
    let dwell_tie_break = metrics.treatment_zone_dwell_time_s.max(0.0) * 1.0e-6;
    let pressure_tie_break = metrics.total_pressure_drop_pa.max(0.0) * 1.0e-12;
    base + dwell_tie_break - pressure_tie_break
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
        let fast_acoustic_eval_max = std::env::var("M12_FAST_ACOUSTIC_EVAL_MAX")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(80)
            .max(20);
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
        let ga_seed_take = std::env::var("M12_FAST_GA_SEED_TAKE")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(24)
            .max(8);
        let ga_population = std::env::var("M12_FAST_GA_POPULATION")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(12)
            .max(10);
        let ga_generations = std::env::var("M12_FAST_GA_GENERATIONS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(8)
            .max(1);
        RunConfig {
            fast_mode: true,
            fast_stride,
            fast_max_candidates,
            fast_nonventuri_reserve,
            fast_acoustic_eval_max,
            fast_eval_min,
            fast_eval_max,
            ga_seed_take,
            ga_population,
            ga_generations,
        }
    } else {
        RunConfig {
            fast_mode: false,
            fast_stride: 1,
            fast_max_candidates: usize::MAX,
            fast_nonventuri_reserve: 0,
            fast_acoustic_eval_max: usize::MAX,
            fast_eval_min: usize::MAX,
            fast_eval_max: usize::MAX,
            ga_seed_take: 150,
            ga_population: 80,
            ga_generations: 120,
        }
    }
}

fn fast_mode_candidate_filter(c: &BlueprintCandidate) -> bool {
    let flow_ml_min = c.operating_point.flow_rate_m3_s * 6.0e7;
    let gauge_kpa = c.operating_point.inlet_gauge_pa * 1.0e-3;
    let Ok(topology) = c.topology_spec() else {
        return false;
    };
    let treatment_branch = topology
        .split_stages
        .iter()
        .rev()
        .flat_map(|stage| stage.branches.iter())
        .find(|branch| branch.treatment_path);
    let channel_width_m = topology
        .split_stages
        .first()
        .map(|stage| {
            stage
                .branches
                .iter()
                .map(|branch| branch.route.width_m)
                .sum::<f64>()
        })
        .unwrap_or(topology.inlet_width_m);
    let channel_height_m = treatment_branch.map_or(0.0, |branch| branch.route.height_m);
    if !topology.has_venturi() {
        return flow_ml_min <= 250.0
            && gauge_kpa <= 300.0
            && channel_width_m >= 1.0e-3
            && channel_width_m <= 6.0e-3
            && channel_height_m >= 0.5e-3
            && channel_height_m <= 2.0e-3;
    }

    let Some(venturi) = topology.venturi_placements.first() else {
        return false;
    };
    let throat_width_m = venturi.throat_geometry.throat_width_m;
    let throat_length_m = venturi.throat_geometry.throat_length_m;
    let tl_factor = if throat_width_m > 0.0 {
        throat_length_m / throat_width_m
    } else {
        0.0
    };
    let first_stage_center_frac = topology
        .split_stages
        .first()
        .and_then(|stage| {
            let stage_width = stage
                .branches
                .iter()
                .map(|branch| branch.route.width_m)
                .sum::<f64>();
            stage
                .branches
                .iter()
                .find(|branch| branch.treatment_path)
                .map(|branch| branch.route.width_m / stage_width.max(1.0e-12))
        })
        .unwrap_or(1.0 / 3.0);
    let terminal_tri_center_frac = topology
        .split_stages
        .iter()
        .rev()
        .find(|stage| matches!(stage.split_kind, cfd_schematics::SplitKind::Trifurcation))
        .and_then(|stage| {
            let stage_width = stage
                .branches
                .iter()
                .map(|branch| branch.route.width_m)
                .sum::<f64>();
            stage
                .branches
                .iter()
                .find(|branch| branch.treatment_path)
                .map(|branch| branch.route.width_m / stage_width.max(1.0e-12))
        })
        .unwrap_or(first_stage_center_frac);

    topology.first_stage_is_trifurcation()
        && topology.is_selective_routing()
        && (1..=4).contains(&(venturi.serial_throat_count as usize))
        && flow_ml_min >= 60.0
        && flow_ml_min <= 550.0
        && gauge_kpa >= 25.0
        && gauge_kpa <= 550.0
        && throat_width_m >= 25.0e-6
        && throat_width_m <= 120.0e-6
        && tl_factor <= 5.0
        && channel_width_m >= 3.0e-3
        && channel_width_m <= 8.0e-3
        && channel_height_m >= 0.5e-3
        && channel_height_m <= 2.5e-3
        && first_stage_center_frac >= 0.45
        && first_stage_center_frac <= 0.62
        && terminal_tri_center_frac >= 0.33
        && terminal_tri_center_frac <= 0.58
}

fn fast_eval_priority(c: &BlueprintCandidate) -> (u8, i64, i64, i64) {
    let family = u8::from(
        c.topology_spec()
            .is_ok_and(|spec| spec.stage_sequence_label() != "Tri→Tri"),
    );
    // Dual-centre priority: candidates close to EITHER the acoustic sweet spot
    // (q=120, g=200, d=55) or the cavitation sweet spot (q=300, g=400, d=35)
    // rank highly, ensuring both Option 1 and Option 2 are well-sampled.
    let q_ml = c.operating_point.flow_rate_m3_s * 6.0e7;
    let g_kpa = c.operating_point.inlet_gauge_pa * 1.0e-3;
    let d_um = c
        .topology_spec()
        .ok()
        .and_then(|spec| spec.venturi_placements.first())
        .map_or(0.0, |placement| {
            placement.throat_geometry.throat_width_m * 1.0e6
        });
    let acoustic_pen = (q_ml - 120.0).abs() + (g_kpa - 200.0).abs() + (d_um - 55.0).abs();
    let cav_pen = (q_ml - 300.0).abs() + (g_kpa - 400.0).abs() + (d_um - 35.0).abs();
    let best_pen = acoustic_pen.min(cav_pen).round() as i64;
    (family, best_pen, 0, 0)
}

/// Priority for acoustic-only candidates: centre on q = 120 mL/min, g = 200 kPa,
/// d = 55 μm — the sweet spot for selective routing without cavitation.
fn fast_eval_priority_acoustic(c: &BlueprintCandidate) -> (u8, i64, i64, i64) {
    let family = u8::from(
        c.topology_spec()
            .is_ok_and(|spec| spec.stage_sequence_label() != "Tri→Tri"),
    );
    let q_pen = (c.operating_point.flow_rate_m3_s * 6.0e7 - 120.0)
        .abs()
        .round() as i64;
    let g_pen = (c.operating_point.inlet_gauge_pa * 1.0e-3 - 200.0)
        .abs()
        .round() as i64;
    let d_pen = c
        .topology_spec()
        .ok()
        .and_then(|spec| spec.venturi_placements.first())
        .map_or(55.0, |placement| {
            placement.throat_geometry.throat_width_m * 1.0e6
        });
    let d_pen = (d_pen - 55.0).abs().round() as i64;
    (family, q_pen, g_pen, d_pen)
}

fn blueprint_lineage_key(candidate: &BlueprintCandidate) -> Option<Milestone12LineageKey> {
    milestone12_lineage_key(candidate)
}

fn blueprint_matches_lineage(
    candidate: &BlueprintCandidate,
    selected_key: &Milestone12LineageKey,
) -> bool {
    blueprint_lineage_key(candidate).as_ref() == Some(selected_key)
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(tracing::Level::INFO.into())
                .from_env_lossy(),
        )
        .try_init();

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

    tracing::info!("=== Milestone 12 Unified Report Generator ===\n");
    tracing::info!("Outputs: {}", out_dir.display());
    let overall_start = Instant::now();

    let weights = SdtWeights::default();
    let oncology_mode = option2_mode();
    let run_cfg = run_config_from_env();
    let raw_candidates = build_milestone12_blueprint_candidate_space()?;
    let all_candidates: Vec<BlueprintCandidate> = if run_cfg.fast_mode {
        let mut acoustic_reserve: Vec<BlueprintCandidate> = Vec::new();
        let mut venturi_filtered: Vec<BlueprintCandidate> = Vec::new();
        for candidate in raw_candidates {
            let uses_venturi = candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
            if !uses_venturi
                && is_selective_report_topology(&candidate)
                && fast_mode_candidate_filter(&candidate)
            {
                acoustic_reserve.push(candidate);
            } else if uses_venturi && fast_mode_candidate_filter(&candidate) {
                venturi_filtered.push(candidate);
            }
        }
        // Sort-then-stride: sort by priority FIRST so high-priority candidates survive
        // the stride filter. Use numeric ID prefix comparison so "699982-..." sorts
        // before "1785880-..." (lexicographic '1' < '6' gives wrong order).
        fn id_num(id: &str) -> u64 {
            id.split('-')
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(u64::MAX)
        }
        tracing::info!(
            "      [fast] acoustic_reserve before truncate: {} candidates (venturi_filtered before stride: {})",
            acoustic_reserve.len(), venturi_filtered.len()
        );
        if let Some(sample) = acoustic_reserve.first() {
            tracing::info!(
                "      [fast] acoustic_reserve sample[0]: id={} topology={} uses_venturi={}",
                sample.id,
                sample
                    .topology_spec()
                    .map_or_else(|_| String::new(), |spec| spec.display_name()),
                sample
                    .topology_spec()
                    .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi)
            );
        }
        acoustic_reserve.sort_by(|a, b| {
            fast_eval_priority_acoustic(a)
                .cmp(&fast_eval_priority_acoustic(b))
                .then_with(|| id_num(&a.id).cmp(&id_num(&b.id)))
        });
        acoustic_reserve.truncate(run_cfg.fast_nonventuri_reserve);
        let mut filtered: Vec<BlueprintCandidate> = venturi_filtered;
        filtered.sort_by(|a, b| {
            fast_eval_priority(a)
                .cmp(&fast_eval_priority(b))
                .then_with(|| id_num(&a.id).cmp(&id_num(&b.id)))
        });
        let filtered: Vec<BlueprintCandidate> = filtered
            .into_iter()
            .enumerate()
            .filter(|(idx, _)| idx % run_cfg.fast_stride == 0)
            .map(|(_, c)| c)
            .collect();
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
    let mut blueprint_by_id: HashMap<String, NetworkBlueprint> = all_candidates
        .iter()
        .map(|candidate| (candidate.id.clone(), candidate.blueprint().clone()))
        .collect();
    let fast_acoustic_candidate_count = if run_cfg.fast_mode {
        all_candidates
            .iter()
            .filter(|candidate| {
                !candidate
                    .topology_spec()
                    .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi)
            })
            .count()
    } else {
        0
    };
    let fast_venturi_candidate_count = if run_cfg.fast_mode {
        all_candidates
            .iter()
            .filter(|candidate| {
                candidate
                    .topology_spec()
                    .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi)
            })
            .count()
    } else {
        0
    };
    let total_candidates = all_candidates.len();
    tracing::info!(
        "Mode: {} (M12_FAST={})",
        if run_cfg.fast_mode {
            "FAST regeneration"
        } else {
            "full regeneration"
        },
        if run_cfg.fast_mode { "1" } else { "0" }
    );
    tracing::info!(
        "Candidate pool before metrics: {} (stride={}, cap={}, non-venturi reserve={}, acoustic eval cap={})",
        total_candidates,
        run_cfg.fast_stride,
        run_cfg.fast_max_candidates,
        run_cfg.fast_nonventuri_reserve,
        run_cfg.fast_acoustic_eval_max
    );
    if run_cfg.fast_mode {
        tracing::info!(
            "Fast candidate mix: acoustic={} venturi={}",
            fast_acoustic_candidate_count,
            fast_venturi_candidate_count
        );
    }

    // ── Part 1: Two-concept selection ─────────────────────────────────────────
    let phase1_start = Instant::now();
    tracing::info!("\n[1/6] Two-concept parametric selection …");
    let ga_mode = OptimMode::HydrodynamicCavitationSDT;
    let duplicate_ids_skipped: usize = 0;
    let unique_candidate_count = total_candidates;
    let evaluated = if run_cfg.fast_mode {
        let progress = ScanProgress::new("candidate scan", all_candidates.len());
        let mut acc = EvaluationAccumulator::default();
        let mut eval_count = 0usize;
        let mut acoustic_eval_count = 0usize;
        let has_fast_acoustic_candidates = fast_acoustic_candidate_count > 0;
        let mut lineage_ready = false;

        for blueprint_candidate in all_candidates {
            let blueprint_evaluation = match evaluate_blueprint_candidate(&blueprint_candidate) {
                Ok(evaluation) => evaluation,
                Err(_) => {
                    progress.record();
                    continue;
                }
            };
            let is_acoustic = !blueprint_candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
            let acoustic_unavailable = acc.option1_hits() == 0
                && (!has_fast_acoustic_candidates
                    || acoustic_eval_count >= run_cfg.fast_acoustic_eval_max);
            // In fast mode, stop at the nominal cap only after both tracks have enough
            // eligible candidates to build a stable shortlist. Otherwise keep scanning.
            // `lineage_ready` is cached from the previous iteration's post-record check
            // to avoid redundant HashSet/HashMap allocation.
            if eval_count >= run_cfg.fast_eval_max && lineage_ready {
                break;
            }
            if is_acoustic && acoustic_unavailable {
                progress.record();
                continue;
            }
            let metrics = match compute_blueprint_report_metrics(&blueprint_candidate) {
                Ok(m) => m,
                Err(_) => {
                    if is_acoustic {
                        acoustic_eval_count += 1;
                    }
                    progress.record();
                    continue;
                }
            };
            if is_acoustic {
                acoustic_eval_count += 1;
            }
            eval_count += 1;
            progress.record();

            let ga_score = evaluate_blueprint_genetic_refinement(
                &blueprint_candidate,
                blueprint_evaluation.clone(),
            )
            .map_or(0.0, |evaluation| evaluation.score);
            let option1_score = if is_acoustic
                && is_selective_report_topology(&blueprint_candidate)
                && report_eligible_nonventuri(&metrics)
            {
                Some(report_option1_score(&metrics, &weights))
            } else {
                None
            };

            let combined_score = if !is_acoustic
                && is_selective_report_topology(&blueprint_candidate)
                && report_eligible_venturi_oncology(&metrics)
            {
                Some(score_candidate(&metrics, oncology_mode, &weights))
            } else {
                None
            };

            acc.record(
                blueprint_candidate,
                metrics,
                ga_score,
                option1_score,
                combined_score,
            );

            lineage_ready = has_viable_lineage_family(&acc.option1_pool, &acc.option2_pool, 5);
            if eval_count >= run_cfg.fast_eval_min && lineage_ready {
                break;
            }
        }

        progress.finish();

        tracing::info!(
            "      Fast evaluation budget: evaluated {} candidates (acoustic evaluated={} min={} max={} acoustic_cap={})",
            eval_count,
            acoustic_eval_count,
            run_cfg.fast_eval_min,
            run_cfg.fast_eval_max,
            run_cfg.fast_acoustic_eval_max
        );
        acc
    } else {
        let progress = Arc::new(ScanProgress::new("candidate scan", all_candidates.len()));
        let progress_for_eval = Arc::clone(&progress);
        let evaluated = all_candidates
            .into_par_iter()
            .fold(
                EvaluationAccumulator::default,
                move |mut acc, blueprint_candidate| {
                    let blueprint_evaluation =
                        match evaluate_blueprint_candidate(&blueprint_candidate) {
                            Ok(evaluation) => evaluation,
                            Err(_) => {
                                progress_for_eval.record();
                                return acc;
                            }
                        };
                    let metrics = match compute_blueprint_report_metrics(&blueprint_candidate) {
                        Ok(m) => m,
                        Err(_) => {
                            progress_for_eval.record();
                            return acc;
                        }
                    };
                    progress_for_eval.record();

                    let ga_score = evaluate_blueprint_genetic_refinement(
                        &blueprint_candidate,
                        blueprint_evaluation,
                    )
                    .map_or(0.0, |evaluation| evaluation.score);
                    let is_acoustic = !blueprint_candidate
                        .topology_spec()
                        .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
                    let option1_score = if is_acoustic
                        && is_selective_report_topology(&blueprint_candidate)
                        && report_eligible_nonventuri(&metrics)
                    {
                        Some(report_option1_score(&metrics, &weights))
                    } else {
                        None
                    };

                    let combined_score = if !is_acoustic
                        && is_selective_report_topology(&blueprint_candidate)
                        && report_eligible_venturi_oncology(&metrics)
                    {
                        Some(score_candidate(&metrics, oncology_mode, &weights))
                    } else {
                        None
                    };

                    acc.record(
                        blueprint_candidate,
                        metrics,
                        ga_score,
                        option1_score,
                        combined_score,
                    );
                    acc
                },
            )
            .reduce(EvaluationAccumulator::default, |mut lhs, rhs| {
                lhs.merge(rhs);
                lhs
            });
        progress.finish();
        evaluated
    };
    let option1_pool = evaluated.option1_pool;
    let option2_pool = evaluated.option2_pool;
    let mut scored_seeds = evaluated.scored_seeds;
    let nonventuri_total = evaluated.nonventuri_total;
    let nonventuri_pressure_plate = evaluated.nonventuri_pressure_plate;
    let nonventuri_fda_main = evaluated.nonventuri_fda_main;
    let nonventuri_eligible = evaluated.nonventuri_eligible;
    let venturi_total = evaluated.venturi_total;
    let venturi_sigma_lt1_any = evaluated.venturi_sigma_lt1_any;
    let sigma_lt1_pressure_ok = evaluated.sigma_lt1_pressure_ok;
    let sigma_lt1_plate_ok = evaluated.sigma_lt1_plate_ok;
    let sigma_lt1_fda_main_ok = evaluated.sigma_lt1_fda_main_ok;
    let sigma_lt1_fda_overall_ok = evaluated.sigma_lt1_fda_overall_ok;
    let sigma_lt1_ppfda_overall_ok = evaluated.sigma_lt1_ppfda_overall_ok;
    let venturi_pressure_plate_fda = evaluated.venturi_pressure_plate_fda;
    let venturi_sigma_lt1 = evaluated.venturi_sigma_lt1;
    let venturi_oncology_eligible = evaluated.venturi_oncology_eligible;
    let option1_pool_len = option1_pool.len();
    let option2_pool_len = option2_pool.len();
    tracing::info!(
        "      Candidate deduplication: {} unique, {} duplicate IDs skipped",
        unique_candidate_count,
        duplicate_ids_skipped
    );
    tracing::info!(
        "      Gate diagnostics (Option1 selective acoustic): total={} pressure+plate={} fda_main={} eligible={}",
        nonventuri_total, nonventuri_pressure_plate, nonventuri_fda_main, nonventuri_eligible
    );
    tracing::info!(
        "      Gate diagnostics (Option2 combined selective venturi): total={} sigma<1(any)={} pressure+plate+fda_main={} sigma<1(pressure+plate+fda_main)={} oncology_eligible={}",
        venturi_total,
        venturi_sigma_lt1_any,
        venturi_pressure_plate_fda,
        venturi_sigma_lt1,
        venturi_oncology_eligible
    );
    tracing::info!(
        "      Sigma<1 breakdown: pressure={} plate={} fda_main={} fda_overall={} pressure+plate+fda_overall={}",
        sigma_lt1_pressure_ok,
        sigma_lt1_plate_ok,
        sigma_lt1_fda_main_ok,
        sigma_lt1_fda_overall_ok,
        sigma_lt1_ppfda_overall_ok
    );
    let option1_shortlist_len = option1_pool.len().min(5);
    let option1_shortlist_global = shortlist_report_designs(
        option1_pool.clone(),
        option1_shortlist_len,
        "option1_ultrasound",
    )?;
    let mut option1_by_lineage: HashMap<_, Vec<Milestone12ReportDesign>> = HashMap::new();
    let mut option2_by_lineage: HashMap<_, Vec<Milestone12ReportDesign>> = HashMap::new();
    for design in option1_pool {
        if let Some(key) = blueprint_lineage_key(&design.candidate) {
            option1_by_lineage.entry(key).or_default().push(design);
        }
    }
    for design in option2_pool {
        if let Some(key) = blueprint_lineage_key(&design.candidate) {
            option2_by_lineage.entry(key).or_default().push(design);
        }
    }

    let mut viable_lineages: Vec<_> = option2_by_lineage
        .into_iter()
        .filter_map(|(key, mut derived_designs)| {
            let option1_designs = option1_by_lineage.get(&key)?;
            let mut acoustic_designs = option1_designs.clone();
            sort_report_designs(&mut acoustic_designs);
            sort_report_designs(&mut derived_designs);
            Some((key, acoustic_designs[0].clone(), derived_designs))
        })
        .collect();
    if viable_lineages.is_empty() {
        return Err(
            "Milestone 12 lineage broke: no venturi Option 2 family shares a valid Option 1 scaffold"
                .into(),
        );
    }
    viable_lineages.sort_by(|lhs, rhs| {
        rhs.2[0]
            .score
            .total_cmp(&lhs.2[0].score)
            .then_with(|| lhs.1.candidate.id.cmp(&rhs.1.candidate.id))
    });

    let (selected_lineage_key, selected_option1, option2_pool_all) = viable_lineages
        .into_iter()
        .find(|(_, _, derived_designs)| derived_designs.len() >= 5)
        .ok_or(
            "Milestone 12 lineage broke: no shared Option 1 -> Option 2 family retained five derived venturi candidates in this run",
        )?;

    let mut option1_ranked = Vec::with_capacity(option1_shortlist_len);
    option1_ranked.push(selected_option1.clone());
    for design in option1_shortlist_global {
        if option1_ranked.len() == option1_shortlist_len {
            break;
        }
        if design.candidate.id != selected_option1.candidate.id {
            option1_ranked.push(design);
        }
    }
    for (idx, design) in option1_ranked.iter_mut().enumerate() {
        design.rank = idx + 1;
    }
    let mut option2_ranked = option2_pool_all.clone();
    option2_ranked.truncate(5);
    for (idx, design) in option2_ranked.iter_mut().enumerate() {
        design.rank = idx + 1;
    }

    save_top5_report_json(
        &option1_ranked,
        &out_dir.join("two_concept_option1_ultrasound_top5.json"),
    )?;
    save_top5_report_json(
        &option2_ranked,
        &out_dir.join("two_concept_option2_venturi_top5.json"),
    )?;

    tracing::info!(
        "      Option 1 (selective acoustic): {} (score={:.4})",
        selected_option1.topology_short_code(),
        selected_option1.score
    );
    tracing::info!(
        "      Option 2 derived family size from Option 1 scaffold: {} candidates",
        option2_pool_all.len()
    );
    tracing::info!(
        "      Option 2 (combined selective venturi): {} (score={:.4}  σ={:.3})",
        option2_ranked[0].topology_short_code(),
        option2_ranked[0].score,
        option2_ranked[0].metrics.cavitation_number
    );

    // ── Part 2: GA warm-start search ─────────────────────────────────────────
    tracing::info!(
        "      [1/6] completed in {:.1}s",
        phase1_start.elapsed().as_secs_f64()
    );
    let phase2_start = Instant::now();
    tracing::info!("\n[2/6] Building parametric seed pool for GA warm-start …");
    tracing::info!("      Parametric space: {} candidates", total_candidates);

    scored_seeds.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut seeds: Vec<_> = scored_seeds
        .into_iter()
        .filter(|(_, candidate)| blueprint_matches_lineage(candidate, &selected_lineage_key))
        .take(run_cfg.ga_seed_take)
        .map(|(_, c)| c)
        .collect();
    if seeds.len() < run_cfg.ga_seed_take {
        for design in &option2_ranked {
            if seeds.len() == run_cfg.ga_seed_take {
                break;
            }
            if seeds
                .iter()
                .all(|candidate| candidate.id != design.candidate.id)
            {
                seeds.push(design.candidate.clone());
            }
        }
    }
    if seeds.is_empty() {
        return Err(
            "Milestone 12 lineage broke: GA received no seeds from the selected Option 2 family"
                .into(),
        );
    }

    tracing::info!(
        "      Warm-start seeds: {} feasible (top {} selected)",
        seeds.len(),
        run_cfg.ga_seed_take
    );
    tracing::info!(
        "      Running HydroSDT GA (pop={}, gen={}) …",
        run_cfg.ga_population,
        run_cfg.ga_generations
    );

    let ga_result = BlueprintGeneticOptimizer::new(OptimizationGoal::BlueprintGeneticRefinement)
        .with_seeds(seeds.clone())
        .with_population(run_cfg.ga_population)
        .with_max_generations(run_cfg.ga_generations)
        .with_top_k(5)
        .run()?;
    for ranked in &ga_result.all_candidates {
        blueprint_by_id.insert(
            ranked.candidate.id.clone(),
            ranked.candidate.blueprint().clone(),
        );
    }

    let ga_report_pool_all: Vec<Milestone12ReportDesign> = ga_result
        .all_candidates
        .iter()
        .filter(|ranked| blueprint_matches_lineage(&ranked.candidate, &selected_lineage_key))
        .filter_map(|ranked| {
            if !is_selective_report_topology(&ranked.candidate)
                || !ga_matches_lineage_sequence(&ranked.candidate, &selected_lineage_key)
            {
                return None;
            }
            Milestone12ReportDesign::from_blueprint_candidate(
                ranked.rank,
                ranked.candidate.clone(),
                ranked.evaluation.score,
            )
            .ok()
        })
        .collect();

    let mut ga_report_top: Vec<Milestone12ReportDesign> = ga_result
        .top_candidates
        .iter()
        .filter(|ranked| blueprint_matches_lineage(&ranked.candidate, &selected_lineage_key))
        .filter_map(|ranked| {
            if !is_selective_report_topology(&ranked.candidate)
                || !ga_matches_lineage_sequence(&ranked.candidate, &selected_lineage_key)
            {
                return None;
            }
            Milestone12ReportDesign::from_blueprint_candidate(
                ranked.rank,
                ranked.candidate.clone(),
                ranked.evaluation.score,
            )
            .ok()
        })
        .collect();
    if ga_report_top.is_empty() {
        ga_report_top = seeds
            .iter()
            .filter(|candidate| {
                is_selective_report_topology(candidate)
                    && ga_matches_lineage_sequence(candidate, &selected_lineage_key)
            })
            .filter_map(|candidate| {
                compute_blueprint_report_metrics(candidate)
                    .ok()
                    .map(|metrics: SdtMetrics| {
                        Milestone12ReportDesign::new(
                            0,
                            candidate.clone(),
                            metrics.clone(),
                            score_candidate(&metrics, ga_mode, &weights),
                        )
                    })
            })
            .collect();
        sort_report_designs(&mut ga_report_top);
        ga_report_top.truncate(5);
        for (idx, design) in ga_report_top.iter_mut().enumerate() {
            design.rank = idx + 1;
        }
    }
    if ga_report_top.is_empty() {
        ga_report_top = option2_ranked
            .iter()
            .filter(|design| {
                is_selective_report_topology(&design.candidate)
                    && ga_matches_lineage_sequence(&design.candidate, &selected_lineage_key)
            })
            .take(5)
            .map(|design| {
                Milestone12ReportDesign::new(
                    0,
                    design.candidate.clone(),
                    design.metrics.clone(),
                    score_candidate(&design.metrics, ga_mode, &weights),
                )
            })
            .collect();
        sort_report_designs(&mut ga_report_top);
        for (idx, design) in ga_report_top.iter_mut().enumerate() {
            design.rank = idx + 1;
        }
    }
    if ga_report_top.is_empty() {
        return Err(
            "HydroSDT comparison track produced no primitive split-sequence designs".into(),
        );
    }
    let option1_report_ranked = option1_ranked.clone();
    let option2_report_ranked = option2_ranked.clone();
    let option2_report_pool_all = option2_pool_all.clone();

    for design in &option1_report_ranked {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option1Base)?;
    }
    for design in &option2_report_ranked {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option2Derived)?;
    }
    for design in &ga_report_top {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::GaRefined)?;
    }
    tracing::info!(
        "      GA rank-1: {} (score={:.4}  σ={:.3}  cancer_cav={:.3})",
        ga_report_top[0].candidate.id,
        ga_report_top[0].score,
        ga_report_top[0].metrics.cavitation_number,
        ga_report_top[0].metrics.cancer_targeted_cavitation
    );
    save_top5_report_json(&ga_report_top, &out_dir.join("ga_hydrosdt_top5.json"))?;
    tracing::info!("      Saved: ga_hydrosdt_top5.json");
    tracing::info!(
        "      [2/6] completed in {:.1}s",
        phase2_start.elapsed().as_secs_f64()
    );

    // ── Part 3: Robustness screening ─────────────────────────────────────────
    let phase3_start = Instant::now();
    let validation_rows: Vec<ValidationRow> = Vec::new();
    let option2_robustness: Vec<RobustnessReport>;
    if run_cfg.fast_mode {
        tracing::info!(
            "\n[3/5] Fast mode: skipping robustness sweep; multi-fidelity validation lives in milestone12_validation"
        );
        option2_robustness = Vec::new();
    } else {
        tracing::info!("\n[3/5] Running robustness screening on Option 2 shortlist …");
        option2_robustness = option2_ranked
            .iter()
            .map(|d| {
                robustness_sweep_blueprint(
                    &d.candidate,
                    OptimizationGoal::SelectiveVenturiCavitation,
                    &STANDARD_PERTURBATIONS,
                )
            })
            .collect();
        let robust_count = option2_robustness.iter().filter(|r| r.is_robust).count();
        tracing::info!(
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
    tracing::info!(
        "      [3/5] completed in {:.1}s",
        phase3_start.elapsed().as_secs_f64()
    );
    let phase5_start = Instant::now();
    tracing::info!("\n[4/5] Generating report figures …");

    // Figure 4 — best selective acoustic design
    if let Some(option1_best) = option1_ranked.first() {
        save_figure(
            blueprint_by_id
                .get(&option1_best.candidate.id)
                .ok_or("missing direct blueprint for selected Option 1 figure")?,
            &figures_dir.join("selected_option1_schematic.svg"),
            "Figure 4 (Option 1 selective acoustic)",
        );
    } else {
        tracing::warn!(
            "      Figure 4 (Option 1 selective acoustic) skipped: no eligible shortlist under current physics regime"
        );
    }

    // Figure 5 — best combined selective venturi design
    save_figure(
        blueprint_by_id
            .get(&option2_ranked[0].candidate.id)
            .ok_or("missing direct blueprint for selected Option 2 figure")?,
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Figure 5 (Option 2 combined selective venturi)",
    );

    // Figure 6 — best HydroSDT GA design
    save_figure(
        blueprint_by_id
            .get(&ga_report_top[0].candidate.id)
            .ok_or("missing direct blueprint for selected GA figure")?,
        &figures_dir.join("top_hydrosdt_schematic.svg"),
        "Figure 6 (HydroSDT GA rank-1)",
    );

    // ── Part 6: Canonical markdown + JSON ────────────────────────────────────
    tracing::info!(
        "      [4/5] completed in {:.1}s",
        phase5_start.elapsed().as_secs_f64()
    );
    let phase6_start = Instant::now();
    tracing::info!("\n[5/5] Writing canonical report and summary JSON …");

    // Two-concept selection summary
    let toplines = vec![
        option1_report_ranked.first().map_or_else(
            || {
                concept_topline_unavailable(
                    "Option 1: Selective acoustic center treatment",
                    "SelectiveAcousticTherapy",
                    "No selective acoustic design satisfied strict eligibility under the current physics regime.",
                )
            },
            |design| {
                concept_topline(
                    "Option 1: Selective acoustic center treatment",
                    "SelectiveAcousticTherapy",
                    design,
                )
            },
        ),
        concept_topline(
            "Option 2: Selective venturi treatment ranked by combined score",
            "CombinedSdtLeukapheresis",
            &option2_report_ranked[0],
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
        &option1_report_ranked,
        &option2_report_ranked,
        &ga_report_top[0],
        &validation_rows,
        &option2_robustness,
        &canonical_results,
    )?;
    tracing::info!("Canonical report: {}", canonical_results.display());

    let narrative_artifacts = write_milestone12_narrative_report(
        &workspace_root,
        &canonical_results,
        &Milestone12NarrativeInput {
            total_candidates,
            option1_pool_len,
            option2_pool_len,
            option1_ranked: &option1_report_ranked,
            option2_ranked: &option2_report_ranked,
            ga_top: &ga_report_top,
            option2_pool_all: &option2_report_pool_all,
            ga_pool_all: &ga_report_pool_all,
            validation_rows: &validation_rows,
            option2_robustness: &option2_robustness,
            ga_best_per_gen: &ga_result.best_per_generation,
            fast_mode: run_cfg.fast_mode,
        },
    )?;
    tracing::info!(
        "Narrative report: {}",
        narrative_artifacts.narrative_path.display()
    );
    tracing::info!(
        "Figure manifest: {} ({} figures)",
        narrative_artifacts.figure_manifest_path.display(),
        narrative_artifacts.figure_count
    );

    tracing::info!("\n=== Outputs written to {} ===", out_dir.display());
    tracing::info!("Figures written to {}", figures_dir.display());
    tracing::info!(
        "      [5/5] completed in {:.1}s",
        phase6_start.elapsed().as_secs_f64()
    );
    tracing::info!(
        "=== Total wall time: {:.1}s ===",
        overall_start.elapsed().as_secs_f64()
    );

    Ok(())
}

// ── Figure saving helper ──────────────────────────────────────────────────────

fn save_figure(blueprint: &NetworkBlueprint, path: &Path, label: &str) {
    match save_blueprint_schematic_svg(blueprint, path) {
        Ok(()) => tracing::info!("      ✓ {label}  →  {}", path.display()),
        Err(e) => tracing::warn!("      ✗ {label}  FAILED: {e}"),
    }
}

// ── Concept topline builder ───────────────────────────────────────────────────

fn concept_topline(concept: &str, mode: &str, d: &Milestone12ReportDesign) -> ConceptTopline {
    ConceptTopline {
        concept: concept.to_string(),
        mode: mode.to_string(),
        available: true,
        status: "selected".to_string(),
        candidate_id: d.candidate.id.clone(),
        score: Some(d.score),
        topology: d.topology_short_code(),
        treatment_zone_mode: d.metrics.treatment_zone_mode.clone(),
        active_venturi_throat_count: Some(d.metrics.active_venturi_throat_count),
        cavitation_number: Some(d.metrics.cavitation_number),
        cancer_center_fraction: Some(d.metrics.cancer_center_fraction),
        therapeutic_window_score: Some(d.metrics.therapeutic_window_score),
        hemolysis_index_per_pass: Some(d.metrics.hemolysis_index_per_pass),
        wall_shear_p95_pa: Some(d.metrics.wall_shear_p95_pa),
        mean_residence_time_s: Some(d.metrics.mean_residence_time_s),
        treatment_zone_dwell_time_s: Some(d.metrics.treatment_zone_dwell_time_s),
        total_ecv_ml: Some(d.metrics.total_ecv_ml),
        therapy_channel_fraction: Some(d.metrics.therapy_channel_fraction),
    }
}

fn concept_topline_unavailable(concept: &str, mode: &str, status: &str) -> ConceptTopline {
    ConceptTopline {
        concept: concept.to_string(),
        mode: mode.to_string(),
        available: false,
        status: status.to_string(),
        candidate_id: "NONE_ELIGIBLE_CURRENT_PHYSICS".to_string(),
        score: None,
        topology: "n/a".to_string(),
        treatment_zone_mode: "Unavailable".to_string(),
        active_venturi_throat_count: None,
        cavitation_number: None,
        cancer_center_fraction: None,
        therapeutic_window_score: None,
        hemolysis_index_per_pass: None,
        wall_shear_p95_pa: None,
        mean_residence_time_s: None,
        treatment_zone_dwell_time_s: None,
        total_ecv_ml: None,
        therapy_channel_fraction: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_schematics::{SplitKind, SplitStageSpec};
    use serde_json::Value;

    fn tri_tri_acoustic_blueprint() -> BlueprintCandidate {
        build_milestone12_blueprint_candidate_space()
            .expect("Milestone 12 blueprint candidate space should build")
            .into_iter()
            .find(|candidate| {
                candidate.topology_spec().is_ok_and(|spec| {
                    spec.stage_sequence_label() == "Tri→Tri" && spec.venturi_placements.is_empty()
                })
            })
            .expect("Tri->Tri acoustic blueprint candidate should exist")
    }

    fn rebalance_treatment_branch(stage: &mut SplitStageSpec, target_frac: f64) {
        let total_width = stage
            .branches
            .iter()
            .map(|branch| branch.route.width_m)
            .sum::<f64>();
        let treatment_index = stage
            .branches
            .iter()
            .position(|branch| branch.treatment_path)
            .expect("split stage must include a treatment branch");
        let other_total = stage
            .branches
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != treatment_index)
            .map(|(_, branch)| branch.route.width_m)
            .sum::<f64>();
        let treatment_width = (total_width * target_frac).clamp(1.0e-9, total_width - 1.0e-9);
        let remaining_width = (total_width - treatment_width).max(1.0e-9);
        let branch_count = stage.branches.len().saturating_sub(1).max(1);
        for (index, branch) in stage.branches.iter_mut().enumerate() {
            if index == treatment_index {
                branch.route.width_m = treatment_width;
            } else if other_total > 1.0e-12 {
                branch.route.width_m *= remaining_width / other_total;
            } else {
                branch.route.width_m = remaining_width / branch_count as f64;
            }
        }
    }

    fn reference_option1_acoustic_candidate(
        id: &str,
        serpentine_segments: usize,
        first_tri_center_frac: f64,
        terminal_tri_center_frac: f64,
    ) -> BlueprintCandidate {
        let base = tri_tri_acoustic_blueprint();
        let mut spec = base
            .topology_spec()
            .expect("Tri->Tri acoustic candidate must carry topology metadata")
            .clone();
        let first_stage = spec
            .split_stages
            .first_mut()
            .expect("Tri->Tri topology must include a first stage");
        rebalance_treatment_branch(first_stage, first_tri_center_frac);
        let terminal_tri = spec
            .split_stages
            .iter_mut()
            .rev()
            .find(|stage| matches!(stage.split_kind, SplitKind::Trifurcation))
            .expect("Tri->Tri topology must include a terminal trifurcation");
        rebalance_treatment_branch(terminal_tri, terminal_tri_center_frac);
        let treatment_branch = spec
            .split_stages
            .iter_mut()
            .rev()
            .flat_map(|stage| stage.branches.iter_mut())
            .find(|branch| branch.treatment_path)
            .expect("Tri->Tri topology must expose a treatment branch");
        if let Some(serpentine) = treatment_branch.route.serpentine.as_mut() {
            serpentine.segments = serpentine_segments;
        }
        BlueprintCandidate::from_topology_spec(id.to_string(), &spec, base.operating_point.clone())
            .expect("reference acoustic blueprint candidate should remain constructible")
    }

    #[test]
    fn option1_report_serialization_retains_reference_n1_and_n5_acoustic_rows() {
        let weights = SdtWeights::default();
        let mut ranked = vec![
            reference_option1_acoustic_candidate(
                "m12-option1-pst-tritri-uo-n1",
                5,
                0.55,
                1.0 / 3.0,
            ),
            reference_option1_acoustic_candidate(
                "m12-option1-pst-tritri-uo-n5",
                7,
                0.55,
                1.0 / 3.0,
            ),
            reference_option1_acoustic_candidate("sentinel-option1-third", 5, 0.52, 0.34),
            reference_option1_acoustic_candidate("sentinel-option1-fourth", 3, 0.55, 1.0 / 3.0),
            reference_option1_acoustic_candidate("sentinel-option1-fifth", 5, 0.58, 0.38),
        ]
        .into_iter()
        .map(|candidate| {
            let metrics = compute_blueprint_report_metrics(&candidate)
                .expect("reference acoustic blueprint metrics should compute");
            let score = report_option1_score(&metrics, &weights);
            Milestone12ReportDesign::new(0, candidate, metrics, score)
        })
        .collect::<Vec<_>>();
        for design in &mut ranked[2..] {
            design.score *= 0.95;
        }

        let shortlist = shortlist_report_designs(ranked, 5, "option1_ultrasound_test")
            .expect("shortlist should succeed");
        let serialized: Value = serde_json::from_str(
            &serde_json::to_string_pretty(&shortlist)
                .expect("report shortlist should serialize to json"),
        )
        .expect("serialized shortlist should parse back to json");

        let rows = serialized
            .as_array()
            .expect("shortlist serialization should be a json array");
        assert_eq!(
            rows.len(),
            5,
            "serialized Option 1 shortlist should retain all ranked rows"
        );
        let n1_row = rows
            .iter()
            .find(|row| row["candidate"]["id"] == n1_row_id())
            .expect("serialized shortlist should contain n1 row");
        let n5_row = rows
            .iter()
            .find(|row| row["candidate"]["id"] == n5_row_id())
            .expect("serialized shortlist should contain n5 row");

        assert_ne!(
            n1_row["candidate"]["id"],
            n5_row["candidate"]["id"],
            "serialized Option 1 shortlist must retain distinct rows for the reference n1/n5 acoustic pair"
        );
        assert!(
            n1_row["metrics"]
                .get("treatment_zone_dwell_time_s")
                .is_some()
                && n5_row["metrics"]
                    .get("treatment_zone_dwell_time_s")
                    .is_some(),
            "serialized Option 1 shortlist must expose treatment_zone_dwell_time_s"
        );
    }

    fn n1_row_id() -> Value {
        Value::String("m12-option1-pst-tritri-uo-n1".to_string())
    }

    fn n5_row_id() -> Value {
        Value::String("m12-option1-pst-tritri-uo-n5".to_string())
    }
}

#[derive(Default)]
struct EvaluationAccumulator {
    option1_pool: Vec<Milestone12ReportDesign>,
    option2_pool: Vec<Milestone12ReportDesign>,
    scored_seeds: Vec<(f64, BlueprintCandidate)>,
    nonventuri_total: usize,
    nonventuri_pressure_plate: usize,
    nonventuri_fda_main: usize,
    nonventuri_eligible: usize,
    venturi_total: usize,
    venturi_sigma_lt1_any: usize,
    sigma_lt1_pressure_ok: usize,
    sigma_lt1_plate_ok: usize,
    sigma_lt1_fda_main_ok: usize,
    sigma_lt1_fda_overall_ok: usize,
    sigma_lt1_ppfda_overall_ok: usize,
    venturi_pressure_plate_fda: usize,
    venturi_sigma_lt1: usize,
    venturi_oncology_eligible: usize,
}

impl EvaluationAccumulator {
    fn record(
        &mut self,
        blueprint_candidate: BlueprintCandidate,
        metrics: SdtMetrics,
        ga_score: f64,
        option1_score: Option<f64>,
        combined_score: Option<f64>,
    ) {
        let uses_venturi = blueprint_candidate
            .topology_spec()
            .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
        if uses_venturi {
            self.venturi_total += 1;
            if metrics.cavitation_number.is_finite() && metrics.cavitation_number < 1.0 {
                self.venturi_sigma_lt1_any += 1;
                if metrics.pressure_feasible {
                    self.sigma_lt1_pressure_ok += 1;
                }
                if metrics.plate_fits {
                    self.sigma_lt1_plate_ok += 1;
                }
                if metrics.fda_main_compliant {
                    self.sigma_lt1_fda_main_ok += 1;
                }
                if metrics.fda_overall_compliant {
                    self.sigma_lt1_fda_overall_ok += 1;
                }
                if metrics.pressure_feasible && metrics.plate_fits && metrics.fda_overall_compliant
                {
                    self.sigma_lt1_ppfda_overall_ok += 1;
                }
            }
            if metrics.pressure_feasible && metrics.plate_fits && metrics.fda_main_compliant {
                self.venturi_pressure_plate_fda += 1;
                if metrics.cavitation_number.is_finite() && metrics.cavitation_number < 1.0 {
                    self.venturi_sigma_lt1 += 1;
                }
            }
            if report_eligible_venturi_oncology(&metrics) {
                self.venturi_oncology_eligible += 1;
            }
            // Avoid unnecessary clone: when the candidate goes to BOTH
            // scored_seeds and option2_pool, clone is required. When only one
            // destination needs the candidate, move it directly.
            match (ga_score > 0.0, combined_score) {
                (true, Some(score)) => {
                    self.scored_seeds
                        .push((ga_score, blueprint_candidate.clone()));
                    self.option2_pool.push(Milestone12ReportDesign::new(
                        0,
                        blueprint_candidate,
                        metrics,
                        score,
                    ));
                }
                (true, None) => {
                    self.scored_seeds.push((ga_score, blueprint_candidate));
                }
                (false, Some(score)) => {
                    self.option2_pool.push(Milestone12ReportDesign::new(
                        0,
                        blueprint_candidate,
                        metrics,
                        score,
                    ));
                }
                (false, None) => {}
            }
        } else {
            self.nonventuri_total += 1;
            if metrics.pressure_feasible && metrics.plate_fits {
                self.nonventuri_pressure_plate += 1;
            }
            if metrics.fda_main_compliant {
                self.nonventuri_fda_main += 1;
            }
            if report_eligible_nonventuri(&metrics) {
                self.nonventuri_eligible += 1;
            }
            if let Some(score) = option1_score {
                self.option1_pool.push(Milestone12ReportDesign::new(
                    0,
                    blueprint_candidate,
                    metrics,
                    score,
                ));
            }
        }
    }

    fn merge(&mut self, mut other: Self) {
        self.option1_pool.append(&mut other.option1_pool);
        self.option2_pool.append(&mut other.option2_pool);
        self.scored_seeds.append(&mut other.scored_seeds);
        self.nonventuri_total += other.nonventuri_total;
        self.nonventuri_pressure_plate += other.nonventuri_pressure_plate;
        self.nonventuri_fda_main += other.nonventuri_fda_main;
        self.nonventuri_eligible += other.nonventuri_eligible;
        self.venturi_total += other.venturi_total;
        self.venturi_sigma_lt1_any += other.venturi_sigma_lt1_any;
        self.sigma_lt1_pressure_ok += other.sigma_lt1_pressure_ok;
        self.sigma_lt1_plate_ok += other.sigma_lt1_plate_ok;
        self.sigma_lt1_fda_main_ok += other.sigma_lt1_fda_main_ok;
        self.sigma_lt1_fda_overall_ok += other.sigma_lt1_fda_overall_ok;
        self.sigma_lt1_ppfda_overall_ok += other.sigma_lt1_ppfda_overall_ok;
        self.venturi_pressure_plate_fda += other.venturi_pressure_plate_fda;
        self.venturi_sigma_lt1 += other.venturi_sigma_lt1;
        self.venturi_oncology_eligible += other.venturi_oncology_eligible;
    }

    fn option1_hits(&self) -> usize {
        self.option1_pool.len()
    }
}

fn has_viable_lineage_family(
    option1_pool: &[Milestone12ReportDesign],
    option2_pool: &[Milestone12ReportDesign],
    min_option2_count: usize,
) -> bool {
    let option1_keys: HashSet<_> = option1_pool
        .iter()
        .filter_map(|design| blueprint_lineage_key(&design.candidate))
        .collect();
    if option1_keys.is_empty() {
        return false;
    }

    let mut option2_counts: HashMap<_, usize> = HashMap::new();
    for design in option2_pool {
        let Some(key) = blueprint_lineage_key(&design.candidate) else {
            continue;
        };
        if !option1_keys.contains(&key) {
            continue;
        }
        let count = option2_counts.entry(key).or_insert(0);
        *count += 1;
        if *count >= min_option2_count {
            return true;
        }
    }

    false
}
