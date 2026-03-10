//! Milestone 12 — Option 1: Selective Acoustic Residence Separation.
//!
//! Evaluates non-venturi selective-routing topologies (e.g. Tri→Tri, Tri→Tri→Tri)
//! and ranks them by `SelectiveAcousticTherapy` score.  Produces:
//!
//! - `report/milestone12/two_concept_option1_ultrasound_top5.json`
//! - `report/figures/selected_option1_schematic.svg`
//! - Regenerates the narrative report (using Option 2 / GA data from prior runs
//!   if available on disk).
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_option1 --no-default-features
//! # fast iteration
//! M12_FAST=1 cargo run -p cfd-optim --example milestone12_option1 --no-default-features
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, compute_blueprint_report_metrics,
    is_milestone12_lineage_topology, load_top5_report_json, save_blueprint_schematic_svg,
    save_json_pretty, save_top5_report_json, score_candidate, shortlist_report_designs,
    write_milestone12_narrative_report, write_milestone12_results, BlueprintCandidate,
    Milestone12NarrativeInput, Milestone12ReportDesign, Milestone12Stage, OptimMode, SdtMetrics,
    SdtWeights,
};
use rayon::prelude::*;
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

// ── Helpers ──────────────────────────────────────────────────────────────────

fn report_eligible_nonventuri(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.therapy_channel_fraction > 0.0
}

fn report_option1_score(metrics: &SdtMetrics, weights: &SdtWeights) -> f64 {
    let base = score_candidate(metrics, OptimMode::SelectiveAcousticTherapy, weights);
    let dwell_tie_break = metrics.treatment_zone_dwell_time_s.max(0.0) * 1.0e-6;
    let pressure_tie_break = metrics.total_pressure_drop_pa.max(0.0) * 1.0e-12;
    base + dwell_tie_break - pressure_tie_break
}

fn is_selective_report_topology(candidate: &BlueprintCandidate) -> bool {
    is_milestone12_lineage_topology(candidate)
}

fn save_figure(blueprint: &cfd_schematics::NetworkBlueprint, path: &Path, label: &str) {
    match save_blueprint_schematic_svg(blueprint, path) {
        Ok(()) => tracing::info!("      \u{2713} {label}  \u{2192}  {}", path.display()),
        Err(e) => tracing::warn!("      \u{2717} {label}  FAILED: {e}"),
    }
}

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
            "      {label}: starting {} candidate evaluations (heartbeat every {})",
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

fn fast_mode() -> bool {
    std::env::var("M12_FAST")
        .ok()
        .map(|v| {
            let s = v.trim().to_ascii_lowercase();
            s == "1" || s == "true" || s == "yes"
        })
        .unwrap_or(false)
}

fn fast_env(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
        .max(1)
}

// ── Main ─────────────────────────────────────────────────────────────────────

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

    tracing::info!("=== Milestone 12 — Option 1: Selective Acoustic ===\n");
    let overall_start = Instant::now();
    let weights = SdtWeights::default();
    let is_fast = fast_mode();

    // ── Build candidate space ────────────────────────────────────────────────
    let t_build = Instant::now();
    let raw_candidates = build_milestone12_blueprint_candidate_space()?;
    tracing::info!(
        "Candidate space: {} blueprints ({:.1}s)",
        raw_candidates.len(),
        t_build.elapsed().as_secs_f64()
    );

    // ── Filter to non-venturi selective topologies ───────────────────────────
    let acoustic_candidates: Vec<BlueprintCandidate> = raw_candidates
        .into_iter()
        .filter(|c| {
            !c.topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi)
                && is_selective_report_topology(c)
        })
        .collect();

    let eval_cap = if is_fast {
        fast_env("M12_FAST_ACOUSTIC_EVAL_MAX", 200)
    } else {
        acoustic_candidates.len()
    };
    let to_evaluate = acoustic_candidates.len().min(eval_cap);
    tracing::info!(
        "Non-venturi selective candidates: {} (evaluating {})",
        acoustic_candidates.len(),
        to_evaluate
    );

    // ── Evaluate ─────────────────────────────────────────────────────────────
    let progress = Arc::new(ScanProgress::new("option1 scan", to_evaluate));
    let progress_ref = Arc::clone(&progress);

    let option1_pool: Vec<Milestone12ReportDesign> = acoustic_candidates
        .into_par_iter()
        .take(to_evaluate)
        .filter_map(move |candidate| {
            let metrics = match compute_blueprint_report_metrics(&candidate) {
                Ok(m) => m,
                Err(_) => {
                    progress_ref.record();
                    return None;
                }
            };
            progress_ref.record();
            if !report_eligible_nonventuri(&metrics) {
                return None;
            }
            let score = report_option1_score(&metrics, &weights);
            Some(Milestone12ReportDesign::new(0, candidate, metrics, score))
        })
        .collect();
    progress.finish();

    tracing::info!("Option 1 eligible: {}", option1_pool.len());

    // ── Rank and shortlist ───────────────────────────────────────────────────
    let option1_pool_len = option1_pool.len();
    let shortlist_len = option1_pool.len().min(5);
    let mut option1_ranked =
        shortlist_report_designs(option1_pool, shortlist_len, "option1_ultrasound")?;
    for (idx, design) in option1_ranked.iter_mut().enumerate() {
        design.rank = idx + 1;
    }

    if let Some(best) = option1_ranked.first() {
        tracing::info!(
            "Option 1 rank-1: {} (score={:.4})",
            best.topology_short_code(),
            best.score
        );
    } else {
        tracing::warn!("No eligible Option 1 candidates under current physics regime");
    }

    // ── Validate ─────────────────────────────────────────────────────────────
    for design in &option1_ranked {
        cfd_optim::validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option1Base)?;
    }

    // ── Save track outputs ───────────────────────────────────────────────────
    save_top5_report_json(
        &option1_ranked,
        &out_dir.join("two_concept_option1_ultrasound_top5.json"),
    )?;
    tracing::info!("Saved: two_concept_option1_ultrasound_top5.json");

    if let Some(topline) = option1_ranked.first().map(|d| {
        concept_topline(
            "Option 1: Selective acoustic center treatment",
            "SelectiveAcousticTherapy",
            d,
        )
    }) {
        save_json_pretty(
            &vec![topline],
            &out_dir.join("option1_concept_topline.json"),
        )?;
    }

    // ── Figure 4 ─────────────────────────────────────────────────────────────
    if let Some(best) = option1_ranked.first() {
        save_figure(
            best.candidate.blueprint(),
            &figures_dir.join("selected_option1_schematic.svg"),
            "Figure 4 (Option 1 selective acoustic)",
        );
    }

    // ── Assemble report (load other tracks if available) ─────────────────────
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;
    let ga_top = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;

    if !option2_ranked.is_empty() && !ga_top.is_empty() {
        let canonical = workspace_root.join("report").join("milestone12_results.md");
        write_milestone12_results(
            0, // total_candidates not tracked in this partial run
            option1_pool_len,
            option2_ranked.len(),
            &option1_ranked,
            &option2_ranked,
            &ga_top[0],
            &[],
            &[],
            &canonical,
        )?;
        tracing::info!("Updated canonical results: {}", canonical.display());

        let narrative_artifacts = write_milestone12_narrative_report(
            &workspace_root,
            &canonical,
            &Milestone12NarrativeInput {
                total_candidates: 0,
                option1_pool_len,
                option2_pool_len: option2_ranked.len(),
                option1_ranked: &option1_ranked,
                option2_ranked: &option2_ranked,
                ga_top: &ga_top,
                option2_pool_all: &option2_ranked,
                ga_pool_all: &ga_top,
                validation_rows: &[],
                option2_robustness: &[],
                ga_best_per_gen: &[],
                fast_mode: is_fast,
            },
        )?;
        tracing::info!(
            "Narrative report: {}",
            narrative_artifacts.narrative_path.display()
        );
    } else {
        tracing::info!(
            "Skipping report assembly — run milestone12_option2 and milestone12_ga first"
        );
    }

    tracing::info!(
        "\n=== Option 1 complete in {:.1}s ===",
        overall_start.elapsed().as_secs_f64()
    );

    Ok(())
}
