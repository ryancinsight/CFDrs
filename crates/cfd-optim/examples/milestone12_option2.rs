//! Milestone 12 — Option 2: Selective Venturi Cavitation.
//!
//! Evaluates venturi-equipped selective-routing topologies, ranks them by
//! `CombinedSdtLeukapheresis` score, and applies robustness screening.
//!
//! The lineage constraint (shared topology family with Option 1) is enforced
//! when Option 1 results exist on disk; otherwise the best venturi family is
//! selected unconstrainedly.
//!
//! Produces:
//! - `report/milestone12/two_concept_option2_venturi_top5.json`
//! - `report/milestone12/option2_combined_robustness_top5.json`
//! - `report/milestone12/option2_scored_seeds.json`  (for GA warm-start)
//! - `report/figures/selected_option2_combined_schematic.svg`
//! - Regenerates the narrative report.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_option2 --no-default-features
//! M12_FAST=1 cargo run -p cfd-optim --example milestone12_option2 --no-default-features
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, compute_blueprint_report_metrics,
    evaluate_blueprint_candidate, evaluate_blueprint_genetic_refinement,
    is_milestone12_lineage_topology, load_top5_report_json, milestone12_lineage_key,
    robustness_sweep_blueprint, save_blueprint_schematic_svg, save_json_pretty,
    save_top5_report_json, score_candidate, sort_report_designs,
    validate_milestone12_candidate, write_milestone12_narrative_report, write_milestone12_results,
    BlueprintCandidate, Milestone12LineageKey, Milestone12NarrativeInput, Milestone12ReportDesign,
    Milestone12Stage, OptimMode, OptimizationGoal, RobustnessReport, SdtMetrics, SdtWeights,
    STANDARD_PERTURBATIONS,
};
use rayon::prelude::*;
use std::collections::HashMap;
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
    OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 70.0,
    }
}

fn is_selective_report_topology(candidate: &BlueprintCandidate) -> bool {
    is_milestone12_lineage_topology(candidate)
}

fn blueprint_lineage_key(candidate: &BlueprintCandidate) -> Option<Milestone12LineageKey> {
    milestone12_lineage_key(candidate)
}

fn save_figure(blueprint: &cfd_schematics::NetworkBlueprint, path: &Path, label: &str) {
    match save_blueprint_schematic_svg(blueprint, path) {
        Ok(()) => tracing::info!("      \u{2713} {label}  \u{2192}  {}", path.display()),
        Err(e) => tracing::warn!("      \u{2717} {label}  FAILED: {e}"),
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
            "      {label}: starting {} evaluations (heartbeat every {})",
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

    tracing::info!("=== Milestone 12 — Option 2: Selective Venturi Cavitation ===\n");
    let overall_start = Instant::now();
    let weights = SdtWeights::default();
    let oncology_mode = option2_mode();
    let is_fast = fast_mode();

    // ── Build candidate space ────────────────────────────────────────────────
    let t_build = Instant::now();
    let all_candidates = build_milestone12_blueprint_candidate_space()?;
    let total_candidates = all_candidates.len();
    tracing::info!(
        "Candidate space: {} blueprints ({:.1}s)",
        total_candidates,
        t_build.elapsed().as_secs_f64()
    );

    // ── Load Option 1 results (from prior run) for lineage matching ─────────
    let option1_from_disk =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let have_option1 = !option1_from_disk.is_empty();
    if have_option1 {
        tracing::info!(
            "Loaded {} Option 1 designs from disk for lineage matching",
            option1_from_disk.len()
        );
    } else {
        tracing::info!("No Option 1 results on disk — will evaluate acoustic candidates inline");
    }

    // ── Evaluate candidates ──────────────────────────────────────────────────
    // We need both venturi candidates (for Option 2) and acoustic candidates
    // (for lineage matching if Option 1 wasn't run separately).
    let selective_candidates: Vec<BlueprintCandidate> = all_candidates
        .into_iter()
        .filter(|c| is_selective_report_topology(c))
        .collect();

    let progress = Arc::new(ScanProgress::new("option2 scan", selective_candidates.len()));
    let progress_ref = Arc::clone(&progress);
    let have_option1_for_par = have_option1;

    struct CandidateResult {
        candidate: BlueprintCandidate,
        metrics: SdtMetrics,
        is_venturi: bool,
        ga_score: f64,
        option1_score: Option<f64>,
        option2_score: Option<f64>,
    }

    let results: Vec<CandidateResult> = selective_candidates
        .into_par_iter()
        .filter_map(move |candidate| {
            let _eval = match evaluate_blueprint_candidate(&candidate) {
                Ok(e) => e,
                Err(_) => {
                    progress_ref.record();
                    return None;
                }
            };
            let metrics = match compute_blueprint_report_metrics(&candidate) {
                Ok(m) => m,
                Err(_) => {
                    progress_ref.record();
                    return None;
                }
            };
            progress_ref.record();

            let is_venturi = candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);

            let ga_score = evaluate_blueprint_genetic_refinement(&candidate, _eval)
                .map_or(0.0, |e| e.score);

            let option1_score = if !is_venturi
                && !have_option1_for_par
                && report_eligible_nonventuri(&metrics)
            {
                Some(report_option1_score(&metrics, &weights))
            } else {
                None
            };

            let option2_score =
                if is_venturi && report_eligible_venturi_oncology(&metrics) {
                    Some(score_candidate(&metrics, oncology_mode, &weights))
                } else {
                    None
                };

            Some(CandidateResult {
                candidate,
                metrics,
                is_venturi,
                ga_score,
                option1_score,
                option2_score,
            })
        })
        .collect();
    progress.finish();

    // ── Split into pools ─────────────────────────────────────────────────────
    let mut option1_pool: Vec<Milestone12ReportDesign> = Vec::new();
    let mut option2_pool: Vec<Milestone12ReportDesign> = Vec::new();
    let mut scored_seeds: Vec<(f64, BlueprintCandidate)> = Vec::new();

    for r in results {
        if r.is_venturi {
            if r.ga_score > 0.0 {
                scored_seeds.push((r.ga_score, r.candidate.clone()));
            }
            if let Some(score) = r.option2_score {
                option2_pool.push(Milestone12ReportDesign::new(0, r.candidate, r.metrics, score));
            }
        } else if let Some(score) = r.option1_score {
            option1_pool.push(Milestone12ReportDesign::new(0, r.candidate, r.metrics, score));
        }
    }

    tracing::info!(
        "Option 2 eligible: {}  (acoustic inline: {})",
        option2_pool.len(),
        option1_pool.len()
    );

    // ── Lineage matching ─────────────────────────────────────────────────────
    // Use loaded Option 1 results if available, otherwise use inline-evaluated acoustic pool.
    let option1_for_lineage: &[Milestone12ReportDesign] = if have_option1 {
        &option1_from_disk
    } else {
        &option1_pool
    };

    let mut option1_by_lineage: HashMap<Milestone12LineageKey, Vec<Milestone12ReportDesign>> =
        HashMap::new();
    let mut option2_by_lineage: HashMap<Milestone12LineageKey, Vec<Milestone12ReportDesign>> =
        HashMap::new();
    for design in option1_for_lineage {
        if let Some(key) = blueprint_lineage_key(&design.candidate) {
            option1_by_lineage.entry(key).or_default().push(design.clone());
        }
    }
    for design in &option2_pool {
        if let Some(key) = blueprint_lineage_key(&design.candidate) {
            option2_by_lineage
                .entry(key)
                .or_default()
                .push(design.clone());
        }
    }

    let mut viable_lineages: Vec<_> = option2_by_lineage
        .into_iter()
        .filter_map(|(key, mut derived)| {
            let option1_designs = option1_by_lineage.get(&key)?;
            let mut acoustic = option1_designs.clone();
            sort_report_designs(&mut acoustic);
            sort_report_designs(&mut derived);
            Some((key, acoustic[0].clone(), derived))
        })
        .collect();

    if viable_lineages.is_empty() {
        // Fall back: use all option2 candidates without lineage constraint.
        tracing::warn!(
            "No shared lineage between Option 1 and Option 2 — using best venturi family"
        );
        sort_report_designs(&mut option2_pool);
        option2_pool.truncate(5);
        for (idx, d) in option2_pool.iter_mut().enumerate() {
            d.rank = idx + 1;
        }
    } else {
        viable_lineages.sort_by(|a, b| {
            b.2[0]
                .score
                .total_cmp(&a.2[0].score)
                .then_with(|| a.1.candidate.id.cmp(&b.1.candidate.id))
        });

        let (_lineage_key, _option1_best, lineage_pool) = viable_lineages
            .into_iter()
            .find(|(_, _, derived)| derived.len() >= 5)
            .or_else(|| {
                // Relax the >=5 constraint if no lineage has enough.
                None
            })
            .ok_or("No shared Option 1→Option 2 family retained five derived venturi candidates")?;

        option2_pool = lineage_pool;
        option2_pool.truncate(5);
        for (idx, d) in option2_pool.iter_mut().enumerate() {
            d.rank = idx + 1;
        }
    }

    let option2_pool_len = option2_pool.len();
    if option2_pool.is_empty() {
        return Err("No eligible Option 2 candidates".into());
    }

    tracing::info!(
        "Option 2 rank-1: {} (score={:.4}  \u{03c3}={:.3})",
        option2_pool[0].topology_short_code(),
        option2_pool[0].score,
        option2_pool[0].metrics.cavitation_number
    );

    // ── Validate ─────────────────────────────────────────────────────────────
    for design in &option2_pool {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option2Derived)?;
    }

    // ── Save track outputs ───────────────────────────────────────────────────
    save_top5_report_json(
        &option2_pool,
        &out_dir.join("two_concept_option2_venturi_top5.json"),
    )?;
    tracing::info!("Saved: two_concept_option2_venturi_top5.json");

    // Save scored seeds for GA warm-start.
    scored_seeds.sort_by(|a, b| b.0.total_cmp(&a.0));
    let seed_candidates: Vec<&BlueprintCandidate> =
        scored_seeds.iter().take(150).map(|(_, c)| c).collect();
    save_json_pretty(
        &seed_candidates,
        &out_dir.join("option2_scored_seeds.json"),
    )?;
    tracing::info!(
        "Saved: option2_scored_seeds.json ({} seeds)",
        seed_candidates.len()
    );

    // ── Robustness screening ─────────────────────────────────────────────────
    let option2_robustness: Vec<RobustnessReport> = if is_fast {
        tracing::info!("Fast mode: skipping robustness sweep");
        Vec::new()
    } else {
        tracing::info!("Running robustness screening on Option 2 shortlist ...");
        let robustness: Vec<RobustnessReport> = option2_pool
            .iter()
            .map(|d| {
                robustness_sweep_blueprint(
                    &d.candidate,
                    OptimizationGoal::SelectiveVenturiCavitation,
                    &STANDARD_PERTURBATIONS,
                )
            })
            .collect();
        let robust_count = robustness.iter().filter(|r| r.is_robust).count();
        tracing::info!("Robust: {}/{}", robust_count, robustness.len());
        let json = serde_json::to_string_pretty(&robustness)?;
        std::fs::write(
            out_dir.join("option2_combined_robustness_top5.json"),
            json,
        )?;
        robustness
    };

    // ── Figure 5 ─────────────────────────────────────────────────────────────
    save_figure(
        option2_pool[0].candidate.blueprint(),
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Figure 5 (Option 2 combined selective venturi)",
    );

    // ── Assemble report ──────────────────────────────────────────────────────
    let option1_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let ga_top = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;

    if !ga_top.is_empty() {
        let canonical = workspace_root.join("report").join("milestone12_results.md");
        write_milestone12_results(
            total_candidates,
            option1_ranked.len(),
            option2_pool_len,
            &option1_ranked,
            &option2_pool,
            &ga_top[0],
            &[],
            &option2_robustness,
            &canonical,
        )?;
        tracing::info!("Updated canonical results: {}", canonical.display());

        let narrative_artifacts = write_milestone12_narrative_report(
            &workspace_root,
            &canonical,
            &Milestone12NarrativeInput {
                total_candidates,
                option1_pool_len: option1_ranked.len(),
                option2_pool_len,
                option1_ranked: &option1_ranked,
                option2_ranked: &option2_pool,
                ga_top: &ga_top,
                option2_pool_all: &option2_pool,
                ga_pool_all: &ga_top,
                validation_rows: &[],
                option2_robustness: &option2_robustness,
                ga_best_per_gen: &[],
                fast_mode: is_fast,
            },
        )?;
        tracing::info!(
            "Narrative report: {}",
            narrative_artifacts.narrative_path.display()
        );
    } else {
        tracing::info!("Skipping full report assembly — run milestone12_ga first for GA track");
    }

    tracing::info!(
        "\n=== Option 2 complete in {:.1}s ===",
        overall_start.elapsed().as_secs_f64()
    );

    Ok(())
}
