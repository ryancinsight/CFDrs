use crate::application::orchestration::{
    ensure_release_reports, fast_env, fast_mode, ga_matches_lineage_sequence, init_tracing,
    is_selective_report_topology, resolve_output_directories, save_figure,
};
use crate::application::search::genetic::{
    promote_option1_candidate_to_ga_seed, BlueprintGeneticOptimizer,
};
use crate::delivery::{load_top5_report_json, save_pareto_points, save_top5_report_json};
use crate::domain::OptimizationGoal;
use crate::reporting::{
    milestone12_lineage_key, validate_milestone12_candidate, write_goal_audit_report,
    GoalAuditEntry, GoalAuditStatus, Milestone12LineageKey, Milestone12ReportDesign,
    Milestone12Stage, ParetoPoint, ParetoTag,
};
use cfd_schematics::{BlueprintTopologyFactory, VenturiPlacementMode};

use super::report::{write_stage_summary, Milestone12GaSummary, GA_SUMMARY_PATH};
use super::types::{Milestone12GaRun, Milestone12StageArtifact};

pub fn run_milestone12_ga() -> Result<Milestone12GaRun, Box<dyn std::error::Error>> {
    init_tracing();
    ensure_release_reports()?;
    let (_, out_dir, figures_dir) = resolve_output_directories()?;
    let is_fast = fast_mode();
    let ga_population = if is_fast {
        fast_env("M12_FAST_GA_POPULATION", 25)
    } else {
        80
    };
    let ga_generations = if is_fast {
        fast_env("M12_FAST_GA_GENERATIONS", 8)
    } else {
        120
    };

    // Load seeds from previous-stage disk artifacts instead of rebuilding
    // the entire 80K candidate space (~1.2 GB). The pipeline guarantees
    // Option 1 and Option 2 run before GA.
    let option1_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;

    let mut seeds = Vec::with_capacity(2);
    if let Some(opt1_best) = option1_ranked.first() {
        if let Ok(promoted) = promote_option1_candidate_to_ga_seed(&opt1_best.candidate) {
            seeds.push(promoted);
        }
    }
    if let Some(opt2_best) = option2_ranked.first() {
        let dean_seed = opt2_best
            .candidate
            .topology_spec()
            .ok()
            .map(|spec| {
                let mut spec = spec.clone();
                for placement in &mut spec.venturi_placements {
                    placement.placement_mode = VenturiPlacementMode::CurvaturePeakDeanNumber;
                }
                spec
            })
            .and_then(|spec| BlueprintTopologyFactory::build(&spec).ok())
            .map(|blueprint| {
                crate::domain::BlueprintCandidate::new(
                    format!("{}-ga-dean-seed", opt2_best.candidate.id),
                    blueprint,
                    opt2_best.candidate.operating_point.clone(),
                )
            })
            .unwrap_or_else(|| opt2_best.candidate.clone());
        seeds.push(dean_seed);
    }
    if seeds.is_empty() {
        return Err("GA seed selection: no Option 1 or Option 2 ranked designs on disk".into());
    }
    drop(option1_ranked); // free the full designs
    let selected_lineage_key: Option<Milestone12LineageKey> = option2_ranked
        .first()
        .and_then(|design| milestone12_lineage_key(&design.candidate));

    let ga_result =
        BlueprintGeneticOptimizer::new(OptimizationGoal::InPlaceDeanSerpentineRefinement)
            .with_seeds(seeds)
            .with_population(ga_population)
            .with_max_generations(ga_generations)
            .with_top_k(5)
            .run()?;

    let mut ranked: Vec<Milestone12ReportDesign> = ga_result
        .top_candidates
        .iter()
        .filter(|ranked_candidate| {
            if let Some(ref lineage_key) = selected_lineage_key {
                is_selective_report_topology(&ranked_candidate.candidate)
                    && ga_matches_lineage_sequence(&ranked_candidate.candidate, lineage_key)
            } else {
                is_selective_report_topology(&ranked_candidate.candidate)
            }
        })
        .filter_map(|ranked_candidate| {
            let score = ranked_candidate.evaluation.score?;
            Milestone12ReportDesign::from_blueprint_candidate(
                ranked_candidate.rank,
                ranked_candidate.candidate.clone(),
                score,
            )
            .ok()
        })
        .collect();

    // Fallback 1: drop lineage filter, use all top_candidates.
    if ranked.is_empty() {
        ranked = ga_result
            .top_candidates
            .iter()
            .filter_map(|ranked_candidate| {
                let score = ranked_candidate.evaluation.score?;
                Milestone12ReportDesign::from_blueprint_candidate(
                    ranked_candidate.rank,
                    ranked_candidate.candidate.clone(),
                    score,
                )
                .ok()
            })
            .collect();
    }

    // Skip GA candidates that fail blueprint validation (e.g. channel
    // crossings from deeply nested mutations) instead of aborting the run.
    ranked.retain(|design| {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::GaRefined).is_ok()
    });

    // Fallback 2: if all top_candidates failed validation, search the
    // full archive (up to 200 candidates) for valid ones.
    if ranked.is_empty() {
        ranked = ga_result
            .all_candidates
            .iter()
            .filter_map(|ranked_candidate| {
                let score = ranked_candidate.evaluation.score?;
                let design = Milestone12ReportDesign::from_blueprint_candidate(
                    ranked_candidate.rank,
                    ranked_candidate.candidate.clone(),
                    score,
                )
                .ok()?;
                validate_milestone12_candidate(&design.candidate, Milestone12Stage::GaRefined)
                    .ok()?;
                Some(design)
            })
            .take(5)
            .collect();
    }
    for (index, design) in ranked.iter_mut().enumerate() {
        design.rank = index + 1;
    }
    if ranked.is_empty() {
        return Err("GA produced no viable designs after validation".into());
    }

    let audit_entries: Vec<GoalAuditEntry> = ga_result
        .all_candidates
        .iter()
        .map(|ranked_candidate| GoalAuditEntry {
            candidate_id: ranked_candidate.candidate.id.clone(),
            blueprint_name: ranked_candidate.candidate.blueprint.name.clone(),
            goal: OptimizationGoal::InPlaceDeanSerpentineRefinement,
            status: match ranked_candidate.evaluation.status {
                crate::application::objectives::BlueprintEvaluationStatus::Eligible => {
                    GoalAuditStatus::Eligible
                }
                crate::application::objectives::BlueprintEvaluationStatus::ScreenedOut => {
                    GoalAuditStatus::ScreenedOut
                }
            },
            reasons: ranked_candidate.evaluation.screening_reasons.clone(),
            score: ranked_candidate.evaluation.score,
            baseline_scores: ranked_candidate.evaluation.baseline_scores.clone(),
            exceeds_all_baselines: ranked_candidate.evaluation.exceeds_all_baselines,
            treatment_residence_time_s: Some(
                ranked_candidate
                    .evaluation
                    .residence
                    .treatment_residence_time_s,
            ),
            separation_efficiency: Some(
                ranked_candidate.evaluation.separation.separation_efficiency,
            ),
            cavitation_selectivity_score: Some(
                ranked_candidate.evaluation.venturi.cavitation_selectivity_score,
            ),
            peak_dean_number: Some(
                ranked_candidate
                    .evaluation
                    .venturi
                    .placements
                    .iter()
                    .map(|placement| placement.dean_number)
                    .fold(0.0_f64, f64::max),
            ),
        })
        .collect();
    let audit = write_goal_audit_report(&out_dir, "ga_audit", &audit_entries)?;

    // Persist lightweight Pareto points for GA pool (~32 bytes each).
    // Extract TTCI/RBC-protection directly from evaluation data — avoids
    // cloning 200 full BlueprintCandidates and computing SdtMetrics for each.
    let ga_pareto_points: Vec<ParetoPoint> = ga_result
        .all_candidates
        .iter()
        .take(200)
        .map(|rc| {
            let eval = &rc.evaluation;
            let cav_potential = eval
                .venturi
                .placements
                .iter()
                .map(|p| (1.0 - p.cavitation_number).clamp(0.0, 1.0))
                .fold(0.0_f64, f64::max);
            let constriction = eval.venturi.cavitation_selectivity_score.clamp(0.0, 1.0);
            let cavitation_intensity = cav_potential * (0.5 + 0.5 * constriction);
            ParetoPoint {
                cancer_targeted_cavitation: eval.separation.cancer_center_fraction.clamp(0.0, 1.0)
                    * cavitation_intensity,
                rbc_venturi_protection: eval
                    .separation
                    .rbc_peripheral_fraction
                    .clamp(0.0, 1.0)
                    * (1.0
                        - cavitation_intensity
                            * eval.venturi.rbc_exposure_fraction.clamp(0.0, 1.0)),
                score: eval.score_or_zero(),
                tag: ParetoTag::Ga,
            }
        })
        .collect();
    let ga_pool_all_path = out_dir.join("ga_pool_all.json");
    save_pareto_points(&ga_pareto_points, &ga_pool_all_path)?;
    drop(ga_pareto_points);

    let top5_path = out_dir.join("ga_hydrosdt_top5.json");
    save_top5_report_json(&ranked, &top5_path)?;
    let figure_path = figures_dir.join("top_hydrosdt_schematic.svg");
    save_figure(
        ranked[0].candidate.blueprint(),
        &figure_path,
        "Figure 6 (HydroSDT GA rank-1)",
    );
    write_stage_summary(
        &out_dir,
        GA_SUMMARY_PATH,
        &Milestone12GaSummary {
            authoritative_run: !crate::application::orchestration::fast_mode(),
            fast_mode: crate::application::orchestration::fast_mode(),
            best_per_generation: ga_result.best_per_generation.clone(),
        },
    )?;

    Ok(Milestone12GaRun {
        ranked,
        best_per_generation: ga_result.best_per_generation,
        audit,
        artifacts: vec![
            Milestone12StageArtifact {
                label: "ga_top5".to_string(),
                path: top5_path,
            },
            Milestone12StageArtifact {
                label: "ga_figure".to_string(),
                path: figure_path,
            },
        ],
    })
}
