use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::application::reporting::figures::{missing_required_figures, FigureManifestEntry};
use crate::application::reporting::narrative::{
    build_component_audit, render_goal_narrative, ComponentAuditEntry,
};
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationEvidence {
    pub robustness_complete: bool,
    pub validation_2d_complete: bool,
    pub validation_3d_complete: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoalEvidence {
    pub goal: OptimizationGoal,
    pub candidate_id: String,
    pub blueprint_name: String,
    pub evaluation: BlueprintObjectiveEvaluation,
    pub figures: Vec<FigureManifestEntry>,
    pub component_audit: Vec<ComponentAuditEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvidenceRunManifest {
    pub run_label: String,
    pub validation: ValidationEvidence,
    pub goals: Vec<GoalEvidence>,
}

pub fn build_goal_evidence(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintObjectiveEvaluation,
    figures: Vec<FigureManifestEntry>,
) -> GoalEvidence {
    GoalEvidence {
        goal: evaluation.goal,
        candidate_id: candidate.id.clone(),
        blueprint_name: candidate.blueprint.name.clone(),
        component_audit: build_component_audit(candidate),
        evaluation,
        figures,
    }
}

pub fn validate_canonical_manifest(manifest: &EvidenceRunManifest) -> Result<(), OptimError> {
    if !manifest.validation.robustness_complete
        || !manifest.validation.validation_2d_complete
        || !manifest.validation.validation_3d_complete
    {
        return Err(OptimError::InvalidParameter(
            "canonical report generation requires robustness, 2D validation, and 3D validation evidence".to_string(),
        ));
    }

    for goal in [
        OptimizationGoal::AsymmetricSplitResidenceSeparation,
        OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        OptimizationGoal::InPlaceDeanSerpentineRefinement,
    ] {
        let evidence = manifest
            .goals
            .iter()
            .find(|entry| entry.goal == goal)
            .ok_or_else(|| {
                OptimError::InvalidParameter(format!(
                    "canonical report is missing evidence for {:?}",
                    goal
                ))
            })?;
        if evidence.component_audit.is_empty() {
            return Err(OptimError::InvalidParameter(format!(
                "canonical report is missing component audit entries for {:?}",
                goal
            )));
        }
        let missing = missing_required_figures(goal, &evidence.figures);
        if !missing.is_empty() {
            return Err(OptimError::InvalidParameter(format!(
                "canonical report is missing required figures for {:?}: {}",
                goal,
                missing.join(", ")
            )));
        }
    }

    Ok(())
}

pub fn render_canonical_report(manifest: &EvidenceRunManifest) -> Result<String, OptimError> {
    validate_canonical_manifest(manifest)?;

    let mut report = String::new();
    report.push_str(&format!("# {}\n\n", manifest.run_label));
    report.push_str("## Evidence Summary\n\n");
    report.push_str(
        "All canonical evidence gates passed: robustness, 2D validation, and 3D validation.\n\n",
    );

    for goal_evidence in &manifest.goals {
        report.push_str(&format!("## {:?}\n\n", goal_evidence.goal));
        report.push_str(&render_goal_narrative(&goal_evidence.evaluation));
        report.push_str("\n\n");
        report.push_str("### Figures\n\n");
        for figure in &goal_evidence.figures {
            report.push_str(&format!(
                "- `{}`: {} ({})\n",
                figure.figure_id,
                figure.caption,
                figure.path.display()
            ));
        }
        report.push_str("\n### Component Audit\n\n");
        for component in &goal_evidence.component_audit {
            report.push_str(&format!(
                "- `{}` [{}] {} Safety: {}\n",
                component.component_id,
                component.component_type,
                component.governing_physics,
                component.safety_relevance
            ));
        }
        report.push('\n');
    }

    Ok(report)
}

#[must_use]
pub fn render_exploratory_report(manifest: &EvidenceRunManifest) -> String {
    let mut report = String::new();
    report.push_str(&format!("# Exploratory {}\n\n", manifest.run_label));
    for goal_evidence in &manifest.goals {
        report.push_str(&format!(
            "- {:?}: score {:.4}, residence {:.4} s, separation {:.4}\n",
            goal_evidence.goal,
            goal_evidence.evaluation.score_or_zero(),
            goal_evidence
                .evaluation
                .residence
                .treatment_residence_time_s,
            goal_evidence.evaluation.separation.separation_efficiency
        ));
    }
    report
}

pub fn write_canonical_report(
    manifest: &EvidenceRunManifest,
    output_path: &Path,
) -> Result<(), OptimError> {
    let report = render_canonical_report(manifest)?;
    fs::write(output_path, report).map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to write canonical report to {}: {error}",
            output_path.display()
        ))
    })
}

pub fn write_exploratory_report(
    manifest: &EvidenceRunManifest,
    output_path: &Path,
) -> Result<(), OptimError> {
    fs::write(output_path, render_exploratory_report(manifest)).map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to write exploratory report to {}: {error}",
            output_path.display()
        ))
    })
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::application::objectives::BlueprintObjectiveEvaluation;
    use crate::application::reporting::figures::{required_figure_ids, FigureManifestEntry};
    use crate::domain::fixtures::{
        canonical_option1_candidate, canonical_option2_candidate, operating_point,
    };
    use crate::domain::{BlueprintCandidate, OptimizationGoal};
    use crate::metrics::{
        BlueprintSafetyMetrics, BlueprintSeparationMetrics, BlueprintVenturiMetrics,
        ResidenceMetrics, VenturiPlacementMetrics,
    };

    use super::{
        build_goal_evidence, render_canonical_report, validate_canonical_manifest,
        EvidenceRunManifest, ValidationEvidence,
    };

    fn sample_candidate(goal: OptimizationGoal) -> BlueprintCandidate {
        let operating_point = operating_point(2.0e-6, 30_000.0, 0.18);
        match goal {
            OptimizationGoal::AsymmetricSplitResidenceSeparation => {
                canonical_option1_candidate(format!("{goal:?}"), operating_point)
            }
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity => {
                canonical_option2_candidate(format!("{goal:?}"), operating_point)
            }
            OptimizationGoal::InPlaceDeanSerpentineRefinement => crate::generate_ga_mutations(
                &canonical_option2_candidate("ga-seed", operating_point),
            )
            .expect("ga mutations")
            .into_iter()
            .next()
            .expect("at least one ga candidate"),
        }
    }

    fn sample_evaluation(goal: OptimizationGoal) -> BlueprintObjectiveEvaluation {
        BlueprintObjectiveEvaluation {
            goal,
            candidate_id: format!("{goal:?}"),
            blueprint_name: format!("{goal:?}"),
            score: Some(0.9),
            status: crate::BlueprintEvaluationStatus::Eligible,
            screening_reasons: Vec::new(),
            baseline_scores: Vec::new(),
            exceeds_all_baselines: None,
            residence: ResidenceMetrics {
                treatment_volume_m3: 1.0e-9,
                treatment_residence_time_s: 0.02,
                treatment_flow_fraction: 0.45,
                mean_treatment_velocity_m_s: 0.3,
            },
            separation: BlueprintSeparationMetrics {
                stage_summaries: Vec::new(),
                cancer_center_fraction: 0.9,
                wbc_center_fraction: 0.2,
                rbc_peripheral_fraction: 0.85,
                separation_efficiency: 0.8,
                center_hematocrit_ratio: 0.6,
            },
            venturi: BlueprintVenturiMetrics {
                placements: Vec::new(),
                cavitation_selectivity_score: 0.7,
                venturi_flow_fraction: 0.45,
                rbc_exposure_fraction: 0.15,
                wbc_exposure_fraction: 0.2,
            },
            safety: BlueprintSafetyMetrics {
                max_main_channel_shear_pa: 80.0,
                max_venturi_shear_pa: 2000.0,
                pressure_drop_pa: 18_000.0,
                pressure_feasible: true,
                main_channel_margin: 0.46,
                cavitation_safety_margin: 1.0,
                mean_device_residence_time_s: 0.03,
            },
        }
    }

    fn figure_entry(goal: OptimizationGoal, figure_id: &str) -> FigureManifestEntry {
        let mut path = std::env::temp_dir();
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        path.push(format!("{goal:?}-{figure_id}-{unique}.svg"));
        std::fs::write(&path, "<svg/>").expect("write figure");
        FigureManifestEntry {
            figure_id: figure_id.to_string(),
            caption: figure_id.to_string(),
            path,
        }
    }

    #[test]
    fn official_report_generation_rejects_missing_evidence() {
        let candidate = sample_candidate(OptimizationGoal::AsymmetricSplitResidenceSeparation);
        let evaluation = sample_evaluation(OptimizationGoal::AsymmetricSplitResidenceSeparation);
        let manifest = EvidenceRunManifest {
            run_label: "missing-evidence".to_string(),
            validation: ValidationEvidence {
                robustness_complete: false,
                validation_2d_complete: true,
                validation_3d_complete: true,
            },
            goals: vec![build_goal_evidence(&candidate, evaluation, Vec::new())],
        };
        assert!(validate_canonical_manifest(&manifest).is_err());
    }

    #[test]
    fn canonical_report_generation_contains_complete_appendices() {
        let goals = [
            OptimizationGoal::AsymmetricSplitResidenceSeparation,
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
            OptimizationGoal::InPlaceDeanSerpentineRefinement,
        ];
        let goal_evidence = goals
            .into_iter()
            .map(|goal| {
                let candidate = sample_candidate(goal);
                let mut evaluation = sample_evaluation(goal);
                if goal != OptimizationGoal::AsymmetricSplitResidenceSeparation {
                    evaluation.venturi.placements = vec![VenturiPlacementMetrics {
                        placement_id: "p0".to_string(),
                        target_channel_id: "stage0_ctc".to_string(),
                        cavitation_number: 0.8,
                        effective_throat_velocity_m_s: 12.0,
                        throat_static_pressure_pa: 20_000.0,
                        diffuser_recovery_pa: 0.0,
                        dean_number: 45.0,
                        curvature_radius_m: 0.002,
                        arc_length_m: 0.01,
                    }];
                }
                let figures = required_figure_ids(goal)
                    .iter()
                    .map(|figure_id| figure_entry(goal, figure_id))
                    .collect::<Vec<_>>();
                build_goal_evidence(&candidate, evaluation, figures)
            })
            .collect::<Vec<_>>();

        let manifest = EvidenceRunManifest {
            run_label: "canonical".to_string(),
            validation: ValidationEvidence {
                robustness_complete: true,
                validation_2d_complete: true,
                validation_3d_complete: true,
            },
            goals: goal_evidence,
        };

        let rendered = render_canonical_report(&manifest).expect("canonical report");
        assert!(rendered.contains("Component Audit"));
        assert!(rendered.contains("VenturiPlacement"));
        assert!(!rendered.to_ascii_lowercase().contains("placeholder"));
    }
}
