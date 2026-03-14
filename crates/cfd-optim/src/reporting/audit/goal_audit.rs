use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::application::objectives::BlueprintEvaluationStatus;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::evaluate_goal;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GoalAuditStatus {
    Eligible,
    ScreenedOut,
    Errored,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoalAuditEntry {
    pub candidate_id: String,
    pub blueprint_name: String,
    pub goal: OptimizationGoal,
    pub status: GoalAuditStatus,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reasons: Vec<String>,
    pub score: Option<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub baseline_scores: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exceeds_all_baselines: Option<bool>,
    pub treatment_residence_time_s: Option<f64>,
    pub separation_efficiency: Option<f64>,
    pub cavitation_selectivity_score: Option<f64>,
    pub peak_dean_number: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct GoalAuditArtifacts {
    pub eligible_path: PathBuf,
    pub screened_out_path: PathBuf,
    pub errored_path: PathBuf,
    pub summary_path: PathBuf,
}

#[must_use]
pub fn audit_goal_candidates(
    candidates: &[BlueprintCandidate],
    goal: OptimizationGoal,
) -> Vec<GoalAuditEntry> {
    candidates
        .iter()
        .map(|candidate| match evaluate_goal(candidate, goal) {
            Ok(evaluation) => GoalAuditEntry {
                candidate_id: candidate.id.clone(),
                blueprint_name: candidate.blueprint.name.clone(),
                goal,
                status: match evaluation.status {
                    BlueprintEvaluationStatus::Eligible => GoalAuditStatus::Eligible,
                    BlueprintEvaluationStatus::ScreenedOut => GoalAuditStatus::ScreenedOut,
                },
                reasons: evaluation.screening_reasons,
                score: evaluation.score,
                baseline_scores: evaluation.baseline_scores,
                exceeds_all_baselines: evaluation.exceeds_all_baselines,
                treatment_residence_time_s: Some(evaluation.residence.treatment_residence_time_s),
                separation_efficiency: Some(evaluation.separation.separation_efficiency),
                cavitation_selectivity_score: Some(evaluation.venturi.cavitation_selectivity_score),
                peak_dean_number: Some(
                    evaluation
                        .venturi
                        .placements
                        .iter()
                        .map(|placement| placement.dean_number)
                        .fold(0.0_f64, f64::max),
                ),
            },
            Err(error) => GoalAuditEntry {
                candidate_id: candidate.id.clone(),
                blueprint_name: candidate.blueprint.name.clone(),
                goal,
                status: GoalAuditStatus::Errored,
                reasons: vec![error.to_string()],
                score: None,
                baseline_scores: Vec::new(),
                exceeds_all_baselines: None,
                treatment_residence_time_s: None,
                separation_efficiency: None,
                cavitation_selectivity_score: None,
                peak_dean_number: None,
            },
        })
        .collect()
}

pub fn write_goal_audit_report(
    output_dir: &Path,
    stem: &str,
    entries: &[GoalAuditEntry],
) -> Result<GoalAuditArtifacts, OptimError> {
    fs::create_dir_all(output_dir).map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to create audit output directory {}: {error}",
            output_dir.display()
        ))
    })?;

    let eligible: Vec<&GoalAuditEntry> = entries
        .iter()
        .filter(|entry| entry.status == GoalAuditStatus::Eligible)
        .collect();
    let screened_out: Vec<&GoalAuditEntry> = entries
        .iter()
        .filter(|entry| entry.status == GoalAuditStatus::ScreenedOut)
        .collect();
    let errored: Vec<&GoalAuditEntry> = entries
        .iter()
        .filter(|entry| entry.status == GoalAuditStatus::Errored)
        .collect();

    let eligible_path = output_dir.join(format!("{stem}_eligible.json"));
    let screened_out_path = output_dir.join(format!("{stem}_screened_out.json"));
    let errored_path = output_dir.join(format!("{stem}_errored.json"));
    let summary_path = output_dir.join(format!("{stem}_summary.md"));

    fs::write(
        &eligible_path,
        serde_json::to_string_pretty(&eligible)
            .map_err(|error| OptimError::InvalidParameter(error.to_string()))?,
    )
    .map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to write eligible audit ledger {}: {error}",
            eligible_path.display()
        ))
    })?;
    fs::write(
        &screened_out_path,
        serde_json::to_string_pretty(&screened_out)
            .map_err(|error| OptimError::InvalidParameter(error.to_string()))?,
    )
    .map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to write screened-out audit ledger {}: {error}",
            screened_out_path.display()
        ))
    })?;
    fs::write(
        &errored_path,
        serde_json::to_string_pretty(&errored)
            .map_err(|error| OptimError::InvalidParameter(error.to_string()))?,
    )
    .map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to write errored audit ledger {}: {error}",
            errored_path.display()
        ))
    })?;

    let summary = format!(
        "# {stem} audit\n\n- total: {}\n- eligible: {}\n- screened_out: {}\n- errored: {}\n",
        entries.len(),
        eligible.len(),
        screened_out.len(),
        errored.len()
    );
    fs::write(&summary_path, summary).map_err(|error| {
        OptimError::InvalidParameter(format!(
            "failed to write audit summary {}: {error}",
            summary_path.display()
        ))
    })?;

    Ok(GoalAuditArtifacts {
        eligible_path,
        screened_out_path,
        errored_path,
        summary_path,
    })
}

#[cfg(test)]
mod tests {
    use super::{audit_goal_candidates, write_goal_audit_report, GoalAuditStatus};
    use crate::domain::fixtures::{canonical_option1_candidate, operating_point};
    use crate::OptimizationGoal;

    #[test]
    fn audit_report_writes_machine_readable_ledgers() {
        let candidate =
            canonical_option1_candidate("audit-opt1", operating_point(2.4e-6, 32_000.0, 0.12));
        let entries = audit_goal_candidates(
            &[candidate],
            OptimizationGoal::AsymmetricSplitResidenceSeparation,
        );
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].status, GoalAuditStatus::Eligible);

        let mut out_dir = std::env::temp_dir();
        out_dir.push(format!(
            "cfd-optim-audit-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("clock")
                .as_nanos()
        ));
        let artifacts =
            write_goal_audit_report(&out_dir, "option1_audit", &entries).expect("audit writes");
        assert!(artifacts.eligible_path.exists());
        assert!(artifacts.screened_out_path.exists());
        assert!(artifacts.errored_path.exists());
        assert!(artifacts.summary_path.exists());
    }
}
