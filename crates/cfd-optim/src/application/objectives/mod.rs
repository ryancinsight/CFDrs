mod ga;
mod option1;
mod option2;

use std::sync::Arc;

use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::{
    evaluate_blueprint_candidate, BlueprintEvaluation, BlueprintSafetyMetrics,
    BlueprintSeparationMetrics, BlueprintVenturiMetrics, ResidenceMetrics,
};
use serde::{Deserialize, Serialize};

pub use ga::evaluate_blueprint_genetic_refinement;
pub use option1::{evaluate_selective_acoustic_residence_separation, score_selective_acoustic_residence_separation};
pub use option2::{evaluate_selective_venturi_cavitation, score_selective_venturi_cavitation};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BlueprintEvaluationStatus {
    Eligible,
    ScreenedOut,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintObjectiveEvaluation {
    pub goal: OptimizationGoal,
    pub candidate_id: String,
    pub blueprint_name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score: Option<f64>,
    pub status: BlueprintEvaluationStatus,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub screening_reasons: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub baseline_scores: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exceeds_all_baselines: Option<bool>,
    pub residence: ResidenceMetrics,
    pub separation: BlueprintSeparationMetrics,
    pub venturi: BlueprintVenturiMetrics,
    pub safety: BlueprintSafetyMetrics,
}

impl BlueprintObjectiveEvaluation {
    #[must_use]
    pub fn from_evaluation(
        goal: OptimizationGoal,
        candidate: &BlueprintCandidate,
        evaluation: BlueprintEvaluation,
        score: f64,
    ) -> Self {
        Self {
            goal,
            candidate_id: candidate.id.clone(),
            blueprint_name: candidate.blueprint.name.clone(),
            score: Some(score),
            status: BlueprintEvaluationStatus::Eligible,
            screening_reasons: Vec::new(),
            baseline_scores: Vec::new(),
            exceeds_all_baselines: None,
            residence: evaluation.residence,
            separation: evaluation.separation,
            venturi: evaluation.venturi,
            safety: evaluation.safety,
        }
    }

    /// Construct from an `Arc<BlueprintEvaluation>`, avoiding deep-cloning
    /// the heap-allocated `Vec` fields inside `BlueprintSeparationMetrics`
    /// and `BlueprintVenturiMetrics` when the Arc has only one strong
    /// reference.  Falls back to cloning from the shared reference when the
    /// Arc is shared.
    #[must_use]
    pub fn from_shared_evaluation(
        goal: OptimizationGoal,
        candidate: &BlueprintCandidate,
        evaluation: Arc<BlueprintEvaluation>,
        score: f64,
    ) -> Self {
        let candidate_id = candidate.id.clone();
        let blueprint_name = candidate.blueprint.name.clone();
        match Arc::try_unwrap(evaluation) {
            Ok(owned) => Self {
                goal,
                candidate_id,
                blueprint_name,
                score: Some(score),
                status: BlueprintEvaluationStatus::Eligible,
                screening_reasons: Vec::new(),
                baseline_scores: Vec::new(),
                exceeds_all_baselines: None,
                residence: owned.residence,
                separation: owned.separation,
                venturi: owned.venturi,
                safety: owned.safety,
            },
            Err(shared) => Self {
                goal,
                candidate_id,
                blueprint_name,
                score: Some(score),
                status: BlueprintEvaluationStatus::Eligible,
                screening_reasons: Vec::new(),
                baseline_scores: Vec::new(),
                exceeds_all_baselines: None,
                residence: shared.residence.clone(),
                separation: shared.separation.clone(),
                venturi: shared.venturi.clone(),
                safety: shared.safety.clone(),
            },
        }
    }

    /// Construct from pre-extracted identity strings and a shared evaluation.
    ///
    /// Used by [`EvaluatedPool`] which stores only the candidate ID and
    /// blueprint name (not the full `BlueprintCandidate`) to avoid OOM at
    /// scale.
    #[must_use]
    pub fn from_identity(
        goal: OptimizationGoal,
        candidate_id: &str,
        blueprint_name: &str,
        evaluation: Arc<BlueprintEvaluation>,
        score: f64,
    ) -> Self {
        let candidate_id = candidate_id.to_owned();
        let blueprint_name = blueprint_name.to_owned();
        match Arc::try_unwrap(evaluation) {
            Ok(owned) => Self {
                goal,
                candidate_id,
                blueprint_name,
                score: Some(score),
                status: BlueprintEvaluationStatus::Eligible,
                screening_reasons: Vec::new(),
                baseline_scores: Vec::new(),
                exceeds_all_baselines: None,
                residence: owned.residence,
                separation: owned.separation,
                venturi: owned.venturi,
                safety: owned.safety,
            },
            Err(shared) => Self {
                goal,
                candidate_id,
                blueprint_name,
                score: Some(score),
                status: BlueprintEvaluationStatus::Eligible,
                screening_reasons: Vec::new(),
                baseline_scores: Vec::new(),
                exceeds_all_baselines: None,
                residence: shared.residence.clone(),
                separation: shared.separation.clone(),
                venturi: shared.venturi.clone(),
                safety: shared.safety.clone(),
            },
        }
    }

    #[must_use]
    pub fn is_eligible(&self) -> bool {
        self.status == BlueprintEvaluationStatus::Eligible
    }

    #[must_use]
    pub fn score_or_zero(&self) -> f64 {
        self.score.unwrap_or(0.0)
    }

    #[must_use]
    pub fn with_screening_reasons(mut self, reasons: Vec<String>) -> Self {
        if !reasons.is_empty() {
            self.status = BlueprintEvaluationStatus::ScreenedOut;
            self.score = None;
            self.screening_reasons = reasons;
        }
        self
    }

    #[must_use]
    pub fn with_baseline_scores(mut self, baseline_scores: &[f64]) -> Self {
        self.baseline_scores = baseline_scores.to_vec();
        if !baseline_scores.is_empty() && self.score.is_some() {
            self.exceeds_all_baselines = Some(
                baseline_scores
                    .iter()
                    .all(|baseline| self.score_or_zero() > *baseline),
            );
        }
        self
    }
}

pub fn evaluate_goal(
    candidate: &BlueprintCandidate,
    goal: OptimizationGoal,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    let evaluation = evaluate_blueprint_candidate(candidate)?;
    match goal {
        OptimizationGoal::AsymmetricSplitResidenceSeparation => Ok(
            evaluate_selective_acoustic_residence_separation(candidate, evaluation),
        ),
        OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity => {
            evaluate_selective_venturi_cavitation(candidate, evaluation)
        }
        OptimizationGoal::InPlaceDeanSerpentineRefinement => {
            evaluate_blueprint_genetic_refinement(candidate, evaluation)
        }
    }
}
