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
pub use option1::evaluate_selective_acoustic_residence_separation;
pub use option2::evaluate_selective_venturi_cavitation;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintObjectiveEvaluation {
    pub goal: OptimizationGoal,
    pub candidate_id: String,
    pub blueprint_name: String,
    pub score: f64,
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
            score,
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
                score,
                residence: owned.residence,
                separation: owned.separation,
                venturi: owned.venturi,
                safety: owned.safety,
            },
            Err(shared) => Self {
                goal,
                candidate_id,
                blueprint_name,
                score,
                residence: shared.residence.clone(),
                separation: shared.separation.clone(),
                venturi: shared.venturi.clone(),
                safety: shared.safety.clone(),
            },
        }
    }
}

pub fn evaluate_goal(
    candidate: &BlueprintCandidate,
    goal: OptimizationGoal,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    let evaluation = evaluate_blueprint_candidate(candidate)?;
    match goal {
        OptimizationGoal::SelectiveAcousticResidenceSeparation => Ok(
            evaluate_selective_acoustic_residence_separation(candidate, evaluation),
        ),
        OptimizationGoal::SelectiveVenturiCavitation => {
            evaluate_selective_venturi_cavitation(candidate, evaluation)
        }
        OptimizationGoal::BlueprintGeneticRefinement => {
            evaluate_blueprint_genetic_refinement(candidate, evaluation)
        }
    }
}
