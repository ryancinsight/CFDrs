use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::BlueprintEvaluation;

pub fn evaluate_selective_venturi_cavitation(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    if evaluation.venturi.placements.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "candidate '{}' has no venturi placements for selective cavitation evaluation",
            candidate.id
        )));
    }
    let cavitation_term = evaluation.venturi.cavitation_selectivity_score.max(0.0);
    let exposure_penalty = 1.0
        - 0.5
            * (evaluation.venturi.rbc_exposure_fraction + evaluation.venturi.wbc_exposure_fraction);
    let safety_term = evaluation.safety.cavitation_safety_margin.max(0.0);
    let score = cavitation_term * exposure_penalty.max(0.0) * safety_term;
    Ok(BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::SelectiveVenturiCavitation,
        candidate,
        evaluation,
        score,
    ))
}
