use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::BlueprintEvaluation;

pub fn evaluate_blueprint_genetic_refinement(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    if evaluation.venturi.placements.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "candidate '{}' is not an Option 2 seed and cannot be ranked as a GA refinement",
            candidate.id
        )));
    }
    let lineage_bonus = candidate
        .blueprint
        .lineage()
        .map_or(0.0, |lineage| lineage.mutations.len() as f64 * 0.05);
    let dean_bonus = evaluation
        .venturi
        .placements
        .iter()
        .map(|placement| placement.dean_number)
        .fold(0.0_f64, f64::max)
        / 100.0;
    let score = evaluation.venturi.cavitation_selectivity_score.max(0.0)
        * evaluation.separation.separation_efficiency.max(0.0)
        * evaluation.safety.main_channel_margin.max(0.0)
        + lineage_bonus
        + dean_bonus;
    Ok(BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::BlueprintGeneticRefinement,
        candidate,
        evaluation,
        score,
    ))
}
