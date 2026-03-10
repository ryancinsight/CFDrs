use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::metrics::BlueprintEvaluation;

pub fn evaluate_selective_acoustic_residence_separation(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> BlueprintObjectiveEvaluation {
    let residence_term = evaluation.residence.treatment_residence_time_s
        * evaluation.residence.treatment_flow_fraction.max(1.0e-9);
    let separation_term = evaluation.separation.separation_efficiency
        * evaluation.separation.cancer_center_fraction.max(1.0e-9);
    let safety_term = evaluation.safety.main_channel_margin.max(0.0);
    let score = residence_term * separation_term * safety_term;
    BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::SelectiveAcousticResidenceSeparation,
        candidate,
        evaluation,
        score,
    )
}
