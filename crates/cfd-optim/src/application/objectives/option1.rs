use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::metrics::BlueprintEvaluation;

/// Score for the asymmetric-split residence-separation objective.
///
/// Evaluates how well asymmetric channel splits (bifurcation/trifurcation)
/// achieve differential flow routing for cell-size/rigidity-based separation
/// of circulating tumor cells from WBCs and RBCs, combined with adequate
/// residence time in the treatment zone.
///
/// Uses an additive weighted sum plus a smooth synergy term so that physically
/// valid designs remain strictly positive while still rewarding coupled
/// behaviour: CTC enrichment into the treatment lane, suppression of WBC and
/// RBC exposure inside the treatment zone, and adequate residence time inside
/// the treatment window.
///
/// Floor: 0.001 — no feasible candidate scores exactly zero.
pub fn evaluate_selective_acoustic_residence_separation(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> BlueprintObjectiveEvaluation {
    let residence_norm = (evaluation.residence.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
    let flow_frac = evaluation.residence.treatment_flow_fraction.clamp(0.0, 1.0);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let cancer = evaluation.separation.cancer_center_fraction.clamp(0.0, 1.0);
    let wbc_exclusion = (1.0 - evaluation.separation.wbc_center_fraction).clamp(0.0, 1.0);
    let rbc_exclusion = evaluation.separation.rbc_peripheral_fraction.clamp(0.0, 1.0);
    let safety = evaluation.safety.main_channel_margin.clamp(0.0, 1.0);

    let base = 0.22 * cancer
        + 0.18 * sep
        + 0.16 * residence_norm
        + 0.12 * flow_frac
        + 0.12 * wbc_exclusion
        + 0.10 * rbc_exclusion
        + 0.10 * safety;
    let healthy_cell_shielding = (wbc_exclusion * rbc_exclusion).sqrt();
    let synergy =
        0.12 * (sep * cancer * residence_norm.max(0.01) * healthy_cell_shielding.max(0.01))
            .powf(0.25);

    let score = (base + synergy).clamp(0.001, 1.0);
    BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::AsymmetricSplitResidenceSeparation,
        candidate,
        evaluation,
        score,
    )
}
