use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::BlueprintEvaluation;

/// Score for the asymmetric-split venturi-cavitation selectivity objective.
///
/// Evaluates how well venturi throat placements in treatment-zone channels
/// achieve selective hydrodynamic cavitation of circulating tumor cells while
/// limiting RBC and WBC exposure.  The key physics: asymmetric splits drive
/// differential flow rates and pressures; channel size determines whether
/// venturi throats can achieve σ < 1 cavitation.  Too many splits may reduce
/// channel widths below the effective cavitation threshold.
///
/// Uses an **additive weighted sum** (90%) plus a geometric-mean synergy (10%)
/// to prevent any single zero factor from eliminating the score.
///
/// Floor: 0.001 — no feasible candidate scores exactly zero.
pub fn evaluate_selective_venturi_cavitation(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    if candidate
        .topology_spec()
        .map_or(true, |topology| topology.venturi_placements.is_empty())
    {
        return Err(OptimError::InvalidParameter(format!(
            "Option 2 requires venturi treatment geometry, but candidate '{}' has no venturi placements",
            candidate.id
        )));
    }

    let cav = evaluation.venturi.cavitation_selectivity_score.clamp(0.0, 1.0);
    let rbc_shield = (1.0 - evaluation.venturi.rbc_exposure_fraction).clamp(0.0, 1.0);
    let wbc_shield = (1.0 - evaluation.venturi.wbc_exposure_fraction).clamp(0.0, 1.0);
    let safety = evaluation.safety.cavitation_safety_margin.clamp(0.0, 1.0);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let routing_support = evaluation.residence.treatment_flow_fraction.clamp(0.0, 1.0);

    let base = 0.32 * cav
        + 0.18 * rbc_shield
        + 0.14 * wbc_shield
        + 0.12 * safety
        + 0.08 * sep
        + 0.06 * routing_support;

    let synergy = 0.10 * (cav * rbc_shield * routing_support.max(0.01)).cbrt();

    let score = (base + synergy).clamp(0.001, 1.0);
    Ok(BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        candidate,
        evaluation,
        score,
    ))
}
