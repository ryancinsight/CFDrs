use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::BlueprintEvaluation;

/// Score for the in-place Dean serpentine refinement objective (GA stage).
///
/// Evaluates how well in-place modifications — converting straight treatment
/// channels to serpentines — enhance the design through Dean vortices
/// (secondary flows from centrifugal forces in curved channels).  Peak Dean
/// number occurs at the outer wall of bends, driving enhanced cell focusing
/// and separation.
///
/// The GA mutates Milestone 12 seeds by:
/// 1. Adjusting branch widths (±8% / ±6%)
/// 2. Adding/extending serpentine segments on treatment paths
/// 3. Introducing or narrowing venturi throats (−8% when already present)
///
/// Uses an **additive weighted sum** (90%) plus a geometric-mean synergy (10%)
/// to prevent any single zero factor from eliminating the score.
///
/// Floor: 0.001 — no feasible candidate scores exactly zero.
pub fn evaluate_blueprint_genetic_refinement(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    if candidate
        .topology_spec()
        .map_or(true, |topology| topology.venturi_placements.is_empty())
    {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement scoring requires venturi treatment geometry, but candidate '{}' has no venturi placements",
            candidate.id
        )));
    }

    let cav = evaluation.venturi.cavitation_selectivity_score.clamp(0.0, 1.0);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let safety = evaluation.safety.main_channel_margin.clamp(0.0, 1.0);

    let dean_norm = evaluation
        .venturi
        .placements
        .iter()
        .map(|placement| placement.dean_number)
        .fold(0.0_f64, f64::max)
        / 100.0;
    let dean_norm = dean_norm.clamp(0.0, 1.0);

    let lineage_norm = candidate
        .blueprint
        .lineage()
        .map_or(0.0, |lineage| (lineage.mutations.len() as f64 / 5.0).clamp(0.0, 1.0));

    let base = 0.30 * cav
        + 0.20 * sep
        + 0.15 * safety
        + 0.15 * dean_norm
        + 0.10 * lineage_norm;

    // Dean vortices arise from curvature-driven secondary flow, with the
    // strongest cross-stream migration at bend apices and the outer wall.
    let synergy = 0.10 * (cav * sep * dean_norm.max(0.01)).cbrt();

    let score = (base + synergy).clamp(0.001, 1.0);
    Ok(BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::InPlaceDeanSerpentineRefinement,
        candidate,
        evaluation,
        score,
    ))
}
