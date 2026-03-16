use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::BlueprintEvaluation;

/// Score for the in-place Dean serpentine refinement objective (GA stage).
///
/// Evaluates how well in-place modifications — converting straight treatment
/// channels to serpentines — enhance the design through Dean vortices
/// (secondary flows from centrifugal forces in curved channels).
///
/// ## Physics — Dean instability and secondary flow
///
/// In a curved channel of hydraulic diameter D_h and bend radius R_c, the
/// Dean number De = Re √(D_h / 2R_c) characterises the competition between
/// centrifugal forcing and viscous dissipation. Above De ≈ 40, two counter-
/// rotating vortex rolls appear in the cross-section. **Centrifugal forces
/// peak at changes in direction** — specifically at bend apices where the
/// local curvature 1/R_c is maximal — driving larger particles (CTCs,
/// ∅ 15–25 µm) toward the outer wall while smaller particles (RBCs ∅ 8 µm,
/// WBCs ∅ 10–12 µm) remain in the recirculation zone.
///
/// This secondary focusing is complementary to Zweifach–Fung enrichment at
/// upstream bifurcations and to hydrodynamic cavitation at venturi throats:
/// Dean vortices pre-concentrate CTCs at the channel wall where they will
/// subsequently experience strongest collapse forces from venturi-nucleated
/// bubbles.
///
/// ## GA mutation operators
///
/// 1. **Serpentine insertion** — converts a treatment channel into a curved
///    Dean-active segment while preserving canonical lineage.
/// 2. **Venturi insertion / retargeting** — applies serial throats to any
///    treatment-zone channel through canonical schematics mutations.
/// 3. **Mirrored split-merge insertion** — inserts a local mirrored selective
///    split followed by canonical remerge inside the treatment zone.
///
/// ## Scoring
///
/// Uses an **additive weighted base** plus a multiplicative synergy term to
/// reward coupled Dean-vortex + cavitation performance while keeping every
/// successful evaluation strictly above zero.
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

    let cav = evaluation
        .venturi
        .cavitation_selectivity_score
        .clamp(0.0, 1.0);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let safety = evaluation.safety.main_channel_margin.clamp(0.0, 1.0);
    let residence_norm = (evaluation.residence.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
    let cancer = evaluation.separation.cancer_center_fraction.clamp(0.0, 1.0);
    let rbc_shield = (1.0 - evaluation.venturi.rbc_exposure_fraction).clamp(0.0, 1.0);
    let wbc_shield = (1.0 - evaluation.venturi.wbc_exposure_fraction).clamp(0.0, 1.0);

    let max_dean = evaluation
        .venturi
        .placements
        .iter()
        .map(|placement| placement.dean_number)
        .fold(0.0_f64, f64::max);

    // Bayat-Rezai (2017) millifluidic Dean enhancement factor.
    // For rectangular channels at Re < 500, the friction factor enhancement
    // from secondary flow is: f/f_straight = 1 + 0.085 * De^0.48
    // This amplifies the Dean score by ~25% for typical millifluidic Dean
    // numbers (De = 20-50), reflecting stronger secondary vortex focusing
    // than the classical Ito (1959) circular-tube prediction.
    let bayat_factor = cfd_1d::bayat_rezai_enhancement(max_dean);
    let dean_norm = (max_dean / 100.0 * bayat_factor.sqrt()).clamp(0.0, 1.0);

    let lineage_norm = candidate.blueprint.lineage().map_or(0.0, |lineage| {
        (lineage.mutations.len() as f64 / 5.0).clamp(0.0, 1.0)
    });

    let base = 0.20 * cav
        + 0.15 * cancer
        + 0.12 * sep
        + 0.12 * residence_norm
        + 0.10 * rbc_shield
        + 0.08 * wbc_shield
        + 0.08 * safety
        + 0.08 * dean_norm
        + 0.07 * lineage_norm;

    // Curvature-driven secondary flow is rewarded only when it coexists with
    // strong treatment-lane enrichment and venturi cavitation support.
    let synergy = 0.12
        * (cav * cancer * residence_norm.max(0.01) * rbc_shield.max(0.01) * dean_norm.max(0.01))
            .powf(0.2);
    let screening_reasons = [(
        evaluation.safety.main_channel_margin <= 0.0,
        "main-channel safety margin must remain positive",
    )]
    .into_iter()
    .filter(|(flag, _)| *flag)
    .map(|(_, reason)| reason.to_string())
    .collect();
    Ok(BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::InPlaceDeanSerpentineRefinement,
        candidate,
        evaluation,
        (base + synergy).clamp(0.001, 1.0),
    )
    .with_screening_reasons(screening_reasons))
}
