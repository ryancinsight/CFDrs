use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::{healthy_cell_protection_index, BlueprintEvaluation};

fn acoustic_residence_support_score(evaluation: &BlueprintEvaluation) -> f64 {
    let residence_norm = (evaluation.residence.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
    let flow_frac = evaluation.residence.treatment_flow_fraction.clamp(0.0, 1.0);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let cancer = evaluation.separation.cancer_center_fraction.clamp(0.0, 1.0);
    let healthy_cell_protection = healthy_cell_protection_index(
        evaluation.separation.wbc_center_fraction,
        evaluation.separation.rbc_peripheral_fraction,
    );
    let safety = evaluation.safety.main_channel_margin.clamp(0.0, 1.0);

    let base = 0.22 * cancer
        + 0.18 * sep
        + 0.16 * residence_norm
        + 0.12 * flow_frac
        + 0.22 * healthy_cell_protection
        + 0.10 * safety;
    let synergy = 0.12 * (sep * cancer * residence_norm * healthy_cell_protection).powf(0.25);
    (base + synergy).clamp(0.0, 1.0)
}

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
/// reward coupled Dean-vortex + cavitation performance while keeping the score
/// non-negative for successful evaluations.
pub fn evaluate_blueprint_genetic_refinement(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> Result<BlueprintObjectiveEvaluation, OptimError> {
    let topology = candidate.topology_spec().map_err(|e| {
        OptimError::InvalidParameter(format!(
            "GA refinement scoring requires canonical topology metadata for '{}': {e}",
            candidate.id
        ))
    })?;

    if topology.venturi_placements.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement scoring requires venturi treatment geometry, but candidate '{}' has no venturi placements",
            candidate.id
        )));
    }

    if !topology.has_serpentine() {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement scoring requires serpentine treatment geometry, but candidate '{}' has no serpentine channels",
            candidate.id
        )));
    }

    if !topology.venturi_placements.iter().any(|placement| {
        placement.placement_mode == cfd_schematics::VenturiPlacementMode::CurvaturePeakDeanNumber
    }) {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement scoring requires curvature-based Dean venturi placement, but candidate '{}' only carries non-Dean throat placements",
            candidate.id
        )));
    }

    let cav = evaluation
        .venturi
        .cavitation_selectivity_score
        .clamp(0.0, 1.0);
    let acoustic_support = acoustic_residence_support_score(&evaluation);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let safety = evaluation.safety.main_channel_margin.clamp(0.0, 1.0);
    let residence_norm = (evaluation.residence.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
    let flow_frac = evaluation.residence.treatment_flow_fraction.clamp(0.0, 1.0);
    let cancer = evaluation.separation.cancer_center_fraction.clamp(0.0, 1.0);
    let healthy_cell_protection = healthy_cell_protection_index(
        evaluation.venturi.wbc_exposure_fraction,
        1.0 - evaluation.venturi.rbc_exposure_fraction,
    );

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
    // Normalize Dean before multiplying to avoid large intermediates.
    let dean_norm = ((max_dean / 100.0).min(1.0) * bayat_factor.sqrt()).clamp(0.0, 1.0);

    let lineage_norm = candidate.blueprint.lineage().map_or(0.0, |lineage| {
        (lineage.mutations.len() as f64 / 5.0).clamp(0.0, 1.0)
    });

    let base = 0.34 * acoustic_support
        + 0.16 * cav
        + 0.10 * cancer
        + 0.08 * sep
        + 0.14 * healthy_cell_protection
        + 0.06 * safety
        + 0.07 * dean_norm
        + 0.05 * lineage_norm;

    // Curvature-driven secondary flow is rewarded only when it coexists with
    // strong treatment-lane enrichment and venturi cavitation support.
    let synergy_base = acoustic_support
        * cav
        * cancer
        * flow_frac
        * residence_norm
        * healthy_cell_protection
        * dean_norm;
    let synergy = 0.18 * synergy_base.powf(0.2);
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
        (base + synergy).clamp(0.0, 1.0),
    )
    .with_screening_reasons(screening_reasons))
}
