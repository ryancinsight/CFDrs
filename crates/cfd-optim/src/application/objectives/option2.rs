use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::BlueprintEvaluation;

/// Score for the asymmetric-split venturi-cavitation selectivity objective.
///
/// Evaluates how well venturi throat placements on mirrored treatment-zone
/// channels derived from the canonical split catalog
/// achieve selective hydrodynamic cavitation of circulating tumor cells while
/// limiting RBC and WBC exposure.
///
/// ## Physics — Rayleigh–Plesset bubble dynamics in venturi constrictions
///
/// The Bernoulli cavitation number σ = (p∞ − pᵥ)/(½ρv²_throat) governs
/// bubble nucleation: σ < 1 indicates that the local static pressure at
/// the vena contracta drops below the blood vapour pressure (pᵥ ≈ 6.3 kPa
/// at 37 °C). The resulting vapour/gas microbubble growth follows the
/// Rayleigh–Plesset ODE:
///
///   R R̈ + (3/2) Ṙ² = (1/ρ)[p_gas(R₀/R)^(3κ) − p∞ − pᵥ − 4μṘ/R − 2γ/R]
///
/// where R is the instantaneous bubble radius, κ ≈ 1.4 (polytropic exponent),
/// and γ ≈ 0.056 N/m (blood–air surface tension). Violent inertial collapse
/// at the divergent exit generates localised shear and sonoluminescence for SDT.
///
/// ## Channel size constraints
///
/// - **Lower bound** (throat ≥ 35 µm): avoids clogging by RBC rouleaux
///   (∅ ≈ 8 µm) and prevents fabrication limits from dominating pressure-
///   drop uncertainty.
/// - **Upper bound** (throat ≤ 120 µm): maintains v_throat sufficient for
///   σ < 1 within ≤ 250 kPa gauge (clinical extracorporeal pump ceiling).
///   Wider throats require impractically high flow rates for cavitation onset.
///
/// ## Scoring
///
/// Uses a **hybrid additive + multiplicative synergy** form:
///   score = base + synergy
/// where:
///   base    = weighted cavitation selectivity, healthy-cell shielding,
///             routing support, and safety support
///   synergy = 0.10 × (cav × rbc_shield × routing_support)^(1/3)
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

    let cav = evaluation
        .venturi
        .cavitation_selectivity_score
        .clamp(0.0, 1.0);
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
    let synergy = 0.10 * (cav * rbc_shield * routing_support).cbrt();
    let screening_reasons = [
        (
            evaluation.safety.cavitation_safety_margin <= 0.0,
            "cavitation safety margin must remain positive",
        ),
        (
            evaluation.residence.treatment_flow_fraction <= 0.0,
            "treatment flow fraction must remain positive",
        ),
        (
            evaluation.venturi.cavitation_selectivity_score <= 0.0,
            "venturi cavitation selectivity must remain positive",
        ),
    ]
    .into_iter()
    .filter(|(flag, _)| *flag)
    .map(|(_, reason)| reason.to_string())
    .collect();
    Ok(BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        candidate,
        evaluation,
        (base + synergy).clamp(0.001, 1.0),
    )
    .with_screening_reasons(screening_reasons))
}

/// Lightweight score-only variant that borrows the evaluation and returns
/// `Some(score)` if eligible, `None` if screened out. Zero allocations.
#[must_use]
pub fn score_selective_venturi_cavitation(
    evaluation: &BlueprintEvaluation,
    has_venturi_placements: bool,
) -> Option<f64> {
    if !has_venturi_placements {
        return None;
    }
    if evaluation.safety.cavitation_safety_margin <= 0.0
        || evaluation.residence.treatment_flow_fraction <= 0.0
        || evaluation.venturi.cavitation_selectivity_score <= 0.0
    {
        return None;
    }
    let cav = evaluation
        .venturi
        .cavitation_selectivity_score
        .clamp(0.0, 1.0);
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
    let synergy = 0.10 * (cav * rbc_shield * routing_support).cbrt();
    Some((base + synergy).clamp(0.001, 1.0))
}
