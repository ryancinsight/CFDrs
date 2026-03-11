use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::metrics::BlueprintEvaluation;

/// Score for the asymmetric-split residence-separation objective.
///
/// Evaluates how well mirrored asymmetric split families (bifurcation,
/// trifurcation, quadfurcation, and pentafurcation)
/// achieve differential flow routing for cell-size/rigidity-based separation
/// of circulating tumor cells from WBCs and RBCs, combined with adequate
/// residence time in the treatment zone.
///
/// ## Physics — Zweifach–Fung effect at asymmetric bifurcations
///
/// When a parent channel (e.g. 6 mm width) splits into branches of unequal
/// width (e.g. 1 mm / 2 mm / 3 mm trifurcation), three cooperating effects
/// route larger cells toward the wider center branch:
///
/// 1. **Fåhræus–Lindqvist plasma skimming** — the cell-free layer (~2 µm)
///    hugs both walls; the wider branch captures proportionally more of the
///    cell-rich core flow, while narrow branches draw from the cell-depleted
///    marginal layer.
/// 2. **Inertial focusing** — at millifluidic Re (10–100), particles with
///    diameter a satisfy a/D_h > 0.07 migrate to equilibrium positions
///    (Segré–Silberberg annulus) that map to the center branch at asymmetric
///    junctions.
/// 3. **Steric exclusion** — CTCs (∅ 15–25 µm) physically cannot enter the
///    narrowest exit branch when its width approaches the cell diameter.
///
/// Each additional splitting stage compounds cancer-center enrichment
/// multiplicatively: a 3-stage Tri→Tri→Tri tree achieves cancer_center_fraction
/// ~0.75+ vs ~0.556 for Tri→Tri.
///
/// ## Scoring
///
/// Uses a **hybrid additive + multiplicative synergy** form:
///   score = base + synergy
/// where:
///   base    = weighted sum of cancer focusing, separation efficiency,
///             residence support, healthy-cell shielding, and safety
///   synergy = 0.12 × (sep × cancer × residence × shielding)^(1/4)
///
/// This preserves a strictly positive score for every successful evaluation
/// while still rewarding coupled physics rather than only independent marginals.
pub fn evaluate_selective_acoustic_residence_separation(
    candidate: &BlueprintCandidate,
    evaluation: BlueprintEvaluation,
) -> BlueprintObjectiveEvaluation {
    let residence_norm = (evaluation.residence.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
    let flow_frac = evaluation.residence.treatment_flow_fraction.clamp(0.0, 1.0);
    let sep = evaluation.separation.separation_efficiency.clamp(0.0, 1.0);
    let cancer = evaluation.separation.cancer_center_fraction.clamp(0.0, 1.0);
    let wbc_exclusion = (1.0 - evaluation.separation.wbc_center_fraction).clamp(0.0, 1.0);
    let rbc_exclusion = evaluation
        .separation
        .rbc_peripheral_fraction
        .clamp(0.0, 1.0);
    let safety = evaluation.safety.main_channel_margin.clamp(0.0, 1.0);

    let base = 0.22 * cancer
        + 0.18 * sep
        + 0.16 * residence_norm
        + 0.12 * flow_frac
        + 0.12 * wbc_exclusion
        + 0.10 * rbc_exclusion
        + 0.10 * safety;
    let healthy_cell_shielding = (wbc_exclusion * rbc_exclusion).sqrt();
    let synergy = 0.12
        * (sep * cancer * residence_norm.max(0.01) * healthy_cell_shielding.max(0.01)).powf(0.25);
    let screening_reasons = [
        (
            evaluation.residence.treatment_flow_fraction <= 0.0,
            "treatment flow fraction must remain positive",
        ),
        (
            evaluation.residence.treatment_residence_time_s <= 0.0,
            "treatment residence time must remain positive",
        ),
        (
            evaluation.safety.main_channel_margin <= 0.0,
            "main-channel safety margin must remain positive",
        ),
    ]
    .into_iter()
    .filter(|(flag, _)| *flag)
    .map(|(_, reason)| reason.to_string())
    .collect();

    BlueprintObjectiveEvaluation::from_evaluation(
        OptimizationGoal::AsymmetricSplitResidenceSeparation,
        candidate,
        evaluation,
        (base + synergy).clamp(0.001, 1.0),
    )
    .with_screening_reasons(screening_reasons)
}
