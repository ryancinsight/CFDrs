//! Mode-specific scoring functions for each [`OptimMode`] variant.
//!
//! Each function maps a [`SdtMetrics`] snapshot to a scalar fitness ∈ [0, 1].
//! The parent module's `score_candidate_impl` applies hard-constraint gates
//! *before* dispatching here, so these functions may assume that pressure
//! feasibility and FDA main-channel compliance have already been checked.

use crate::constraints::{HI_PASS_LIMIT, PAI_PASS_LIMIT};
use crate::metrics::SdtMetrics;

use super::types::{OptimMode, SdtWeights};

// ── Combined-mode gate constants ──────────────────────────────────────────────

/// Combined-mode gate target: minimum leukapheresis sub-score for full credit.
pub(super) const COMBINED_LEUKA_MIN_SCORE: f64 = 0.10;
/// Combined-mode gate target: minimum WBC recovery for full credit.
pub(super) const COMBINED_LEUKA_MIN_WBC_RECOVERY: f64 = 0.20;
/// Combined-mode floor so GA/search retains a weak gradient signal.
pub(super) const COMBINED_LEUKA_GATE_FLOOR: f64 = 0.02;
/// Combined-mode pediatric cumulative hemolysis target over 15 min.
pub(super) const COMBINED_PEDIATRIC_HI15_LIMIT: f64 = 0.01;
/// Floor for the pediatric cumulative-HI gate to preserve search gradient.
pub(super) const COMBINED_PEDIATRIC_HI15_GATE_FLOOR: f64 = 0.05;
/// Combined-mode minimum cancer-targeted cavitation for full oncology credit.
pub(super) const COMBINED_ONCOLOGY_MIN_CANCER_CAV: f64 = 0.20;
/// Combined-mode minimum oncology selectivity index for full oncology credit.
pub(super) const COMBINED_ONCOLOGY_MIN_SELECTIVITY: f64 = 0.10;
/// Floor for oncology gate to preserve optimization gradient.
pub(super) const COMBINED_ONCOLOGY_GATE_FLOOR: f64 = 0.05;
/// Combined-mode selective-remerge proximity target score for full credit.
pub(super) const COMBINED_SELECTIVE_REMERGE_MIN_SCORE: f64 = 0.70;
/// Floor for the selective-remerge gate to preserve optimization gradient.
pub(super) const COMBINED_SELECTIVE_REMERGE_GATE_FLOOR: f64 = 0.20;

// ── Shared hemolysis helpers ──────────────────────────────────────────────────

/// Hill-function hemolysis safety factor with half-max at HI = 0.5 × limit.
///
/// Maps hemolysis_index_per_pass to a [0, 1] safety score via a Hill function
/// of order 2: `1 / (1 + (HI / (0.5 × limit))²)`.  Returns 1.0 when HI = 0,
/// 0.5 when HI = 0.5 × limit, and decays quadratically above.
///
/// Used by cavitation, cell separation, and three-population separation modes
/// where the half-max threshold is set at 50% of the regulatory limit.
fn hemolysis_safety_factor(metrics: &SdtMetrics) -> f64 {
    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    (1.0_f64 / (1.0 + (hi_ratio / 0.5).powi(2))).clamp(0.0, 1.0)
}

/// Hemolysis compliance score with half-max at HI = limit.
///
/// Simpler Hill function `1 / (1 + (HI / limit)²)` for modes where the
/// half-max penalty sits at the regulatory limit itself rather than at 50%.
/// Returns 1.0 at HI = 0, 0.5 at HI = limit.
fn hemolysis_compliance_score(metrics: &SdtMetrics) -> f64 {
    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    (1.0_f64 / (1.0_f64 + hi_ratio.powi(2))).clamp(0.0, 1.0)
}

/// Clotting safety score: 1.0 = no risk, 0.0 = maximum stasis clotting.
fn clotting_safety_score(metrics: &SdtMetrics) -> f64 {
    (1.0 - metrics.clotting_risk_index.clamp(0.0, 1.0)).clamp(0.0, 1.0)
}

/// Apply a soft-floor gate: `floor + (1 − floor) × raw_value`.
///
/// Guarantees the gate never drops below `floor`, preserving gradient signal
/// for the optimizer while still penalizing designs that miss the target.
fn apply_gate_floor(raw_value: f64, floor: f64) -> f64 {
    floor + (1.0 - floor) * raw_value.clamp(0.0, 1.0)
}

// ── Dispatch ──────────────────────────────────────────────────────────────────

pub(super) fn score_mode_raw(metrics: &SdtMetrics, mode: OptimMode, weights: &SdtWeights) -> f64 {
    match mode {
        OptimMode::SdtCavitation => score_cavitation(metrics, weights),
        OptimMode::UniformExposure => score_exposure(metrics, weights),
        OptimMode::Combined {
            cavitation_weight,
            exposure_weight,
        } => {
            let s_cav = score_cavitation(metrics, weights);
            let s_exp = score_exposure(metrics, weights);
            (cavitation_weight * s_cav + exposure_weight * s_exp)
                / (cavitation_weight + exposure_weight).max(1e-12)
        }
        OptimMode::CellSeparation => score_cell_separation(metrics, weights),
        OptimMode::ThreePopSeparation => score_three_pop_separation(metrics, weights),
        OptimMode::SdtTherapy => score_sdt_therapy(metrics),
        OptimMode::PediatricLeukapheresis { patient_weight_kg } => {
            score_pediatric_leukapheresis(metrics, patient_weight_kg)
        }
        OptimMode::HydrodynamicCavitationSDT => score_hydrodynamic_cavitation_sdt(metrics, weights),
        OptimMode::CombinedSdtLeukapheresis {
            leuka_weight,
            sdt_weight,
            patient_weight_kg,
        } => score_combined_sdt_leukapheresis(
            metrics,
            weights,
            leuka_weight,
            sdt_weight,
            patient_weight_kg,
        ),
        OptimMode::RbcProtectedSdt => score_rbc_protected_sdt(metrics, weights),
    }
}

// ── Penalty helpers ───────────────────────────────────────────────────────────

/// Combined coagulation penalty (high-shear platelet activation + low-flow stasis clotting).
pub(super) fn coagulation_penalty(metrics: &SdtMetrics, weights: &SdtWeights) -> f64 {
    let pai_term = (metrics.platelet_activation_index / PAI_PASS_LIMIT).min(1.0);
    let stasis_term = metrics.clotting_risk_index.clamp(0.0, 1.0);
    // Platelet activation remains primary; stasis clotting adds a complementary penalty.
    weights.coag_weight * (0.70 * pai_term + 0.30 * stasis_term).min(1.0)
}

// ── Mode-specific scorers ─────────────────────────────────────────────────────

/// Score for the SDT cavitation objective.
///
/// Non-cavitating designs (σ ≥ 1) receive a small gradient signal
/// proportional to 1/σ, providing the GA with a smooth path toward
/// cavitation rather than a hard zero cliff.  The signal approaches
/// the full additive score as σ → 1⁺.
fn score_cavitation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let hi_factor = hemolysis_safety_factor(metrics);
    let coverage = metrics.well_coverage_fraction.clamp(0.0, 1.0);

    if metrics.cavitation_number >= 1.0 {
        // Non-cavitating: small gradient signal based on proximity to σ = 1
        // and the HI / coverage sub-scores.  A design at σ = 1.01 scores
        // ~0.10 × full, providing smooth gradient toward cavitation.
        let proximity = (1.0 / metrics.cavitation_number.max(1.0)).clamp(0.0, 1.0);
        let non_cav_base = w.cav_hemolysis * hi_factor + w.cav_coverage * coverage;
        return (0.10 * proximity * (non_cav_base + 0.01)).clamp(0.001, 0.10);
    }

    let cav = metrics.cavitation_potential;
    w.cav_potential * cav + w.cav_hemolysis * hi_factor + w.cav_coverage * coverage
}

/// Score for the uniform exposure objective.
fn score_exposure(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let uniformity = metrics.flow_uniformity.clamp(0.0, 1.0);
    let coverage = metrics.well_coverage_fraction.clamp(0.0, 1.0);
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);
    let optical_405 = metrics.blue_light_delivery_index_405nm.clamp(0.0, 1.0);

    let target_s = 30.0_f64;
    let res_norm = if metrics.mean_residence_time_s > 0.0 {
        (metrics.mean_residence_time_s.ln() - 1e-3_f64.ln()) / (target_s.ln() - 1e-3_f64.ln())
    } else {
        0.0
    }
    .clamp(0.0, 1.0);

    let core =
        w.exp_uniformity * uniformity + w.exp_coverage * coverage + w.exp_residence * res_norm;
    (0.70 * core + 0.15 * sep3 + 0.15 * optical_405).clamp(0.0, 1.0)
}

/// Score for the cell separation objective.
///
/// Uses **adaptive weighting** to avoid a floor effect when the achievable
/// separation is inherently low (< 0.25).
fn score_cell_separation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let separation = metrics.cell_separation_efficiency.clamp(0.0, 1.0);
    let cav = metrics.cavitation_potential;
    let hi_factor = hemolysis_safety_factor(metrics);

    let (sep_w, cav_w, hi_w) = if separation < 0.25 {
        (0.30, 0.50, 0.20)
    } else {
        (w.sep_efficiency, w.sep_cavitation, w.sep_hemolysis)
    };

    sep_w * separation + cav_w * cav + hi_w * hi_factor
}

/// Score for the three-population separation objective (WBC+cancer→center, RBC→wall).
fn score_three_pop_separation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);
    let wbc_center = metrics.wbc_center_fraction.clamp(0.0, 1.0);
    let cav = metrics.cavitation_potential;
    let hi_factor = hemolysis_safety_factor(metrics);

    w.sep3_efficiency * sep3
        + w.sep3_wbc_center * wbc_center
        + w.sep3_cavitation * cav
        + w.sep3_hemolysis * hi_factor
}

fn score_selective_routing_base(metrics: &SdtMetrics) -> f64 {
    let cancer_center = metrics.cancer_center_fraction.clamp(0.0, 1.0);
    let wbc_center = metrics.wbc_center_fraction.clamp(0.0, 1.0);
    let rbc_peripheral = metrics.rbc_peripheral_fraction_three_pop.clamp(0.0, 1.0);
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);
    // Use cancer_dose_fraction (cancer cells actually treated) rather than
    // therapy_channel_fraction (total blood through treatment zone).  For PST
    // designs the deliberate flow reduction is the selectivity mechanism, so
    // penalising it with therapy_fraction unfairly ranks PST3 below PST2.
    let cancer_dose = metrics.cancer_dose_fraction.clamp(0.0, 1.0);
    let residence = (metrics.mean_residence_time_s / 1.0).clamp(0.0, 1.0);
    let optical_405 = metrics.blue_light_delivery_index_405nm.clamp(0.0, 1.0);

    let hi_score = hemolysis_compliance_score(metrics);
    let clot_score = clotting_safety_score(metrics);

    (0.20 * cancer_center
        + 0.15 * wbc_center
        + 0.15 * rbc_peripheral
        + 0.15 * sep3
        + 0.12 * cancer_dose
        + 0.08 * residence
        + 0.06 * hi_score
        + 0.05 * clot_score
        + 0.04 * optical_405)
        .clamp(0.0, 1.0)
}

/// Score for the combined selective SDT therapy objective.
///
/// Rewards designs that achieve:
/// - High three-population separation (cancer → center, RBC → periphery)
/// - Low hemolysis index (safe for blood cells)
/// - High cavitation potential (effective SDT treatment)
/// - High cancer-targeted dose fraction
/// - RBC protection in peripheral bypass
/// - Adequate treatment residence time (throat exposure ≥ 1 ms target)
/// - Low cumulative hemolysis over a 15-minute pediatric therapy window
/// - Uniform wall-shear distribution (FDA spatial compliance)
fn score_sdt_therapy(metrics: &SdtMetrics) -> f64 {
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);

    // SDT therapy uses a milder Hill function (order 1) for hemolysis:
    // `1/(1 + HI/limit)` — half-max at HI = limit, linear decay, providing
    // softer gradient than the squared version in other modes.
    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    let hi_score = (1.0_f64 / (1.0_f64 + hi_ratio)).clamp(0.0_f64, 1.0_f64);

    let cav = metrics.cavitation_potential;
    let dose = metrics.cancer_dose_fraction.clamp(0.0, 1.0);
    let rbc_periph = metrics.rbc_peripheral_fraction_three_pop.clamp(0.0, 1.0);

    // Sonoluminescence proxy: sonosensitiser activation energy is proportional
    // to the adiabatic collapse temperature ratio.  This captures the
    // theraputically relevant dimension of cavitation intensity that
    // cavitation_potential alone does not — higher inlet pressure drives
    // more energetic collapse and stronger 5-ALA/Ce6 activation.
    let sono = metrics.sonoluminescence_proxy.clamp(0.0, 1.0);
    let oncology_selective = metrics.oncology_selectivity_index.clamp(0.0, 1.0);
    let optical_405 = metrics.blue_light_delivery_index_405nm.clamp(0.0, 1.0);

    // Throat residence time factor: reward designs where the cancer-enriched
    // stream spends at least ~1 ms in the cavitating throat (enough for R-P
    // bubble growth/collapse).  Uses throat_transit_time_s (throat-only transit,
    // ~1–5 ms) rather than mean_residence_time_s (full-chip, ~150–500 ms) so
    // that the factor discriminates between designs instead of always saturating.
    let res_target_s = 1e-3; // 1 ms target throat residence time
    let res_factor = if metrics.throat_transit_time_s > 0.0 {
        (metrics.throat_transit_time_s / res_target_s).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Cumulative hemolysis gate: penalise designs whose projected 15-minute
    // pediatric hemolysis exceeds the 1% clinical limit.  The gate goes from
    // 1.0 (when cumulative HI = 0) to 0.0 (when HI ≥ 2× the 1% limit).
    let cumul_hi_ratio =
        metrics.projected_hemolysis_15min_pediatric_3kg / COMBINED_PEDIATRIC_HI15_LIMIT.max(1e-18);
    let cumul_hi_gate = (1.0 / (1.0 + cumul_hi_ratio.powi(2))).clamp(0.0, 1.0);

    // Wall shear uniformity bonus: reward designs with low shear CV (per ASTM
    // F1841-20 spatial distribution requirement).  CV < 0.3 → full bonus;
    // CV > 1.0 → zero bonus.
    let shear_uniformity = (1.0 - (metrics.wall_shear_cv - 0.3).max(0.0) / 0.7).clamp(0.0, 1.0);

    let base = 0.16_f64 * sep3
        + 0.14_f64 * hi_score
        + 0.16_f64 * cav
        + 0.12_f64 * dose
        + 0.08_f64 * rbc_periph
        + 0.06_f64 * res_factor
        + 0.06_f64 * sono
        + 0.04_f64 * metrics.therapeutic_window_score.clamp(0.0, 1.0)
        + 0.04_f64 * oncology_selective
        + 0.04_f64 * shear_uniformity
        + 0.04_f64 * cumul_hi_gate
        + 0.06_f64 * optical_405;

    let synergy = if dose > 0.20 && cav > 0.50 && oncology_selective > 0.10 {
        0.05
    } else {
        0.0
    };
    let targeting_gate = (0.5 * dose + 0.5 * cav).clamp(0.0, 1.0);
    let non_target_penalty = if dose < 0.05 { 0.25 } else { 1.0 };

    ((base + synergy) * (0.65 + 0.35 * targeting_gate) * non_target_penalty).min(1.0_f64)
}

/// Score for the paediatric leukapheresis objective.
///
/// Weights: 40% WBC recovery + 30% RBC removal + 20% WBC purity + 10% throughput.
fn score_pediatric_leukapheresis(metrics: &SdtMetrics, patient_weight_kg: f64) -> f64 {
    if metrics.total_ecv_ml <= 0.0 {
        return 0.001;
    }

    let wbc_rec = metrics.wbc_recovery.clamp(0.0, 1.0);
    let has_wbc_capture = wbc_rec > 1.0e-9;
    let rbc_rem = if has_wbc_capture {
        (1.0_f64 - metrics.rbc_pass_fraction).clamp(0.0_f64, 1.0_f64)
    } else {
        0.0
    };
    let purity = if has_wbc_capture {
        metrics.wbc_purity.clamp(0.0, 1.0)
    } else {
        0.0
    };

    let throughput_score = (metrics.flow_rate_ml_min / 10.0).min(1.0);

    let patient_bv_ml = patient_weight_kg * 85.0;
    let max_ecv_ml = patient_bv_ml * 0.10;
    let ecv_ok_factor = if metrics.total_ecv_ml <= max_ecv_ml {
        1.0
    } else {
        (max_ecv_ml / metrics.total_ecv_ml).clamp(0.0, 1.0)
    };

    let raw = 0.40 * wbc_rec + 0.30 * rbc_rem + 0.20 * purity + 0.10 * throughput_score;
    (raw * ecv_ok_factor).max(0.001)
}

/// Score for the hydrodynamic cavitation SDT objective.
///
/// Rewards designs that route cancer cells into the cavitating venturi center
/// stream while keeping RBCs in peripheral bypass arms.
fn score_hydrodynamic_cavitation_sdt(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let cancer_cav = metrics.cancer_targeted_cavitation.clamp(0.0, 1.0);
    let oncology_selective = metrics.oncology_selectivity_index.clamp(0.0, 1.0);
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);
    let rbc_prot = metrics.rbc_venturi_protection.clamp(0.0, 1.0);
    let sono = metrics.sonoluminescence_proxy.clamp(0.0, 1.0);
    let wbc_cav = metrics.wbc_targeted_cavitation.clamp(0.0, 1.0);
    let cancer_term = (0.70 * cancer_cav + 0.30 * oncology_selective).clamp(0.0, 1.0);

    // WBC exclusion: reward designs where WBCs are EXCLUDED from the
    // venturi treatment zone.  Higher wbc_cav means more WBC exposure
    // to cavitation (bad for healthy cells), so we invert it.
    let wbc_exclusion = (1.0 - metrics.wbc_recovery.clamp(0.0, 1.0)).clamp(0.0, 1.0);

    let base = w.hydro_cancer_cav * cancer_term
        + w.hydro_sep3 * sep3
        + w.hydro_rbc_protection * rbc_prot
        + w.hydro_sonolum * sono
        + w.hydro_wbc_cav * wbc_exclusion;

    // Synergy bonus: strong cancer targeting AND effective cell separation.
    let synergy = if cancer_term > 0.30 && sep3 > 0.30 {
        0.05
    } else {
        0.0
    };
    let selective_delivery = metrics.selective_cavitation_delivery_index.clamp(0.0, 1.0);
    let cav_bias = metrics.cancer_rbc_cavitation_bias_index.clamp(0.0, 1.0);
    let optical_405 = metrics.blue_light_delivery_index_405nm.clamp(0.0, 1.0);
    let bias_gate = 0.15 + 0.85 * cav_bias;
    let remerge_gate = if metrics.cif_outlet_tail_length_mm > 0.0 {
        0.35 + 0.65 * metrics.cif_remerge_proximity_score.clamp(0.0, 1.0)
    } else {
        1.0
    };
    let optical_gate = 0.40 + 0.60 * optical_405;
    let stasis_gate = 0.25 + 0.75 * clotting_safety_score(metrics);

    ((base + 0.08 * selective_delivery + synergy)
        * bias_gate
        * remerge_gate
        * optical_gate
        * stasis_gate)
        .clamp(0.001, 1.0)
}

/// Score for the RBC-protected SDT objective.
///
/// Maximises the therapeutic window — the ratio of cancer treatment intensity
/// to RBC lysis risk — while rewarding FDA compliance and therapy-zone coverage.
fn score_rbc_protected_sdt(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let window = metrics.therapeutic_window_score.clamp(0.0, 1.0);
    let cancer_cav = metrics.cancer_targeted_cavitation.clamp(0.0, 1.0);

    // Lysis penalty: lysis_risk_index = 0 → score 1.0; index = 0.001 (limit) → 0.0.
    let lysis_penalty = (1.0 - metrics.lysis_risk_index * 1000.0).clamp(0.0, 1.0);

    // FDA bonus: full compliance (main + throat transit exception) → 1.0;
    // throat-exceeds-FDA without exception → 0.5 (partial credit: main is OK).
    let fda_bonus = if metrics.fda_overall_compliant {
        1.0_f64
    } else {
        0.5_f64
    };

    let coverage = metrics.therapy_channel_fraction.clamp(0.0, 1.0);

    let base = w.rbc_protected_window * window
        + w.rbc_protected_cav * cancer_cav
        + w.rbc_protected_lysis * lysis_penalty
        + w.rbc_protected_fda * fda_bonus
        + w.rbc_protected_coverage * coverage;

    // Synergy: strong cancer treatment AND low lysis risk — best-case scenario.
    let synergy = if cancer_cav > 0.25 && metrics.lysis_risk_index < 0.001 {
        0.05
    } else {
        0.0
    };

    (base + synergy).min(1.0)
}

/// Score for the combined SDT + paediatric leukapheresis objective.
///
/// Blends `score_pediatric_leukapheresis` and `score_hydrodynamic_cavitation_sdt`
/// with caller-supplied weights, enabling joint optimisation that satisfies both
/// the leukapheresis WBC-recovery goal and the cancer-targeted cavitation goal.
///
/// The oncology gate suppresses the venturi-SDT term when cavitation metrics
/// are weak, without penalising routing or leukapheresis value.  The
/// cancer_rbc_cavitation_bias gate is intentionally NOT applied here because
/// it is already incorporated inside `score_hydrodynamic_cavitation_sdt`;
/// double-gating would compound the penalty spuriously.
///
/// Safety gates (leukapheresis recovery, cumulative hemolysis, selective remerge)
/// apply to the entire blended score.
fn score_combined_sdt_leukapheresis(
    metrics: &SdtMetrics,
    weights: &SdtWeights,
    leuka_weight: f64,
    sdt_weight: f64,
    patient_weight_kg: f64,
) -> f64 {
    let s_l = score_pediatric_leukapheresis(metrics, patient_weight_kg);
    let selective_base = score_selective_routing_base(metrics);
    let venturi_sdt = score_hydrodynamic_cavitation_sdt(metrics, weights);

    // Oncology gate: suppress venturi-SDT when cavitation or selectivity is
    // too weak to deliver meaningful cancer treatment.
    let cav_gate_raw =
        (metrics.cancer_targeted_cavitation / COMBINED_ONCOLOGY_MIN_CANCER_CAV).clamp(0.0, 1.0);
    let sel_gate_raw =
        (metrics.oncology_selectivity_index / COMBINED_ONCOLOGY_MIN_SELECTIVITY).clamp(0.0, 1.0);
    let oncology_gate = apply_gate_floor(cav_gate_raw.min(sel_gate_raw), COMBINED_ONCOLOGY_GATE_FLOOR);

    let gated_venturi_sdt = venturi_sdt * oncology_gate;
    let s_s = (0.60 * selective_base + 0.40 * gated_venturi_sdt).clamp(0.0, 1.0);
    let blended = (leuka_weight * s_l + sdt_weight * s_s) / (leuka_weight + sdt_weight).max(1e-12);

    // Safety gates: apply to the entire blended score.
    let wbc_gate = (metrics.wbc_recovery / COMBINED_LEUKA_MIN_WBC_RECOVERY).clamp(0.0, 1.0);
    let leuka_gate = (s_l / COMBINED_LEUKA_MIN_SCORE).clamp(0.0, 1.0);
    let leuka_gate_total = apply_gate_floor(wbc_gate.min(leuka_gate), COMBINED_LEUKA_GATE_FLOOR);
    // Cumulative hemolysis gate: use patient-weight-appropriate metric.
    // For adult patients (≥ 40 kg), the pediatric 3 kg projection is
    // irrelevant — an adult at 70 kg has 5950 mL blood volume vs 255 mL for
    // a neonate, so pass count and projected HI are drastically different.
    let projected_hi15 = if patient_weight_kg >= 40.0 {
        metrics.projected_hemolysis_15min_adult
    } else {
        metrics.projected_hemolysis_15min_pediatric_3kg
    };
    let hi15_gate_raw = (1.0 - projected_hi15 / COMBINED_PEDIATRIC_HI15_LIMIT).clamp(0.0, 1.0);
    let hi15_gate = apply_gate_floor(hi15_gate_raw, COMBINED_PEDIATRIC_HI15_GATE_FLOOR);
    let cif_remerge_gate = if metrics.cif_outlet_tail_length_mm > 0.0 {
        let raw = (metrics.cif_remerge_proximity_score / COMBINED_SELECTIVE_REMERGE_MIN_SCORE)
            .clamp(0.0, 1.0);
        apply_gate_floor(raw, COMBINED_SELECTIVE_REMERGE_GATE_FLOOR)
    } else {
        1.0
    };

    // Floor the gated product: the individual gate floors (0.02, 0.05, 0.20)
    // multiply to ~0.0002, which destroys relative ranking when all candidates
    // land at the caller's INFEASIBILITY_FLOOR.  Flooring here preserves
    // meaningful ordering among gated candidates.
    (blended * leuka_gate_total * hi15_gate * cif_remerge_gate).max(0.001)
}
