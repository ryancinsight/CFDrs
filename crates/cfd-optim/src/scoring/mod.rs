//! Multi-objective scoring functions for SDT design candidates.
//!
//! Two optimisation modes are provided:
//!
//! - **[`OptimMode::SdtCavitation`]** — maximise cavitation intensity at the
//!   venturi throat while maintaining FDA haemolysis compliance in the main
//!   channels.
//!
//! - **[`OptimMode::UniformExposure`]** — maximise spatial uniformity of flow
//!   across all 36 treatment wells and maximise residence time in the exposure
//!   zone.
//!
//! - **[`OptimMode::SelectiveAcousticTherapy`]** — maximise selective
//!   center-lane enrichment for acoustically treated treatment channels.
//!
//! - **[`OptimMode::Combined`]** — weighted combination of both objectives.
//!
//! Any candidate that violates either hard constraint — non-feasible pressure
//! drop **or** FDA main-channel shear exceedance — receives a score of `0.0`.

mod types;

pub use types::{OptimMode, ScoreMode, SdtWeights};

use crate::constraints::{HI_PASS_LIMIT, PAI_PASS_LIMIT, THERAPEUTIC_HI_PASS_LIMIT};
use crate::metrics::SdtMetrics;

/// Combined-mode gate target: minimum leukapheresis sub-score for full credit.
const COMBINED_LEUKA_MIN_SCORE: f64 = 0.10;
/// Combined-mode gate target: minimum WBC recovery for full credit.
const COMBINED_LEUKA_MIN_WBC_RECOVERY: f64 = 0.20;
/// Combined-mode floor so GA/search retains a weak gradient signal.
const COMBINED_LEUKA_GATE_FLOOR: f64 = 0.02;
/// Combined-mode pediatric cumulative hemolysis target over 15 min.
const COMBINED_PEDIATRIC_HI15_LIMIT: f64 = 0.01;
/// Floor for the pediatric cumulative-HI gate to preserve search gradient.
const COMBINED_PEDIATRIC_HI15_GATE_FLOOR: f64 = 0.05;
/// Combined-mode minimum cancer-targeted cavitation for full oncology credit.
const COMBINED_ONCOLOGY_MIN_CANCER_CAV: f64 = 0.20;
/// Combined-mode minimum oncology selectivity index for full oncology credit.
const COMBINED_ONCOLOGY_MIN_SELECTIVITY: f64 = 0.10;
/// Floor for oncology gate to preserve optimization gradient.
const COMBINED_ONCOLOGY_GATE_FLOOR: f64 = 0.05;
/// Combined-mode selective-remerge proximity target score for full credit.
const COMBINED_SELECTIVE_REMERGE_MIN_SCORE: f64 = 0.70;
/// Floor for the selective-remerge gate to preserve optimization gradient.
const COMBINED_SELECTIVE_REMERGE_GATE_FLOOR: f64 = 0.20;

// ── Score functions ──────────────────────────────────────────────────────────

/// Compute the score for a single candidate given a mode and weights.
///
/// Returns a value in **[0.0, 1.0]**.  Higher scores indicate better designs.
/// Returns `0.0` for any candidate that violates hard constraints
/// (hard-constraint mode) or a tiny gradient signal (smooth-penalty mode).
///
/// # Theorem
/// In `HardConstraint` mode, `score_candidate(...) = 0` whenever any hard
/// feasibility predicate fails (`pressure_feasible == false`,
/// `fda_main_compliant == false`, or `plate_fits == false`).
///
/// **Proof sketch**
/// `score_candidate` delegates to `score_candidate_impl` with
/// `ScoreMode::HardConstraint`.  The first branch performs these feasibility
/// checks and returns `0.0` before any objective-specific scoring is evaluated.
/// Therefore no infeasible candidate can receive non-zero score in hard mode.
#[must_use]
pub fn score_candidate(metrics: &SdtMetrics, mode: OptimMode, weights: &SdtWeights) -> f64 {
    score_candidate_impl(metrics, mode, weights, ScoreMode::HardConstraint, 0.0)
}

/// Variant of [`score_candidate`] used by the genetic algorithm.
///
/// Uses smooth feasibility penalties instead of a hard `0.0` cliff, allowing the
/// GA to navigate toward the feasible region.  Pass `inlet_gauge_pa` from the
/// candidate so the pressure-margin sigmoid can be computed.
#[must_use]
pub fn score_candidate_ga(
    metrics: &SdtMetrics,
    mode: OptimMode,
    weights: &SdtWeights,
    inlet_gauge_pa: f64,
) -> f64 {
    score_candidate_impl(
        metrics,
        mode,
        weights,
        ScoreMode::SmoothPenalty,
        inlet_gauge_pa,
    )
}

fn score_candidate_impl(
    metrics: &SdtMetrics,
    mode: OptimMode,
    weights: &SdtWeights,
    constraint_mode: ScoreMode,
    inlet_gauge_pa: f64,
) -> f64 {
    // ── Feasibility check ─────────────────────────────────────────────────
    match constraint_mode {
        ScoreMode::HardConstraint => {
            // Disqualify any infeasible candidate.
            // `plate_fits` rejects candidates whose footprint exceeds the plate boundary.
            if !metrics.pressure_feasible || !metrics.fda_main_compliant || !metrics.plate_fits {
                return 0.0;
            }
            // Acoustic-only modes: strict 0.1 % per-pass limit.
            if matches!(
                mode,
                OptimMode::SdtTherapy | OptimMode::SelectiveAcousticTherapy
            ) && metrics.hemolysis_index_per_pass > HI_PASS_LIMIT
            {
                return 0.0;
            }
            // Hydrodynamic cavitation modes: relaxed 0.8 % per-pass limit
            // (FDA Class II extracorporeal therapeutic device).
            if matches!(
                mode,
                OptimMode::HydrodynamicCavitationSDT | OptimMode::CombinedSdtLeukapheresis { .. }
            ) && metrics.hemolysis_index_per_pass > THERAPEUTIC_HI_PASS_LIMIT
            {
                return 0.0;
            }
        }
        ScoreMode::SmoothPenalty => {
            // Smooth sigmoid multiplier: provides a non-zero gradient for the GA
            // to navigate from the infeasible region toward feasibility.
            let pressure_margin = if inlet_gauge_pa > 0.0 {
                (inlet_gauge_pa - metrics.total_pressure_drop_pa) / inlet_gauge_pa
            } else if metrics.pressure_feasible {
                1.0
            } else {
                -0.5
            };
            let fda_margin = (150.0 - metrics.max_main_channel_shear_pa) / 150.0;
            let hi_margin = if matches!(
                mode,
                OptimMode::HydrodynamicCavitationSDT | OptimMode::CombinedSdtLeukapheresis { .. }
            ) {
                (THERAPEUTIC_HI_PASS_LIMIT - metrics.hemolysis_index_per_pass)
                    / THERAPEUTIC_HI_PASS_LIMIT
            } else if matches!(
                mode,
                OptimMode::SdtTherapy | OptimMode::SelectiveAcousticTherapy
            ) {
                (HI_PASS_LIMIT - metrics.hemolysis_index_per_pass) / HI_PASS_LIMIT
            } else {
                0.2
            };
            let plate_ok = if metrics.plate_fits { 1.0_f64 } else { 0.0 };

            let feasibility = sigmoid_penalty(pressure_margin)
                * sigmoid_penalty(fda_margin)
                * sigmoid_penalty(hi_margin)
                * plate_ok;

            if feasibility < 0.05 {
                // Deeply infeasible: return a tiny gradient signal so the GA can escape.
                return feasibility * 0.1;
            }
            // Moderately feasible: continue with normal scoring, then multiply.
            let raw = score_mode_raw(metrics, mode, weights);
            let coag_penalty = coagulation_penalty(metrics, weights);
            return (raw * feasibility - coag_penalty).max(0.0);
        }
    }

    let raw = score_mode_raw(metrics, mode, weights);
    let coag_penalty = coagulation_penalty(metrics, weights);
    (raw - coag_penalty).max(0.0)
}

/// Compute the raw mode-specific score (without constraint checks or coag penalty).
fn score_mode_raw(metrics: &SdtMetrics, mode: OptimMode, weights: &SdtWeights) -> f64 {
    match mode {
        OptimMode::SdtCavitation => score_cavitation(metrics, weights),
        OptimMode::UniformExposure => score_exposure(metrics, weights),
        OptimMode::SelectiveAcousticTherapy => score_selective_acoustic_therapy(metrics),
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

/// Combined coagulation penalty (high-shear platelet activation + low-flow stasis clotting).
fn coagulation_penalty(metrics: &SdtMetrics, weights: &SdtWeights) -> f64 {
    let pai_term = (metrics.platelet_activation_index / PAI_PASS_LIMIT).min(1.0);
    let stasis_term = metrics.clotting_risk_index.clamp(0.0, 1.0);
    // Platelet activation remains primary; stasis clotting adds a complementary penalty.
    weights.coag_weight * (0.70 * pai_term + 0.30 * stasis_term).min(1.0)
}

/// Score for the SDT cavitation objective.
fn score_cavitation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let cav = metrics.cavitation_potential;

    // Smooth sigmoid HI penalty: avoids the clamp-at-2 plateau where designs
    // with HI=0.002 and HI=0.010 receive identical scores.
    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    let hi_factor = (1.0_f64 / (1.0 + (hi_ratio / 0.5).powi(2))).clamp(0.0_f64, 1.0_f64);

    let coverage = metrics.well_coverage_fraction.clamp(0.0, 1.0);

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

    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    let hi_factor = (1.0_f64 / (1.0 + (hi_ratio / 0.5).powi(2))).clamp(0.0_f64, 1.0_f64);

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

    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    let hi_factor = (1.0_f64 / (1.0 + (hi_ratio / 0.5).powi(2))).clamp(0.0_f64, 1.0_f64);

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
    let therapy_fraction = metrics.therapy_channel_fraction.clamp(0.0, 1.0);
    let residence = (metrics.mean_residence_time_s / 1.0).clamp(0.0, 1.0);
    let optical_405 = metrics.blue_light_delivery_index_405nm.clamp(0.0, 1.0);

    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    let hi_score = (1.0_f64 / (1.0_f64 + hi_ratio.powi(2))).clamp(0.0, 1.0);
    let clot_score = (1.0 - metrics.clotting_risk_index.clamp(0.0, 1.0)).clamp(0.0, 1.0);

    (0.20 * cancer_center
        + 0.15 * wbc_center
        + 0.15 * rbc_peripheral
        + 0.15 * sep3
        + 0.12 * therapy_fraction
        + 0.08 * residence
        + 0.06 * hi_score
        + 0.05 * clot_score
        + 0.04 * optical_405)
        .clamp(0.0, 1.0)
}

fn score_selective_acoustic_therapy(metrics: &SdtMetrics) -> f64 {
    let base = score_selective_routing_base(metrics);
    let therapy_fraction = metrics.therapy_channel_fraction.clamp(0.0, 1.0);
    let residence = (metrics.mean_residence_time_s / 1.0).clamp(0.0, 1.0);
    let well_coverage = metrics.well_coverage_fraction.clamp(0.0, 1.0);
    let optical_405 = metrics.blue_light_delivery_index_405nm.clamp(0.0, 1.0);

    (0.55 * base
        + 0.20 * therapy_fraction
        + 0.15 * residence
        + 0.05 * well_coverage
        + 0.05 * optical_405)
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
        (metrics.throat_transit_time_s / res_target_s)
            .min(1.0)
            .max(0.0)
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
        return 0.0;
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
    raw * ecv_ok_factor
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

    let base = w.hydro_cancer_cav * cancer_term
        + w.hydro_sep3 * sep3
        + w.hydro_rbc_protection * rbc_prot
        + w.hydro_sonolum * sono
        + w.hydro_wbc_cav * wbc_cav;

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
    let stasis_guard = (1.0 - metrics.clotting_risk_index.clamp(0.0, 1.0)).clamp(0.0, 1.0);
    let stasis_gate = 0.25 + 0.75 * stasis_guard;

    ((base + 0.08 * selective_delivery + synergy)
        * bias_gate
        * remerge_gate
        * optical_gate
        * stasis_gate)
        .min(1.0)
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
    let oncology_gate = COMBINED_ONCOLOGY_GATE_FLOOR
        + (1.0 - COMBINED_ONCOLOGY_GATE_FLOOR) * cav_gate_raw.min(sel_gate_raw);

    let gated_venturi_sdt = venturi_sdt * oncology_gate;
    let s_s = (0.60 * selective_base + 0.40 * gated_venturi_sdt).clamp(0.0, 1.0);
    let blended = (leuka_weight * s_l + sdt_weight * s_s) / (leuka_weight + sdt_weight).max(1e-12);

    // Safety gates: apply to the entire blended score.
    let wbc_gate = (metrics.wbc_recovery / COMBINED_LEUKA_MIN_WBC_RECOVERY).clamp(0.0, 1.0);
    let leuka_gate = (s_l / COMBINED_LEUKA_MIN_SCORE).clamp(0.0, 1.0);
    let leuka_gate_total =
        COMBINED_LEUKA_GATE_FLOOR + (1.0 - COMBINED_LEUKA_GATE_FLOOR) * wbc_gate.min(leuka_gate);
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
    let hi15_gate = COMBINED_PEDIATRIC_HI15_GATE_FLOOR
        + (1.0 - COMBINED_PEDIATRIC_HI15_GATE_FLOOR) * hi15_gate_raw;
    let cif_remerge_gate = if metrics.cif_outlet_tail_length_mm > 0.0 {
        let raw = (metrics.cif_remerge_proximity_score / COMBINED_SELECTIVE_REMERGE_MIN_SCORE)
            .clamp(0.0, 1.0);
        COMBINED_SELECTIVE_REMERGE_GATE_FLOOR + (1.0 - COMBINED_SELECTIVE_REMERGE_GATE_FLOOR) * raw
    } else {
        1.0
    };

    blended * leuka_gate_total * hi15_gate * cif_remerge_gate
}

// ── Constraint helpers ────────────────────────────────────────────────────────

/// Smooth sigmoid penalty based on a signed feasibility margin.
///
/// Returns `1.0` when the margin is `≥ +0.1` (well inside the feasible region),
/// `0.5` at the boundary (`margin = 0`), and `0.0` when margin `≤ −0.1`
/// (deeply infeasible).
///
/// # Theorem
/// `sigmoid_penalty(m)` is bounded and monotone:
/// 1. `0 ≤ sigmoid_penalty(m) ≤ 1` for all real `m`
/// 2. if `m1 ≤ m2`, then `sigmoid_penalty(m1) ≤ sigmoid_penalty(m2)`
///
/// **Proof sketch**
/// The function is affine (`0.5 + 5m`) followed by `clamp(0, 1)`.
/// Clamping maps all inputs into `[0,1]` and preserves monotonicity because
/// both the affine map and clamp are monotone non-decreasing.
#[inline]
#[must_use]
pub fn sigmoid_penalty(margin: f64) -> f64 {
    (0.5 + margin * 5.0).clamp(0.0, 1.0)
}

// ── Utility ──────────────────────────────────────────────────────────────────

/// Return a human-readable summary line for a scored candidate.
#[must_use]
pub fn score_description(mode: OptimMode) -> &'static str {
    match mode {
        OptimMode::SdtCavitation => "SDT Cavitation",
        OptimMode::UniformExposure => "Uniform Exposure",
        OptimMode::SelectiveAcousticTherapy => {
            "Selective Acoustic Therapy (Center Enrichment + Peripheral RBC Routing)"
        }
        OptimMode::Combined { .. } => "Combined (Cavitation + Exposure)",
        OptimMode::CellSeparation => "Cell Separation + SDT",
        OptimMode::ThreePopSeparation => "Three-Pop Separation (WBC+Cancer→Center, RBC→Wall) + SDT",
        OptimMode::SdtTherapy => "SDT Therapy (Selective Sep + HI + Cav + Dose)",
        OptimMode::PediatricLeukapheresis { .. } => {
            "Paediatric Leukapheresis (WBC Recovery + RBC Removal + Purity + ECV)"
        }
        OptimMode::HydrodynamicCavitationSDT => {
            "Hydrodynamic Cavitation SDT (Cancer-Targeted Cav + Sep3 + RBC Protection + Sono)"
        }
        OptimMode::CombinedSdtLeukapheresis { .. } => {
            "Combined SDT + Leukapheresis (WBC Recovery + RBC Removal + Cancer-Targeted Cavitation)"
        }
        OptimMode::RbcProtectedSdt => {
            "RBC-Protected SDT (Therapeutic Window + Cancer-Cav + Lysis Safety + FDA Compliance)"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn base_metrics() -> SdtMetrics {
        SdtMetrics {
            pressure_feasible: true,
            fda_main_compliant: true,
            plate_fits: true,
            hemolysis_index_per_pass: 2.0e-4,
            projected_hemolysis_15min_pediatric_3kg: 0.002,
            wbc_recovery: 0.70,
            rbc_pass_fraction: 0.30,
            wbc_purity: 0.70,
            total_ecv_ml: 0.2,
            flow_rate_ml_min: 120.0,
            three_pop_sep_efficiency: 0.30,
            cancer_targeted_cavitation: 0.20,
            oncology_selectivity_index: 0.10,
            cancer_rbc_cavitation_bias_index: 0.70,
            selective_cavitation_delivery_index: 0.30,
            rbc_venturi_protection: 0.40,
            sonoluminescence_proxy: 0.50,
            wbc_targeted_cavitation: 0.25,
            ..SdtMetrics::default()
        }
    }

    #[test]
    fn combined_mode_prefers_high_oncology_selectivity() {
        let mut low = base_metrics();
        low.cancer_targeted_cavitation = 0.04;
        low.oncology_selectivity_index = 0.01;

        let mut high = base_metrics();
        high.cancer_targeted_cavitation = 0.55;
        high.oncology_selectivity_index = 0.40;

        let mode = OptimMode::CombinedSdtLeukapheresis {
            leuka_weight: 0.5,
            sdt_weight: 0.5,
            patient_weight_kg: 3.0,
        };
        let w = SdtWeights::default();
        let s_low = score_candidate(&low, mode, &w);
        let s_high = score_candidate(&high, mode, &w);
        assert!(
            s_high > s_low,
            "combined mode should reward oncology-selective SDT candidates"
        );
    }

    #[test]
    fn hydrosdt_uses_oncology_selectivity_signal() {
        let mut low = base_metrics();
        low.oncology_selectivity_index = 0.05;
        low.cancer_targeted_cavitation = 0.35;

        let mut high = base_metrics();
        high.oncology_selectivity_index = 0.40;
        high.cancer_targeted_cavitation = 0.35;

        let w = SdtWeights::default();
        let mode = OptimMode::HydrodynamicCavitationSDT;
        let s_low = score_candidate(&low, mode, &w);
        let s_high = score_candidate(&high, mode, &w);
        assert!(
            s_high > s_low,
            "HydroSDT score should increase with oncology selectivity at fixed cavitation"
        );
    }

    #[test]
    fn hydrosdt_prefers_higher_cancer_rbc_cavitation_bias() {
        let mut low_bias = base_metrics();
        low_bias.cancer_rbc_cavitation_bias_index = 0.20;
        low_bias.selective_cavitation_delivery_index = 0.12;

        let mut high_bias = base_metrics();
        high_bias.cancer_rbc_cavitation_bias_index = 0.85;
        high_bias.selective_cavitation_delivery_index = 0.52;

        let w = SdtWeights::default();
        let mode = OptimMode::HydrodynamicCavitationSDT;
        let s_low = score_candidate(&low_bias, mode, &w);
        let s_high = score_candidate(&high_bias, mode, &w);
        assert!(
            s_high > s_low,
            "HydroSDT score should increase with stronger cancer-vs-RBC cavitation bias"
        );
    }

    #[test]
    fn combined_mode_prefers_cif_remerge_near_outlet() {
        let mode = OptimMode::CombinedSdtLeukapheresis {
            leuka_weight: 0.5,
            sdt_weight: 0.5,
            patient_weight_kg: 3.0,
        };
        let w = SdtWeights::default();

        let mut far_remerge = base_metrics();
        far_remerge.cif_outlet_tail_length_mm = 5.5;
        far_remerge.cif_remerge_proximity_score = 0.10;

        let mut near_remerge = base_metrics();
        near_remerge.cif_outlet_tail_length_mm = 0.8;
        near_remerge.cif_remerge_proximity_score = 0.92;

        let s_far = score_candidate(&far_remerge, mode, &w);
        let s_near = score_candidate(&near_remerge, mode, &w);
        assert!(
            s_near > s_far,
            "combined mode should reward selective remerge-near-outlet layouts"
        );
    }

    #[test]
    fn clotting_risk_penalizes_score() {
        let mode = OptimMode::CombinedSdtLeukapheresis {
            leuka_weight: 0.5,
            sdt_weight: 0.5,
            patient_weight_kg: 3.0,
        };
        let w = SdtWeights::default();

        let mut low_risk = base_metrics();
        low_risk.clotting_risk_index = 0.0;

        let mut high_risk = base_metrics();
        high_risk.clotting_risk_index = 1.0;

        let s_low = score_candidate(&low_risk, mode, &w);
        let s_high = score_candidate(&high_risk, mode, &w);
        assert!(
            s_low > s_high,
            "higher low-flow clotting risk should reduce score via coagulation penalty"
        );
    }

    #[test]
    fn selective_acoustic_prefers_selective_center_routing() {
        let mut broad = base_metrics();
        broad.cancer_center_fraction = 0.22;
        broad.wbc_center_fraction = 0.18;
        broad.rbc_peripheral_fraction_three_pop = 0.24;
        broad.three_pop_sep_efficiency = 0.14;
        broad.therapy_channel_fraction = 0.16;
        broad.mean_residence_time_s = 0.20;

        let mut selective = base_metrics();
        selective.cancer_center_fraction = 0.76;
        selective.wbc_center_fraction = 0.66;
        selective.rbc_peripheral_fraction_three_pop = 0.82;
        selective.three_pop_sep_efficiency = 0.61;
        selective.therapy_channel_fraction = 0.34;
        selective.mean_residence_time_s = 1.30;

        let mode = OptimMode::SelectiveAcousticTherapy;
        let w = SdtWeights::default();
        let broad_score = score_candidate(&broad, mode, &w);
        let selective_score = score_candidate(&selective, mode, &w);
        assert!(
            selective_score > broad_score,
            "selective acoustic mode should reward center-lane enrichment and peripheral RBC routing"
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn prop_hard_constraints_always_zero_score(
            cav in 0.0_f64..1.0,
            sep3 in 0.0_f64..1.0,
            wbc_recovery in 0.0_f64..1.0,
            optical_405 in 0.0_f64..1.0
        ) {
            let modes = [
                OptimMode::SdtCavitation,
                OptimMode::UniformExposure,
                OptimMode::SelectiveAcousticTherapy,
                OptimMode::Combined { cavitation_weight: 0.5, exposure_weight: 0.5 },
                OptimMode::CellSeparation,
                OptimMode::ThreePopSeparation,
                OptimMode::SdtTherapy,
                OptimMode::HydrodynamicCavitationSDT,
                OptimMode::CombinedSdtLeukapheresis { leuka_weight: 0.5, sdt_weight: 0.5, patient_weight_kg: 3.0 },
                OptimMode::RbcProtectedSdt,
            ];
            let weights = SdtWeights::default();

            for mode in modes {
                let mut m = base_metrics();
                m.cavitation_potential = cav;
                m.three_pop_sep_efficiency = sep3;
                m.wbc_recovery = wbc_recovery;
                m.blue_light_delivery_index_405nm = optical_405;

                m.pressure_feasible = false;
                prop_assert_eq!(score_candidate(&m, mode, &weights), 0.0);

                let mut m = base_metrics();
                m.cavitation_potential = cav;
                m.three_pop_sep_efficiency = sep3;
                m.wbc_recovery = wbc_recovery;
                m.blue_light_delivery_index_405nm = optical_405;
                m.fda_main_compliant = false;
                prop_assert_eq!(score_candidate(&m, mode, &weights), 0.0);

                let mut m = base_metrics();
                m.cavitation_potential = cav;
                m.three_pop_sep_efficiency = sep3;
                m.wbc_recovery = wbc_recovery;
                m.blue_light_delivery_index_405nm = optical_405;
                m.plate_fits = false;
                prop_assert_eq!(score_candidate(&m, mode, &weights), 0.0);
            }
        }

        #[test]
        fn prop_score_candidate_bounded_on_feasible_inputs(
            cav in 0.0_f64..1.0,
            sep3 in 0.0_f64..1.0,
            wbc_recovery in 0.0_f64..1.0,
            rbc_pass in 0.0_f64..1.0,
            optical_405 in 0.0_f64..1.0,
            clotting_risk in 0.0_f64..1.0
        ) {
            let modes = [
                OptimMode::SdtCavitation,
                OptimMode::UniformExposure,
                OptimMode::SelectiveAcousticTherapy,
                OptimMode::Combined { cavitation_weight: 0.5, exposure_weight: 0.5 },
                OptimMode::CellSeparation,
                OptimMode::ThreePopSeparation,
                OptimMode::SdtTherapy,
                OptimMode::HydrodynamicCavitationSDT,
                OptimMode::CombinedSdtLeukapheresis { leuka_weight: 0.5, sdt_weight: 0.5, patient_weight_kg: 3.0 },
                OptimMode::RbcProtectedSdt,
            ];
            let weights = SdtWeights::default();

            let mut m = base_metrics();
            m.pressure_feasible = true;
            m.fda_main_compliant = true;
            m.plate_fits = true;
            m.cavitation_potential = cav;
            m.three_pop_sep_efficiency = sep3;
            m.wbc_recovery = wbc_recovery;
            m.rbc_pass_fraction = rbc_pass;
            m.blue_light_delivery_index_405nm = optical_405;
            m.clotting_risk_index = clotting_risk;

            for mode in modes {
                let score = score_candidate(&m, mode, &weights);
                prop_assert!(score.is_finite());
                prop_assert!(score >= 0.0, "score = {score}");
                prop_assert!(score <= 1.0, "score = {score}");
            }
        }

        #[test]
        fn prop_sigmoid_penalty_bounded_and_monotone(
            m1 in -2.0_f64..2.0,
            m2 in -2.0_f64..2.0
        ) {
            let s1 = sigmoid_penalty(m1);
            let s2 = sigmoid_penalty(m2);
            prop_assert!((0.0..=1.0).contains(&s1));
            prop_assert!((0.0..=1.0).contains(&s2));
            if m1 <= m2 {
                prop_assert!(s1 <= s2 + 1e-12, "m1={m1}, m2={m2}, s1={s1}, s2={s2}");
            } else {
                prop_assert!(s2 <= s1 + 1e-12, "m1={m1}, m2={m2}, s1={s1}, s2={s2}");
            }
        }
    }
}
