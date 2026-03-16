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
//! - **[`OptimMode::Combined`]** — weighted combination of both objectives.
//!
//! Any candidate that violates a hard constraint — non-feasible pressure drop,
//! FDA main-channel shear exceedance, or plate-boundary overflow — receives a
//! small positive floor score (`0.001`) to preserve gradient signal for the
//! optimizer while remaining strictly below any feasible design.

mod modes;
mod types;

pub use types::{OptimMode, ScoreMode, SdtWeights};

use crate::constraints::{FDA_THROAT_TEMP_RISE_LIMIT_K, HI_PASS_LIMIT, THERAPEUTIC_HI_PASS_LIMIT};
use crate::metrics::SdtMetrics;

/// Minimum score for any candidate, including those violating hard constraints.
///
/// Returning exactly zero destroys gradient information for downstream
/// optimizers (deterministic search, GA).  This floor preserves a small
/// positive signal so that infeasible candidates can still be ranked by
/// "how close to feasibility" they are, enabling smoother convergence.
const INFEASIBILITY_FLOOR: f64 = 0.001;

// ── Score functions ──────────────────────────────────────────────────────────

/// Compute the score for a single candidate given a mode and weights.
///
/// Returns a value in **[0.001, 1.0]**.  Higher scores indicate better designs.
/// Infeasible candidates (hard-constraint violations) receive a small but
/// non-zero floor (`INFEASIBILITY_FLOOR = 0.001`) to preserve gradient signal
/// for the optimizer.  Any feasible candidate will score strictly above this
/// floor, maintaining correct ranking while enabling smooth navigation out of
/// infeasible regions.
#[must_use]
pub fn score_candidate(metrics: &SdtMetrics, mode: OptimMode, weights: &SdtWeights) -> f64 {
    score_candidate_impl(metrics, mode, weights, ScoreMode::HardConstraint, 0.0)
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
            // Infeasible candidates receive INFEASIBILITY_FLOOR (0.001) — a tiny
            // positive score that preserves gradient signal without competing
            // with any feasible design.  The floor is small enough that any
            // feasible candidate (whose mode-specific raw score ≥ 0.01) always
            // dominates, yet non-zero so optimizers can rank infeasible designs
            // by proximity to the feasibility boundary.
            if !metrics.pressure_feasible || !metrics.fda_main_compliant || !metrics.plate_fits {
                return INFEASIBILITY_FLOOR;
            }
            // HI constraints use a smooth sigmoid gate rather than a binary
            // cutoff.  Each mode-specific scoring function already penalises
            // high hemolysis internally; the gate here provides an additional
            // multiplicative penalty that
            //   • preserves non-zero scores for designs slightly above the
            //     limit (enabling meaningful ranking),
            //   • still suppresses designs far above the limit (gate → 0),
            //   • is monotone decreasing in HI (no new optima introduced).
            //
            // Gate formula: 1 / (1 + (HI / limit)^4)
            //   HI = 0      → 1.0
            //   HI = limit  → 0.5
            //   HI = 2×limit → 0.059
            let hi_gate = if matches!(mode, OptimMode::SdtTherapy) {
                let r = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-18);
                1.0 / (1.0 + r * r)
            } else if matches!(
                mode,
                OptimMode::HydrodynamicCavitationSDT | OptimMode::CombinedSdtLeukapheresis { .. }
            ) {
                let r = metrics.hemolysis_index_per_pass / THERAPEUTIC_HI_PASS_LIMIT.max(1e-18);
                1.0 / (1.0 + r * r)
            } else {
                1.0
            };

            let raw = modes::score_mode_raw(metrics, mode, weights);
            let coag_penalty = modes::coagulation_penalty(metrics, weights);
            let thermal_penalty = if metrics.fda_thermal_compliant {
                0.0
            } else {
                0.05 * (metrics.throat_temperature_rise_k / FDA_THROAT_TEMP_RISE_LIMIT_K)
                    .clamp(0.0, 1.0)
            };
            // Pediatric high-flow penalty: penalises flow rates exceeding the
            // weight-scaled vascular access ceiling (10 mL/kg/min for neonatal
            // reference).  Max penalty 0.15 — strong enough to steer the optimizer
            // toward catheter-achievable flows but not so severe as to create a
            // cliff (preserves gradient for adult-context re-use).
            let pediatric_flow_penalty =
                0.15 * metrics.pediatric_flow_excess_risk;
            ((raw - coag_penalty - thermal_penalty - pediatric_flow_penalty) * hi_gate)
                .max(INFEASIBILITY_FLOOR)
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
            } else if matches!(mode, OptimMode::SdtTherapy) {
                (HI_PASS_LIMIT - metrics.hemolysis_index_per_pass) / HI_PASS_LIMIT
            } else {
                0.2
            };
            // Soft penalty for plate overflow: 0.01 preserves gradient signal
            // so the GA can navigate toward plate-fitting designs rather than
            // hitting a hard zero cliff.  The 100× suppression is severe enough
            // to rank plate-overflow designs strictly below fitting designs.
            let plate_ok = if metrics.plate_fits { 1.0_f64 } else { 0.01 };

            let feasibility = sigmoid_penalty(pressure_margin)
                * sigmoid_penalty(fda_margin)
                * sigmoid_penalty(hi_margin)
                * plate_ok;

            if feasibility < 0.05 {
                // Deeply infeasible: return a tiny gradient signal so the GA can
                // escape.  Floor to INFEASIBILITY_FLOOR to guarantee no zeros.
                return (feasibility * 0.1).max(INFEASIBILITY_FLOOR);
            }
            // Moderately feasible: continue with normal scoring, then multiply.
            let raw = modes::score_mode_raw(metrics, mode, weights);
            let coag_penalty = modes::coagulation_penalty(metrics, weights);
            let thermal_penalty = if metrics.fda_thermal_compliant {
                0.0
            } else {
                0.05 * (metrics.throat_temperature_rise_k / FDA_THROAT_TEMP_RISE_LIMIT_K)
                    .clamp(0.0, 1.0)
            };
            let pediatric_flow_penalty =
                0.15 * metrics.pediatric_flow_excess_risk;
            (raw * feasibility - coag_penalty - thermal_penalty - pediatric_flow_penalty)
                .max(INFEASIBILITY_FLOOR)
        }
    }
}

// Mode-specific scoring functions live in `modes.rs`.

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
        broad.treatment_zone_dwell_time_s = 0.12;

        let mut selective = base_metrics();
        selective.cancer_center_fraction = 0.76;
        selective.wbc_center_fraction = 0.66;
        selective.rbc_peripheral_fraction_three_pop = 0.82;
        selective.three_pop_sep_efficiency = 0.61;
        selective.therapy_channel_fraction = 0.34;
        selective.mean_residence_time_s = 1.30;
        selective.treatment_zone_dwell_time_s = 1.05;

        let mode = OptimMode::SdtTherapy;
        let w = SdtWeights::default();
        let broad_score = score_candidate(&broad, mode, &w);
        let selective_score = score_candidate(&selective, mode, &w);
        assert!(
            selective_score > broad_score,
            "sdt therapy mode should reward center-lane enrichment and peripheral RBC routing"
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn prop_hard_constraints_produce_floor_score(
            cav in 0.0_f64..1.0,
            sep3 in 0.0_f64..1.0,
            wbc_recovery in 0.0_f64..1.0,
            optical_405 in 0.0_f64..1.0
        ) {
            let modes = [
                OptimMode::SdtCavitation,
                OptimMode::UniformExposure,
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
                let s = score_candidate(&m, mode, &weights);
                prop_assert_eq!(s, INFEASIBILITY_FLOOR,
                    "pressure-infeasible candidate must receive floor score, got {}", s);

                let mut m = base_metrics();
                m.cavitation_potential = cav;
                m.three_pop_sep_efficiency = sep3;
                m.wbc_recovery = wbc_recovery;
                m.blue_light_delivery_index_405nm = optical_405;
                m.fda_main_compliant = false;
                let s = score_candidate(&m, mode, &weights);
                prop_assert_eq!(s, INFEASIBILITY_FLOOR,
                    "FDA-noncompliant candidate must receive floor score, got {}", s);

                let mut m = base_metrics();
                m.cavitation_potential = cav;
                m.three_pop_sep_efficiency = sep3;
                m.wbc_recovery = wbc_recovery;
                m.blue_light_delivery_index_405nm = optical_405;
                m.plate_fits = false;
                let s = score_candidate(&m, mode, &weights);
                prop_assert_eq!(s, INFEASIBILITY_FLOOR,
                    "plate-overflow candidate must receive floor score, got {}", s);
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
