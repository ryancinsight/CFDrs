use super::modes;
use super::sigmoid_penalty;
use super::{INFEASIBILITY_SCORE, OptimMode, ScoreMode, SdtWeights};
use crate::constraints::{FDA_THROAT_TEMP_RISE_LIMIT_K, HI_PASS_LIMIT, THERAPEUTIC_HI_PASS_LIMIT};
use crate::metrics::SdtMetrics;

// ── Score functions ──────────────────────────────────────────────────────────

/// Compute the score for a single candidate given a mode and weights.
///
/// Returns a value in **[0.0, 1.0]**.  Higher scores indicate better designs.
/// Infeasible candidates (hard-constraint violations) receive an exact score
/// of `0.0`; the smooth penalty mode provides the non-zero gradient variant.
#[must_use]
pub fn score_candidate(metrics: &SdtMetrics, mode: OptimMode, weights: &SdtWeights) -> f64 {
    score_candidate_impl(metrics, mode, weights, ScoreMode::HardConstraint, 0.0)
}

pub(super) fn score_candidate_impl(
    metrics: &SdtMetrics,
    mode: OptimMode,
    weights: &SdtWeights,
    constraint_mode: ScoreMode,
    inlet_gauge_pa: f64,
) -> f64 {
    // ── Feasibility check ─────────────────────────────────────────────────
    match constraint_mode {
        ScoreMode::HardConstraint => {
            if !metrics.pressure_feasible || !metrics.fda_main_compliant || !metrics.plate_fits {
                return INFEASIBILITY_SCORE;
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
            let pediatric_flow_penalty = 0.15 * metrics.pediatric_flow_excess_risk;
            // Channel overlap penalty: the 1D independent-resistance model
            // loses accuracy when channels physically merge.  However, when
            // overlapping channels have different widths (ratio > 2), the
            // velocity mismatch at the merge zone creates passive inertial
            // cell sorting that may enhance separation.  The penalty is
            // therefore reduced for high width ratios (potential sorting
            // benefit) and full strength for symmetric merges (ratio ~ 1).
            let overlap_raw = ((metrics.channel_overlap_fraction - 0.3) / 0.7).clamp(0.0, 1.0);
            // Width-ratio discount: ratio 1.0 -> factor 1.0 (full penalty);
            // ratio 3.0+ -> factor 0.3 (70% discount for sorting benefit).
            let ratio_discount = if metrics.overlap_width_ratio > 1.0 {
                (1.0 / metrics.overlap_width_ratio.min(3.0)).max(0.3)
            } else {
                1.0
            };
            let overlap_penalty = 0.05 * overlap_raw * ratio_discount;
            ((raw - coag_penalty - thermal_penalty - pediatric_flow_penalty - overlap_penalty)
                * hi_gate)
                .max(INFEASIBILITY_SCORE)
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
            // Plate overflow remains an exact zero-score condition because the
            // physical design cannot be fabricated within the plate boundary.
            let plate_ok = if metrics.plate_fits { 1.0_f64 } else { 0.0 };

            let feasibility = sigmoid_penalty(pressure_margin)
                * sigmoid_penalty(fda_margin)
                * sigmoid_penalty(hi_margin)
                * plate_ok;

            if feasibility < 0.05 {
                // Deeply infeasible: return a small feasibility-proportional signal
                // so the GA can escape while keeping the penalty monotone.
                return feasibility * 0.1;
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
            let pediatric_flow_penalty = 0.15 * metrics.pediatric_flow_excess_risk;
            let overlap_raw = ((metrics.channel_overlap_fraction - 0.3) / 0.7).clamp(0.0, 1.0);
            let ratio_discount = if metrics.overlap_width_ratio > 1.0 {
                (1.0 / metrics.overlap_width_ratio.min(3.0)).max(0.3)
            } else {
                1.0
            };
            let overlap_penalty = 0.05 * overlap_raw * ratio_discount;
            (raw * feasibility
                - coag_penalty
                - thermal_penalty
                - pediatric_flow_penalty
                - overlap_penalty)
                .max(INFEASIBILITY_SCORE)
        }
    }
}

// Mode-specific scoring functions live in `modes.rs`.
