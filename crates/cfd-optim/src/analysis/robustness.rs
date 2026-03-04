//! Design robustness analysis via parametric perturbation sweeps.
//!
//! Evaluates how sensitive a candidate's score is to ±10 % / ±20 % changes in
//! key operating parameters (`flow_rate_m3_s`, `inlet_gauge_pa`,
//! `throat_diameter_m`).  The coefficient of variation (CV = σ/μ) across all
//! perturbed evaluations quantifies overall robustness.
//!
//! A design is considered *robust* when `CV < 0.10` (< 10 % relative variation).
//!
//! # Usage
//! ```ignore
//! use cfd_optim::analysis::robustness::{robustness_sweep, STANDARD_PERTURBATIONS};
//! use cfd_optim::scoring::{OptimMode, SdtWeights};
//!
//! let report = robustness_sweep(&candidate, OptimMode::RbcProtectedSdt,
//!                               &SdtWeights::default(),
//!                               &STANDARD_PERTURBATIONS);
//! println!("Robust: {}", report.is_robust);
//! ```

use serde::{Deserialize, Serialize};

use crate::{
    design::DesignCandidate,
    metrics::compute_metrics,
    scoring::{score_candidate, OptimMode, SdtWeights},
};

// ── Constants ─────────────────────────────────────────────────────────────────

/// Default perturbation fractions: ±10 % and ±20 %.
pub const STANDARD_PERTURBATIONS: [f64; 4] = [-0.20, -0.10, 0.10, 0.20];

/// CV threshold below which a design is classified as "robust".
pub const ROBUSTNESS_CV_THRESHOLD: f64 = 0.10;

// ── Output types ──────────────────────────────────────────────────────────────

/// Robustness summary for a single design candidate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RobustnessReport {
    /// Identifier of the evaluated candidate.
    pub candidate_id: String,
    /// Score under nominal (unperturbed) parameters.
    pub score_nominal: f64,
    /// Minimum score across all perturbations.
    pub score_min: f64,
    /// Maximum score across all perturbations.
    pub score_max: f64,
    /// Coefficient of variation `σ/μ` across all perturbed evaluations.
    ///
    /// Zero when the mean score is essentially zero (`< 1e-9`).
    pub score_cv: f64,
    /// Label of the perturbation that produced the worst (lowest) score.
    ///
    /// Format: `"flow_rate+20%"`, `"inlet_pressure-10%"`, `"throat_diameter+20%"`, or `"nominal"`.
    pub worst_case_param: String,
    /// `true` if `score_cv < ROBUSTNESS_CV_THRESHOLD` (< 10 % relative variation).
    pub is_robust: bool,
}

// ── Core function ─────────────────────────────────────────────────────────────

/// Run parametric perturbation sweeps and produce a [`RobustnessReport`].
///
/// Perturbs `flow_rate_m3_s`, `inlet_gauge_pa`, and (if the topology has a
/// venturi) `throat_diameter_m` by each fraction in `perturbation_fracs`.
///
/// `perturbation_fracs`: relative perturbations to apply, e.g.
/// `&[-0.20, -0.10, 0.10, 0.20]`.  Use [`STANDARD_PERTURBATIONS`] for the
/// canonical ±10 %/±20 % sweep.
///
/// # Returns
///
/// [`RobustnessReport`] summarising nominal score, min/max, CV, and the
/// worst-case perturbed parameter.
#[must_use]
pub fn robustness_sweep(
    candidate: &DesignCandidate,
    mode: OptimMode,
    weights: &SdtWeights,
    perturbation_fracs: &[f64],
) -> RobustnessReport {
    let score_nominal = compute_metrics(candidate)
        .map(|m| score_candidate(&m, mode, weights))
        .unwrap_or(0.0);

    let mut scored: Vec<(f64, String)> = vec![(score_nominal, "nominal".into())];

    // ── Flow rate perturbations ───────────────────────────────────────────────
    for &frac in perturbation_fracs {
        let mut c = candidate.clone();
        c.flow_rate_m3_s *= 1.0 + frac;
        let s = compute_metrics(&c)
            .map(|m| score_candidate(&m, mode, weights))
            .unwrap_or(0.0);
        scored.push((s, format!("flow_rate{:+.0}%", frac * 100.0)));
    }

    // ── Inlet pressure perturbations ─────────────────────────────────────────
    for &frac in perturbation_fracs {
        let mut c = candidate.clone();
        c.inlet_gauge_pa *= 1.0 + frac;
        let s = compute_metrics(&c)
            .map(|m| score_candidate(&m, mode, weights))
            .unwrap_or(0.0);
        scored.push((s, format!("inlet_pressure{:+.0}%", frac * 100.0)));
    }

    // ── Throat diameter perturbations (venturi topologies only) ──────────────
    if candidate.topology.has_venturi() {
        for &frac in perturbation_fracs {
            let mut c = candidate.clone();
            c.throat_diameter_m *= 1.0 + frac;
            let s = compute_metrics(&c)
                .map(|m| score_candidate(&m, mode, weights))
                .unwrap_or(0.0);
            scored.push((s, format!("throat_diameter{:+.0}%", frac * 100.0)));
        }
    }

    // ── Statistics ────────────────────────────────────────────────────────────
    let vals: Vec<f64> = scored.iter().map(|(s, _)| *s).collect();
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let variance = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();
    let score_cv = if mean.abs() > 1e-9 {
        std_dev / mean
    } else {
        0.0
    };

    let (score_min, worst_label) = scored
        .iter()
        .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .cloned()
        .unwrap_or((0.0, "none".into()));

    let score_max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    RobustnessReport {
        candidate_id: candidate.id.clone(),
        score_nominal,
        score_min,
        score_max,
        score_cv,
        worst_case_param: worst_label,
        is_robust: score_cv < ROBUSTNESS_CV_THRESHOLD,
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::{CrossSectionShape, DesignCandidate, DesignTopology, TreatmentZoneMode};
    use crate::scoring::{OptimMode, SdtWeights};

    /// Minimal candidate with a venturi — used to verify perturbation logic.
    fn venturi_candidate() -> DesignCandidate {
        DesignCandidate {
            id: "rob-test-sv".into(),
            topology: DesignTopology::SingleVenturi,
            flow_rate_m3_s: 5.0e-6, // 300 mL/min
            inlet_gauge_pa: 300_000.0,
            throat_diameter_m: 100e-6,
            inlet_diameter_m: 4e-3,
            throat_length_m: 200e-6,
            channel_width_m: 4e-3,
            channel_height_m: 1e-3,
            serpentine_segments: 6,
            segment_length_m: 0.045,
            bend_radius_m: 4.5e-3,
            feed_hematocrit: 0.45,
            trifurcation_center_frac: 1.0 / 3.0,
            cif_pretri_center_frac: 1.0 / 3.0,
            cif_terminal_tri_center_frac: 1.0 / 3.0,
            cif_terminal_bi_treat_frac: 0.68,
            asymmetric_narrow_frac: 0.5,
            trifurcation_left_frac: 1.0 / 3.0,
            cross_section_shape: CrossSectionShape::Rectangular,
            treatment_zone_mode: TreatmentZoneMode::VenturiThroats,
            centerline_venturi_throat_count: 1,
        }
    }

    #[test]
    fn report_returns_candidate_id() {
        let c = venturi_candidate();
        let w = SdtWeights::default();
        let report = robustness_sweep(&c, OptimMode::RbcProtectedSdt, &w, &STANDARD_PERTURBATIONS);
        assert_eq!(report.candidate_id, "rob-test-sv");
    }

    #[test]
    fn score_min_le_nominal_le_max() {
        let c = venturi_candidate();
        let w = SdtWeights::default();
        let report = robustness_sweep(&c, OptimMode::RbcProtectedSdt, &w, &STANDARD_PERTURBATIONS);
        // score_min ≤ score_nominal ≤ score_max (allowing floating-point ties)
        assert!(
            report.score_min <= report.score_nominal + 1e-9,
            "score_min {} > score_nominal {}",
            report.score_min,
            report.score_nominal
        );
        assert!(
            report.score_nominal <= report.score_max + 1e-9,
            "score_nominal {} > score_max {}",
            report.score_nominal,
            report.score_max
        );
    }

    #[test]
    fn cv_is_non_negative() {
        let c = venturi_candidate();
        let w = SdtWeights::default();
        let report = robustness_sweep(&c, OptimMode::RbcProtectedSdt, &w, &STANDARD_PERTURBATIONS);
        assert!(
            report.score_cv >= 0.0,
            "CV must be non-negative: {}",
            report.score_cv
        );
    }

    #[test]
    fn empty_perturbations_gives_nominal_only() {
        let c = venturi_candidate();
        let w = SdtWeights::default();
        let report = robustness_sweep(&c, OptimMode::RbcProtectedSdt, &w, &[]);
        // With only nominal, CV = 0
        assert_eq!(report.score_cv, 0.0);
        assert_eq!(report.score_min, report.score_nominal);
        assert_eq!(report.score_max, report.score_nominal);
    }

    #[test]
    fn is_robust_flag_reflects_cv_threshold() {
        let c = venturi_candidate();
        let w = SdtWeights::default();
        let report = robustness_sweep(&c, OptimMode::RbcProtectedSdt, &w, &STANDARD_PERTURBATIONS);
        // Verify the flag matches the threshold
        assert_eq!(report.is_robust, report.score_cv < ROBUSTNESS_CV_THRESHOLD);
    }
}
