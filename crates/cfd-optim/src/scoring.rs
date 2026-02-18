//! Multi-objective scoring functions for SDT design candidates.
//!
//! Two optimisation modes are provided:
//!
//! - **[`OptimMode::SdtCavitation`]** — maximise cavitation intensity at the
//!   venturi throat while maintaining FDA haemolysis compliance in the main
//!   channels.  Designed for **hydrodynamic cavitation sonodynamic therapy**
//!   where the mechanical energy of collapsing bubbles drives the therapeutic
//!   effect.
//!
//! - **[`OptimMode::UniformExposure`]** — maximise spatial uniformity of flow
//!   across all 36 treatment wells and maximise residence time in the exposure
//!   zone.  Designed for **light / ultrasound** modalities where every well
//!   must receive equal treatment duration.
//!
//! - **[`OptimMode::Combined`]** — weighted combination of both objectives,
//!   useful for devices that couple cavitation initiation with uniform
//!   distribution (e.g. `VenturiSerpentine` topologies).
//!
//! Any candidate that violates either hard constraint — non-feasible pressure
//! drop **or** FDA main-channel shear exceedance — receives a score of `0.0`.

use crate::constraints::HI_PASS_LIMIT;
use crate::metrics::SdtMetrics;
use serde::{Deserialize, Serialize};

// ── Optimisation mode ────────────────────────────────────────────────────────

/// The optimisation objective to use when ranking design candidates.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum OptimMode {
    /// Maximise hydrodynamic cavitation potential at the venturi throat.
    ///
    /// Hard constraints (score = 0 if violated):
    /// - Total ΔP ≤ available inlet gauge pressure
    /// - Main-channel shear ≤ 150 Pa
    ///
    /// Soft objectives (weighted sum):
    /// - Cavitation potential (σ < 1, higher potential = lower σ)
    /// - Low haemolysis index per pass
    /// - Well coverage fraction (cavitation should reach tissue)
    SdtCavitation,

    /// Maximise uniform light / ultrasound exposure across the 6×6 well zone.
    ///
    /// Hard constraints: same as `SdtCavitation`.
    ///
    /// Soft objectives:
    /// - Flow uniformity at outlets
    /// - Well coverage fraction
    /// - Residence time in treatment zone
    UniformExposure,

    /// Weighted combination of cavitation and exposure objectives.
    Combined {
        /// Weight on the cavitation sub-score (0.0–1.0).
        cavitation_weight: f64,
        /// Weight on the exposure sub-score (0.0–1.0).
        exposure_weight: f64,
    },
}

// ── Scoring weights ──────────────────────────────────────────────────────────

/// Per-objective scoring weights.
///
/// All weight groups must sum to ≤ 1.0; the remainder is unused.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SdtWeights {
    // ── Cavitation mode weights ──────────────────────────────────────────
    /// Weight on the cavitation potential (σ < 1, maximise 1 − σ).
    pub cav_potential: f64,
    /// Weight on low haemolysis index (penalise HI > threshold).
    pub cav_hemolysis: f64,
    /// Weight on well coverage fraction in cavitation mode.
    pub cav_coverage: f64,

    // ── Exposure mode weights ────────────────────────────────────────────
    /// Weight on outlet flow uniformity (CV of outlet flows).
    pub exp_uniformity: f64,
    /// Weight on well coverage fraction in exposure mode.
    pub exp_coverage: f64,
    /// Weight on normalised residence time in the treatment zone.
    pub exp_residence: f64,
}

impl Default for SdtWeights {
    fn default() -> Self {
        Self {
            // Cavitation: emphasis on achieving cavitation; haemolysis is a
            // secondary concern since throat exposure is brief.
            cav_potential: 0.60,
            cav_hemolysis: 0.20,
            cav_coverage: 0.20,

            // Exposure: uniformity and coverage are equally important;
            // residence time is a useful secondary criterion.
            exp_uniformity: 0.45,
            exp_coverage: 0.35,
            exp_residence: 0.20,
        }
    }
}

// ── Score functions ──────────────────────────────────────────────────────────

/// Compute the score for a single candidate given a mode and weights.
///
/// Returns a value in **[0.0, 1.0]**.  Higher scores indicate better designs.
/// Returns `0.0` for any candidate that violates hard constraints.
pub fn score_candidate(metrics: &SdtMetrics, mode: OptimMode, weights: &SdtWeights) -> f64 {
    // ── Hard constraints (disqualify) ──────────────────────────────────────
    if !metrics.pressure_feasible || !metrics.fda_main_compliant {
        return 0.0;
    }

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
    }
}

/// Score for the SDT cavitation objective.
fn score_cavitation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    // Cavitation potential is already in [0, 1]; σ → 0 gives potential → 1.
    let cav = metrics.cavitation_potential;

    // Haemolysis penalty: normalise HI against acceptable threshold.
    // hi_factor = 1.0 when HI = 0; 0.5 when HI = threshold; 0 when HI = 2× threshold.
    let hi_ratio = (metrics.hemolysis_index_per_pass / HI_PASS_LIMIT).min(2.0);
    let hi_factor = (1.0 - 0.5 * hi_ratio).max(0.0);

    // Coverage of the treatment zone
    let coverage = metrics.well_coverage_fraction.clamp(0.0, 1.0);

    w.cav_potential * cav + w.cav_hemolysis * hi_factor + w.cav_coverage * coverage
}

/// Score for the uniform exposure objective.
fn score_exposure(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let uniformity = metrics.flow_uniformity.clamp(0.0, 1.0);
    let coverage = metrics.well_coverage_fraction.clamp(0.0, 1.0);

    // Normalise residence time on a logarithmic scale with a 30-second target.
    // This prevents extremely long paths (and huge ΔP) from dominating.
    let target_s = 30.0_f64;
    let res_norm = if metrics.mean_residence_time_s > 0.0 {
        (metrics.mean_residence_time_s.ln() - 1e-3_f64.ln())
            / (target_s.ln() - 1e-3_f64.ln())
    } else {
        0.0
    }
    .clamp(0.0, 1.0);

    w.exp_uniformity * uniformity + w.exp_coverage * coverage + w.exp_residence * res_norm
}

// ── Utility ──────────────────────────────────────────────────────────────────

/// Return a human-readable summary line for a scored candidate.
pub fn score_description(mode: OptimMode) -> &'static str {
    match mode {
        OptimMode::SdtCavitation => "SDT Cavitation",
        OptimMode::UniformExposure => "Uniform Exposure",
        OptimMode::Combined { .. } => "Combined (Cavitation + Exposure)",
    }
}
