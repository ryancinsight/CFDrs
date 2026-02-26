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

    /// Maximise separation efficiency of cancer cells from healthy cells
    /// while maintaining cavitation and FDA compliance.
    ///
    /// Hard constraints: same as `SdtCavitation`.
    ///
    /// Soft objectives:
    /// - Separation efficiency (maximise |x̃_cancer - x̃_healthy|)
    /// - Cavitation potential (secondary goal)
    /// - Low haemolysis index
    CellSeparation,

    /// Maximise three-population separation (WBC + cancer → center, RBC → periphery)
    /// while delivering SDT cavitation through the venturi throat.
    ///
    /// Hard constraints: same as `SdtCavitation`.
    ///
    /// Soft objectives:
    /// - Three-population separation efficiency (`x̃_rbc − max(x̃_cancer, x̃_wbc)`)
    /// - WBC center fraction (WBCs co-focused with cancer cells)
    /// - Cavitation potential (SDT delivery in center channel)
    /// - Low haemolysis index
    ThreePopSeparation,

    /// Combined SDT therapy: cell separation, haemolysis minimisation, cavitation
    /// potential, and flow uniformity.  Intended as the primary clinical scoring
    /// mode that balances all therapy-relevant objectives.
    ///
    /// Hard constraints: same as `SdtCavitation`.
    ///
    /// Soft objectives (fixed weights):
    /// - 35% three-population separation efficiency
    /// - 45% haemolysis minimisation (`1 − HI/HI_limit`, clamped to [0, 1])
    /// - 20% cavitation potential
    ///
    /// Note: flow_uniformity is 1.0 by construction for all symmetric topologies;
    /// its former 15% weight has been folded into haemolysis minimisation.
    SdtTherapy,

    /// Optimise for **paediatric leukapheresis**:
    /// WBC recovery + RBC removal + WBC purity + throughput feasibility.
    ///
    /// Hard constraints:
    /// - `pressure_feasible` — device can be driven at the available gauge pressure.
    /// - Wall shear ≤ 150 Pa (`fda_main_compliant`) — protect fragile neonatal cells.
    ///
    /// Soft objectives (fixed weights):
    /// - 40% WBC recovery (`wbc_recovery`)
    /// - 30% RBC removal (`1 − rbc_pass_fraction`)
    /// - 20% WBC purity (`wbc_purity`)
    /// - 10% throughput feasibility (`min(1.0, Q_total_mL_min / 10.0)`)
    ///
    /// Design for neonates: `patient_weight_kg` is stored for downstream ECV
    /// constraint checking (`ECV ≤ 10% × patient_blood_volume`).
    /// `patient_blood_volume_ml ≈ patient_weight_kg × 85`.
    PediatricLeukapheresis {
        /// Patient weight [kg]; used to check ECV ≤ 10% of total blood volume.
        patient_weight_kg: f64,
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

    // ── Cell separation weights (2-population: cancer vs RBC) ───────────
    /// Weight on separation efficiency.
    pub sep_efficiency: f64,
    /// Weight on cavitation potential in separation mode.
    pub sep_cavitation: f64,
    /// Weight on low haemolysis index in separation mode.
    pub sep_hemolysis: f64,

    // ── Three-population separation weights (WBC+cancer→center, RBC→wall) ─
    /// Weight on three-population separation efficiency (`x̃_rbc − max(x̃_cancer, x̃_wbc)`).
    pub sep3_efficiency: f64,
    /// Weight on WBC center fraction (WBCs co-focused with cancer cells).
    pub sep3_wbc_center: f64,
    /// Weight on cavitation potential in three-population mode.
    pub sep3_cavitation: f64,
    /// Weight on low haemolysis index in three-population mode.
    pub sep3_hemolysis: f64,
}

impl Default for SdtWeights {
    fn default() -> Self {
        Self {
            // Cavitation: emphasis on achieving cavitation; haemolysis is a
            // secondary concern since throat exposure is brief.
            cav_potential: 0.60,
            cav_hemolysis: 0.20,
            cav_coverage: 0.20,

            // Exposure: flow_uniformity is always 1.0 by construction for all
            // current symmetric topologies (balanced binary trees + serpentines).
            // Assigning weight to a constant term wastes discriminating power.
            // Redistribute the 0.45 uniformity weight to coverage and residence.
            exp_uniformity: 0.0,
            exp_coverage: 0.55,
            exp_residence: 0.45,

            // Cell Separation: primary goal is separation; cavitation is secondary
            // but still required for therapy.
            sep_efficiency: 0.50,
            sep_cavitation: 0.30,
            sep_hemolysis: 0.20,

            // Three-population separation: efficiency and WBC co-focus are primary;
            // cavitation (SDT delivery to co-focused cells) and HI are secondary.
            sep3_efficiency: 0.40,
            sep3_wbc_center: 0.20,
            sep3_cavitation: 0.25,
            sep3_hemolysis: 0.15,
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
    // `plate_fits` rejects ParallelMicrochannelArray candidates whose
    // n_channels × 2×width footprint exceeds the 96-well plate height (85.47 mm).
    // This eliminates the degenerate GA solution where n=500 channels at 6 mm width
    // (500×12 mm = 6000 mm) scores perfectly despite being physically impossible.
    if !metrics.pressure_feasible || !metrics.fda_main_compliant || !metrics.plate_fits {
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
        OptimMode::CellSeparation => score_cell_separation(metrics, weights),
        OptimMode::ThreePopSeparation => score_three_pop_separation(metrics, weights),
        OptimMode::SdtTherapy => score_sdt_therapy(metrics),
        OptimMode::PediatricLeukapheresis { patient_weight_kg } => {
            score_pediatric_leukapheresis(metrics, patient_weight_kg)
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

/// Score for the cell separation objective.
fn score_cell_separation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    let separation = metrics.cell_separation_efficiency.clamp(0.0, 1.0);
    let cav = metrics.cavitation_potential;

    // Haemolysis penalty (similar to cavitation mode)
    let hi_ratio = (metrics.hemolysis_index_per_pass / HI_PASS_LIMIT).min(2.0);
    let hi_factor = (1.0 - 0.5 * hi_ratio).max(0.0);

    w.sep_efficiency * separation + w.sep_cavitation * cav + w.sep_hemolysis * hi_factor
}

/// Score for the three-population separation objective (WBC+cancer→center, RBC→wall).
///
/// Rewards designs where:
/// - `three_pop_sep_efficiency` is high (RBCs pushed to wall relative to both WBC and cancer)
/// - `wbc_center_fraction` is high (WBCs co-focus with cancer cells in center channel)
/// - Cavitation potential is non-zero (SDT is delivered to the co-focused cells)
/// - Haemolysis index is low
fn score_three_pop_separation(metrics: &SdtMetrics, w: &SdtWeights) -> f64 {
    // Clamp sep_efficiency to [0, 1]; negative efficiency (cells mixed up) scores 0.
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);
    let wbc_center = metrics.wbc_center_fraction.clamp(0.0, 1.0);
    let cav = metrics.cavitation_potential;

    let hi_ratio = (metrics.hemolysis_index_per_pass / HI_PASS_LIMIT).min(2.0);
    let hi_factor = (1.0 - 0.5 * hi_ratio).max(0.0);

    w.sep3_efficiency * sep3
        + w.sep3_wbc_center * wbc_center
        + w.sep3_cavitation * cav
        + w.sep3_hemolysis * hi_factor
}

/// Score for the combined SDT therapy objective.
///
/// Fixed weights: 35% separation, 45% haemolysis, 20% cavitation.
///
/// Note: flow_uniformity is 1.0 by construction for all current symmetric
/// topologies, so its former 15% weight has been folded into haemolysis
/// (the next most clinically critical differentiating objective).
fn score_sdt_therapy(metrics: &SdtMetrics) -> f64 {
    // Three-population separation efficiency: clamp negative (bad mixing) to 0.
    let sep3 = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);

    // Haemolysis minimisation: 1 at HI=0, 0 at HI=HI_limit.
    let hi_score = (1.0 - metrics.hemolysis_index_per_pass / HI_PASS_LIMIT).clamp(0.0, 1.0);

    // Cavitation potential (already in [0, 1]).
    let cav = metrics.cavitation_potential;

    0.35 * sep3 + 0.45 * hi_score + 0.20 * cav
}

/// Score for the paediatric leukapheresis objective.
///
/// Weights: 40% WBC recovery + 30% RBC removal + 20% WBC purity + 10% throughput.
/// Additionally penalises designs whose ECV exceeds 10% of the patient's blood volume
/// (`patient_blood_volume_ml ≈ patient_weight_kg × 85 mL/kg`).
fn score_pediatric_leukapheresis(metrics: &SdtMetrics, patient_weight_kg: f64) -> f64 {
    // Primary objectives
    let wbc_rec = metrics.wbc_recovery.clamp(0.0, 1.0);
    let rbc_rem = (1.0 - metrics.rbc_pass_fraction).clamp(0.0, 1.0);
    let purity  = metrics.wbc_purity.clamp(0.0, 1.0);

    // Throughput score: normalised to 10 mL/min clinical leukapheresis target.
    // flow_rate_ml_min is stored in SdtMetrics (m³/s × 60 × 1e6).
    // A single chip at 1 mL/min scores 0.1; six chips in parallel reach ~6 mL/min (0.6).
    // Ten chips at 1 mL/min = 10 mL/min → saturates to 1.0.
    let throughput_score = (metrics.flow_rate_ml_min / 10.0).min(1.0);

    // ECV penalty: if ECV > 10% of patient blood volume, apply a fractional penalty.
    let patient_bv_ml = patient_weight_kg * 85.0; // 85 mL/kg neonatal blood volume
    let max_ecv_ml    = patient_bv_ml * 0.10;
    let ecv_ok_factor = if metrics.total_ecv_ml > 0.0 && metrics.total_ecv_ml <= max_ecv_ml {
        1.0
    } else if metrics.total_ecv_ml > max_ecv_ml {
        // Graded penalty: score × (max_ecv / actual_ecv)
        (max_ecv_ml / metrics.total_ecv_ml).clamp(0.0, 1.0)
    } else {
        1.0 // total_ecv_ml = 0 means non-leukapheresis topology — no ECV constraint
    };

    let raw = 0.40 * wbc_rec + 0.30 * rbc_rem + 0.20 * purity + 0.10 * throughput_score;
    raw * ecv_ok_factor
}

// ── Utility ──────────────────────────────────────────────────────────────────

/// Return a human-readable summary line for a scored candidate.
pub fn score_description(mode: OptimMode) -> &'static str {
    match mode {
        OptimMode::SdtCavitation => "SDT Cavitation",
        OptimMode::UniformExposure => "Uniform Exposure",
        OptimMode::Combined { .. } => "Combined (Cavitation + Exposure)",
        OptimMode::CellSeparation => "Cell Separation + SDT",
        OptimMode::ThreePopSeparation => "Three-Pop Separation (WBC+Cancer→Center, RBC→Wall) + SDT",
        OptimMode::SdtTherapy => "SDT Therapy (Sep + HI + Cav + Uniformity)",
        OptimMode::PediatricLeukapheresis { .. } => {
            "Paediatric Leukapheresis (WBC Recovery + RBC Removal + Purity + ECV)"
        }
    }
}
