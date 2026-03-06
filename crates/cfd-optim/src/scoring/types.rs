//! Optimisation-mode and scoring-weight types.

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
    /// - Three-population separation support (avoid simple-geometry bias)
    /// - 405 nm optical-delivery proxy (reward thin treatment channels)
    UniformExposure,

    /// Selective acoustic therapy on center-enriched treatment lanes.
    ///
    /// Uses the same selective routing objective as venturi-based therapy:
    /// cancer cells and WBCs should remain in the center treatment stream
    /// while RBCs are pushed into peripheral bypass channels. The treatment
    /// mechanism is externally applied ultrasound rather than venturi throats.
    SelectiveAcousticTherapy,

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
    /// - Separation efficiency (maximise |`x̃_cancer` - `x̃_healthy`|)
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

    /// Combined SDT therapy: cell separation, selective cavitation delivery,
    /// and haemolysis safety floor.  Intended as the primary clinical scoring
    /// mode that rewards topologies that place venturi throats on cancer-
    /// enriched treatment streams.
    ///
    /// Hard constraints: same as `SdtCavitation`.
    ///
    /// Soft objectives (revised weights):
    /// - 16% three-population separation efficiency (WBC/cancer center vs RBC wall)
    /// - 14% haemolysis minimisation (uses selective venturi composition)
    /// - 16% cavitation potential
    /// - 12% targeted cancer dose fraction (`cancer_dose_fraction`)
    /// - 8% RBC peripheral routing (`rbc_peripheral_fraction_three_pop`)
    /// - 6% throat-residence adequacy
    /// - 6% sonoluminescence proxy
    /// - 6% 405 nm optical-delivery proxy (Lumidox II support)
    /// - +5% synergy bonus when both selective dose and cavitation are strong
    /// - Additional gating: low-dose layouts (`dose < 0.05`) are down-weighted
    ///   to prevent non-targeted topologies from dominating SDT-therapy ranking
    ///
    /// This formulation prevents long-dwell low-shear layouts from outscoring
    /// true selective-cavitation designs where RBC exposure to venturi throats
    /// is intentionally reduced.
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

    /// Maximise **cancer-targeted hydrodynamic cavitation SDT**.
    ///
    /// This mode is the primary objective for millifluidic SDT designs that route
    /// cancer cells into a cavitating venturi center stream while shielding RBCs
    /// in peripheral bypass arms.  It combines:
    ///
    /// 1. Strong cavitation at the venturi throat (σ ≪ 1, narrow constriction).
    /// 2. High enrichment of cancer cells in the selected treatment arm of a primitive split tree.
    /// 3. RBC peripheral routing — minimising red cell exposure to cavitation.
    /// 4. High-energy bubble collapse (sonoluminescence proxy) for sonosensitiser
    ///    activation.
    ///
    /// Hard constraints: same as `SdtCavitation`.
    ///
    /// Soft objectives (configurable via [`SdtWeights`]):
    /// - **35%** cancer-targeted cavitation (`cancer_targeted_cavitation`)
    /// - **20%** three-population separation efficiency
    /// - **20%** RBC venturi protection (`rbc_venturi_protection`)
    /// - **15%** sonoluminescence proxy (collapse temperature)
    /// - **10%** WBC-targeted cavitation (`wbc_targeted_cavitation`)
    /// - Additional multiplicative optical gate from 405 nm blood attenuation proxy
    ///   (prefers thinner treatment channels for light-assisted protocols).
    /// - **+5%** synergy bonus when both `cancer_targeted_cavitation > 0.30`
    ///   and `three_pop_sep_efficiency > 0.30`.
    HydrodynamicCavitationSDT,

    /// Combined SDT + paediatric leukapheresis — addresses Milestone 12 joint-feasibility gap.
    ///
    /// Rewards topologies that simultaneously achieve:
    /// 1. WBC recovery, RBC removal, purity (leukapheresis objective).
    /// 2. Cancer-targeted hydrodynamic cavitation SDT.
    ///
    /// The two sub-scores are blended by `leuka_weight` and `sdt_weight`.
    /// An additional leukapheresis gate (WBC recovery + pediatric sub-score)
    /// suppresses SDT-only candidates in this combined mode.
    ///
    /// Hard constraints: `pressure_feasible`, `fda_main_compliant`, `plate_fits`.
    CombinedSdtLeukapheresis {
        /// Weight on the leukapheresis sub-score (default 0.5).
        leuka_weight: f64,
        /// Weight on the hydrodynamic SDT sub-score (default 0.5).
        sdt_weight: f64,
        /// Patient weight [kg]; used to check ECV ≤ 10% of total blood volume (neonate: 3.0).
        patient_weight_kg: f64,
    },

    /// Maximise therapeutic window: cancer-targeted SDT with bounded RBC lysis.
    ///
    /// Rewards designs that deliver strong cancer-cell cavitation treatment while
    /// minimising RBC haemolysis risk through selective flow routing and channel
    /// size asymmetry.  Uses the composite [`SdtMetrics::lysis_risk_index`] and
    /// [`SdtMetrics::therapeutic_window_score`] metrics.
    ///
    /// Hard constraints: `fda_main_compliant`, `pressure_feasible`, `plate_fits`.
    ///
    /// Soft objectives (configurable via [`SdtWeights`]):
    /// - **35%** therapeutic window score (`cancer_targeted_cav / lysis_risk`)
    /// - **25%** cancer-targeted cavitation (`cancer_targeted_cavitation`)
    /// - **20%** low lysis risk (`1 − lysis_risk_index × 1000`, clamped [0, 1])
    /// - **10%** full FDA compliance bonus (`fda_overall_compliant` 1.0 / 0.5)
    /// - **10%** therapy zone channel utilisation (`therapy_channel_fraction`)
    /// - **+5%** synergy when `cancer_targeted_cavitation > 0.25`
    ///   **and** `lysis_risk_index < 0.001`
    RbcProtectedSdt,
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

    // ── Coagulation / platelet-activation penalty (all modes) ───────────
    /// Penalty weight for platelet activation index (PAI).
    ///
    /// Applied after the mode-specific score as:
    /// `score -= coag_weight × clamp(PAI / PAI_PASS_LIMIT, 0, 1)`
    ///
    /// Default: 0.10.  Set to 0.0 to disable the penalty.
    pub coag_weight: f64,

    // ── HydrodynamicCavitationSDT mode weights ───────────────────────────
    /// Weight on cancer-targeted cavitation (`cancer_center_fraction × cavitation_intensity`).
    pub hydro_cancer_cav: f64,
    /// Weight on three-population separation efficiency in `HydroSDT` mode.
    pub hydro_sep3: f64,
    /// Weight on RBC venturi protection (`rbc_peripheral_fraction × (1 − cav_intensity × rbc_exposure)`).
    pub hydro_rbc_protection: f64,
    /// Weight on sonoluminescence proxy (adiabatic collapse temperature).
    pub hydro_sonolum: f64,
    /// Weight on WBC-targeted cavitation (`wbc_center_fraction × cavitation_intensity`).
    pub hydro_wbc_cav: f64,

    // ── RbcProtectedSdt mode weights ─────────────────────────────────────
    /// Weight on therapeutic window score
    /// (`cancer_targeted_cavitation / lysis_risk_index / THERAPEUTIC_WINDOW_REF`).
    pub rbc_protected_window: f64,
    /// Weight on cancer-targeted cavitation in `RbcProtectedSdt` mode.
    pub rbc_protected_cav: f64,
    /// Weight on lysis risk penalty (`1 − lysis_risk_index × 1000`, clamped [0, 1]).
    pub rbc_protected_lysis: f64,
    /// Weight on the overall FDA compliance bonus (1.0 = fully compliant, 0.5 = throat exception).
    pub rbc_protected_fda: f64,
    /// Weight on therapy channel fraction (fraction of chip in active therapy zone).
    pub rbc_protected_coverage: f64,
}

impl Default for SdtWeights {
    fn default() -> Self {
        Self {
            cav_potential: 0.60,
            cav_hemolysis: 0.20,
            cav_coverage: 0.20,

            exp_uniformity: 0.20,
            exp_coverage: 0.45,
            exp_residence: 0.35,

            sep_efficiency: 0.50,
            sep_cavitation: 0.30,
            sep_hemolysis: 0.20,

            sep3_efficiency: 0.40,
            sep3_wbc_center: 0.20,
            sep3_cavitation: 0.25,
            sep3_hemolysis: 0.15,

            coag_weight: 0.10,

            hydro_cancer_cav: 0.35,
            hydro_sep3: 0.20,
            hydro_rbc_protection: 0.20,
            hydro_sonolum: 0.15,
            hydro_wbc_cav: 0.10,

            rbc_protected_window: 0.35,
            rbc_protected_cav: 0.25,
            rbc_protected_lysis: 0.20,
            rbc_protected_fda: 0.10,
            rbc_protected_coverage: 0.10,
        }
    }
}

// ── Constraint mode ───────────────────────────────────────────────────────────

/// Controls how physical-feasibility constraints are applied during scoring.
///
/// * [`ScoreMode::HardConstraint`] — disqualifies infeasible candidates with score `0.0`.
/// * [`ScoreMode::SmoothPenalty`] — replaces the hard `0.0` cliff with a smooth sigmoid
///   gradient so the genetic algorithm can navigate toward the feasible region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum ScoreMode {
    /// Return `0.0` for any infeasible candidate. Default.
    #[default]
    HardConstraint,
    /// Apply a smooth feasibility multiplier instead of a hard kill.
    SmoothPenalty,
}
