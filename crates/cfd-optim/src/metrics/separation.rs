//! Cell-separation physics helpers for three-population and leukapheresis models.

use cfd_core::physics::fluid::blood::CassonBlood;

use crate::constraints::{BLOOD_DENSITY_KG_M3, TREATMENT_HEIGHT_MM};
use crate::design::{DesignCandidate, DesignTopology};

// ── Three-population separation ──────────────────────────────────────────────

/// Intermediate results from the three-population inertial separation model.
///
/// Computed by [`three_population_separation`] for topologies with curved or
/// high-Dean-number channel sections.
#[derive(Debug, Clone)]
pub(super) struct ThreePopMetrics {
    /// `x̃_rbc − max(x̃_cancer, x̃_wbc)` clamped to `[−1, 1]`.
    /// Positive = RBCs are more peripheral than both WBCs and cancer cells.
    pub sep_efficiency: f64,
    /// Fraction of WBCs focused into the center channel (`x̃ < 0.3`).
    pub wbc_center_fraction: f64,
    /// Fraction of cancer cells focused into the center channel (`x̃ < 0.3`).
    pub cancer_center_fraction: f64,
    /// Fraction of RBCs in the peripheral channels (`x̃ > 0.3`).
    pub rbc_peripheral_fraction: f64,
    /// WBC normalised equilibrium position `x̃ ∈ [0, 1]`.
    pub wbc_eq_pos: f64,
    /// Cancer cell normalised equilibrium position `x̃ ∈ [0, 1]`.
    pub cancer_eq_pos: f64,
    /// RBC normalised equilibrium position `x̃ ∈ [0, 1]`.
    pub rbc_eq_pos: f64,
}

impl Default for ThreePopMetrics {
    fn default() -> Self {
        Self {
            sep_efficiency: 0.0,
            wbc_center_fraction: 0.0,
            cancer_center_fraction: 0.0,
            rbc_peripheral_fraction: 0.0,
            wbc_eq_pos: 0.5,
            cancer_eq_pos: 0.5,
            rbc_eq_pos: 0.5,
        }
    }
}

/// Compute three-population inertial focusing equilibrium positions for
/// WBC + cancer → center and RBC → periphery.
///
/// Delegates to [`cfd_1d::three_population_equilibria`] which evaluates all
/// three cell types simultaneously under identical flow conditions (Di Carlo 2009
/// lift + Gossett 2009 Dean drag), with CFL correction for WBCs (Fedosov 2012).
///
/// The κ > 0.07 inertial-focusing gate is NOT applied: the underlying solver
/// always returns a valid force-balance result via bisection, valid across all Re.
/// Straight venturi channels pass `bend_radius_m = None` (no Dean flow).
pub(super) fn three_population_separation(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> ThreePopMetrics {
    let w = candidate.channel_width_m;
    let h = candidate.channel_height_m;

    // Per-arm flow: parallel topologies split Q across arms.
    let n_pv = candidate.topology.parallel_venturi_count();
    let n_ser = candidate.topology.serpentine_arm_count();
    let n_arms = n_pv.max(1).max(n_ser);
    let q_arm = candidate.flow_rate_m3_s / n_arms as f64;

    let mean_v = q_arm / (w * h);
    let shear_est = 6.0 * mean_v / h;
    let mu = blood.apparent_viscosity(shear_est);
    let rho = BLOOD_DENSITY_KG_M3;

    // Dean flow only for topologies that contain curved (serpentine) sections.
    let dh = 2.0 * w * h / (w + h);
    let bend_opt: Option<f64> = match candidate.topology {
        DesignTopology::SerpentineGrid
        | DesignTopology::VenturiSerpentine
        | DesignTopology::BifurcationSerpentine
        | DesignTopology::TrifurcationSerpentine
        | DesignTopology::AsymmetricBifurcationSerpentine => {
            Some(candidate.bend_radius_m.max(dh * 5.0))
        }
        DesignTopology::SpiralSerpentine { n_turns } => {
            let turn_len = TREATMENT_HEIGHT_MM * 1e-3 / n_turns as f64;
            let r = turn_len / (2.0 * std::f64::consts::PI);
            Some(r.max(dh * 5.0))
        }
        _ => None,
    };

    let eq = cfd_1d::three_population_equilibria(
        w,
        h,
        q_arm,
        rho,
        mu,
        candidate.feed_hematocrit,
        bend_opt,
    );

    ThreePopMetrics {
        sep_efficiency: eq.rbc_vs_center_sep,
        wbc_center_fraction: eq.wbc_center_fraction,
        cancer_center_fraction: eq.cancer_center_fraction,
        rbc_peripheral_fraction: eq.rbc_peripheral_fraction,
        wbc_eq_pos: eq.wbc_eq,
        cancer_eq_pos: eq.cancer_eq,
        rbc_eq_pos: eq.rbc_eq,
    }
}

// ── Leukapheresis separation ─────────────────────────────────────────────────

/// Intermediate results from the leukapheresis cell separation model.
///
/// Computed by [`leukapheresis_separation`] for
/// `ConstrictionExpansionArray`, `SpiralSerpentine`, and
/// `ParallelMicrochannelArray` topologies.
#[derive(Debug, Clone)]
pub(super) struct LeukapheresisMetrics {
    /// Fraction of WBCs focused into the center outlet channel.
    pub wbc_recovery: f64,
    /// Fraction of RBCs that pass into the center outlet channel (want low).
    pub rbc_pass_fraction: f64,
    /// WBC purity: `wbc_center / (wbc_center + rbc_center)`.
    pub wbc_purity: f64,
    /// Total extracorporeal volume [mL].
    pub total_ecv_ml: f64,
}

impl Default for LeukapheresisMetrics {
    fn default() -> Self {
        Self {
            wbc_recovery: 0.0,
            rbc_pass_fraction: 1.0,
            wbc_purity: 0.0,
            total_ecv_ml: 0.0,
        }
    }
}

/// Compute leukapheresis separation metrics for WBC and neonatal RBC using
/// the cell-free layer corrected lateral equilibrium model.
///
/// Uses [`enhanced_lateral_equilibrium`] at the candidate's `feed_hematocrit`
/// (CFL model from Fedosov 2012) rather than the plain inertial model.
///
/// # Outlet partition model
/// Leukapheresis chips expose two outlet streams (collection vs waste), but the
/// absolute streamline split is topology-specific.  To avoid imposing a brittle
/// fixed `x̃` threshold across spiral and straight microchannels, we infer a
/// per-candidate split plane at the midpoint between the WBC and RBC equilibrium
/// positions and map each population to the collection stream with a smooth
/// logistic transition.
pub(super) fn leukapheresis_separation(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> LeukapheresisMetrics {
    use cfd_1d::physics::cell_separation::{enhanced_lateral_equilibrium, CellProperties};

    let wbc = CellProperties::white_blood_cell();
    let rbc = CellProperties::neonatal_rbc();

    let w = candidate.channel_width_m;
    let h = candidate.channel_height_m;

    let n_arms = candidate.topology.serpentine_arm_count().max(1);
    let q_arm = candidate.flow_rate_m3_s / n_arms as f64;

    let mean_v = q_arm / (w * h);
    let shear_est = 6.0 * mean_v / h;
    let mu = blood.apparent_viscosity(shear_est);
    let rho = BLOOD_DENSITY_KG_M3;
    let dh = 2.0 * w * h / (w + h);

    let bend_opt: Option<f64> = match candidate.topology {
        DesignTopology::SpiralSerpentine { n_turns } => {
            let turn_len = TREATMENT_HEIGHT_MM * 1e-3 / n_turns as f64;
            let r = turn_len / (2.0 * std::f64::consts::PI);
            Some(r.max(dh * 5.0))
        }
        _ => None,
    };

    const DISPERSED: f64 = 0.5;
    let partition_sigma = match candidate.topology {
        DesignTopology::SpiralSerpentine { .. } => 0.008,
        _ => 0.03,
    };

    let eq = |cell: &CellProperties| -> f64 {
        enhanced_lateral_equilibrium(
            cell,
            rho,
            mu,
            mean_v,
            w,
            h,
            bend_opt,
            candidate.feed_hematocrit,
        )
        .map_or(DISPERSED, |r| r.x_tilde_eq)
    };

    let wbc_x = eq(&wbc);
    let rbc_x = eq(&rbc);
    let sep_gap = rbc_x - wbc_x;

    let (wbc_center, rbc_center) = if sep_gap <= 1e-9 {
        (0.0, 1.0)
    } else {
        let x_split = 0.5 * (wbc_x + rbc_x);
        let center_fraction =
            |x: f64| -> f64 { 1.0 / (1.0 + ((x - x_split) / partition_sigma).exp()) };
        (center_fraction(wbc_x), center_fraction(rbc_x))
    };

    let center_capture = wbc_center + rbc_center;
    let (rbc_pass_fraction, wbc_purity) = if center_capture > 1e-9 {
        (rbc_center, (wbc_center / center_capture).clamp(0.0, 1.0))
    } else {
        (1.0, 0.0)
    };

    let total_ecv_ml = match candidate.topology {
        DesignTopology::ParallelMicrochannelArray { n_channels } => {
            let ch_len = TREATMENT_HEIGHT_MM * 1e-3;
            n_channels as f64 * ch_len * w * h * 1e6
        }
        DesignTopology::ConstrictionExpansionArray { n_cycles } => {
            let wide_w = w;
            let narrow_w = w * 0.40;
            let seg_len = TREATMENT_HEIGHT_MM * 1e-3 / (n_cycles as f64 * 2.0);
            let narrow_len = seg_len * 0.5;
            n_cycles as f64 * (seg_len * wide_w * h + narrow_len * narrow_w * h) * 1e6
        }
        DesignTopology::SpiralSerpentine { n_turns } => {
            let turn_len = TREATMENT_HEIGHT_MM * 1e-3 / n_turns as f64;
            n_turns as f64 * turn_len * w * h * 1e6
        }
        _ => 0.0,
    };

    LeukapheresisMetrics {
        wbc_recovery: wbc_center,
        rbc_pass_fraction,
        wbc_purity,
        total_ecv_ml,
    }
}
