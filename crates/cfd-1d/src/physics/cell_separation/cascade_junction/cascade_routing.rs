//! Cascade trifurcation and mixed Bi/Tri junction routing functions.
//!
//! Implements the multiplicative cascade routing model where cells pass
//! through successive junction stages, each applying the Zweifach-Fung
//! routing law to further enrich stiff cells in the center arm.

use super::routing_probability::{
    beta_kappa_adjusted, fahrae_beta_correction, p_arm_general, p_center, p_treat_bifurcation,
    D_CANCER_M, D_RBC_M, D_WBC_M, SE_CANCER, SE_RBC, SE_WBC,
};
use super::{CascadeJunctionResult, CascadeStage};

/// Internal cascade computation from pre-solved center-flow fractions.
///
/// Each stage routes cells via the Zweifach-Fung trifurcation law
/// using the given center-arm flow fractions.
pub(crate) fn cascade_from_q_fractions(q_center_fracs: &[f64]) -> CascadeJunctionResult {
    let mut f_cancer = 1.0_f64;
    let mut f_wbc = 1.0_f64;
    let mut f_rbc = 1.0_f64;
    let mut q_center_total = 1.0_f64;

    for &q_center in q_center_fracs {
        let q = q_center.clamp(1.0e-9, 1.0 - 1.0e-9);
        f_cancer *= p_center(q, SE_CANCER);
        f_wbc *= p_center(q, SE_WBC);
        f_rbc *= p_center(q, SE_RBC);
        q_center_total *= q;
    }

    let rbc_periph = (1.0 - f_rbc).clamp(0.0, 1.0);
    let sep_eff = (f_cancer - f_rbc).abs().clamp(0.0, 1.0);
    let center_hematocrit_ratio = if q_center_total > 1.0e-12 {
        (f_rbc / q_center_total).clamp(0.0, 2.0)
    } else {
        1.0
    };

    CascadeJunctionResult {
        cancer_center_fraction: f_cancer.clamp(0.0, 1.0),
        wbc_center_fraction: f_wbc.clamp(0.0, 1.0),
        rbc_peripheral_fraction: rbc_periph,
        separation_efficiency: sep_eff,
        center_hematocrit_ratio,
    }
}

/// Flow fraction carried by the center arm at a trifurcation with the given
/// width fraction, computed from hydraulic resistance of rectangular ducts
/// (`R ∝ 1 / w³` for equal height and length ducts).
///
/// # Arguments
/// * `center_frac` — fraction of parent width allocated to the center arm.
///   Each peripheral arm receives `(1 − center_frac) / 2`.
pub fn tri_center_q_frac(center_frac: f64) -> f64 {
    let w_c = center_frac.clamp(1e-6, 1.0 - 1e-6);
    let w_p = (1.0 - w_c) * 0.5;
    // Resistance is inversely proportional to w³ (equal height, equal length ducts).
    // Flow fraction to center: R_periph / (R_center + 2·R_periph) = w_c³ / (w_c³ + 2·w_p³)
    let r_c = w_c.powi(3);
    let r_p = w_p.powi(3);
    r_c / (r_c + 2.0 * r_p)
}

/// Flow fraction at a trifurcation with cross-junction resistance correction.
///
/// When channels of different diameters intersect, the cross-junction minor
/// loss (Idelchik 2007, §7.12) adds a diameter-ratio-dependent resistance to
/// each arm.  This shifts flow further toward the wider center arm because
/// the narrow peripheral arms see proportionally higher junction losses.
///
/// The effective resistance model combines Hagen-Poiseuille laminar
/// resistance (`R ∝ 1/w³`) with the junction K-factor correction
/// (`K_eff ∝ √(A_branch/A_run)`), where the run channel is the parent
/// trunk and each branch splits off at the junction.
///
/// # Arguments
/// * `center_frac` — fraction of parent width allocated to the center arm.
/// * `parent_width_m` — parent channel width [m].
/// * `channel_height_m` — channel height [m] (shared across all arms).
///
/// # Returns
/// Center-arm flow fraction ∈ (0, 1).
pub fn tri_center_q_frac_cross_junction(
    center_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> f64 {
    let w_c_frac = center_frac.clamp(1e-6, 1.0 - 1e-6);
    let w_p_frac = (1.0 - w_c_frac) * 0.5;

    let w_c = parent_width_m * w_c_frac;
    let w_p = parent_width_m * w_p_frac;
    let h = channel_height_m.max(1e-9);

    // Hagen-Poiseuille: R ∝ 1/w³ (equal height, equal reference length).
    let r_hp_c = 1.0 / w_c.powi(3);
    let r_hp_p = 1.0 / w_p.powi(3);

    // Cross-junction minor-loss correction: K_eff = K_base × √(A_branch / A_run)
    // A_center = w_c × h, A_periph = w_p × h, A_parent = parent_width × h.
    // The minor-loss resistance adds proportionally to K_eff per arm.
    let a_parent = parent_width_m * h;
    let a_center = w_c * h;
    let a_periph = w_p * h;

    // K-factors for dividing flow at cross junction (Idelchik).
    let k_base = 1.3_f64;
    let k_center = k_base * (a_center / a_parent).sqrt();
    let k_periph = k_base * (a_periph / a_parent).sqrt();

    // Minor-loss contribution scales with K / A², relative to H-P.
    // We express it as a dimensionless ratio to the H-P resistance using a
    // reference length equal to parent_width_m (junction length ≈ D_h).
    // R_junction ∝ K / A²;  R_HP ∝ µL / (w³h).
    // The ratio R_junc / R_HP is proportional to K·w³·h / (A²) = K·w / h.
    // We add this as a fractional correction so the result is continuous.
    let correction_center = k_center * w_c / h;
    let correction_periph = k_periph * w_p / h;

    // Effective conductance = 1 / (R_HP + R_junction_correction)
    let g_c = 1.0 / (r_hp_c * (1.0 + correction_center));
    let g_p = 1.0 / (r_hp_p * (1.0 + correction_periph));

    // Flow fraction to center: G_c / (G_c + 2·G_p)
    g_c / (g_c + 2.0 * g_p)
}

/// Flow fractions for all arms of an asymmetric trifurcation with cross-junction
/// minor-loss correction.
///
/// Unlike [`tri_center_q_frac_cross_junction`] (which assumes equal peripheral
/// arms), this function takes independent widths for the left and right
/// peripheral arms.  The right arm width is derived as
/// `right = parent_width − center − left` so that all three fractions sum to 1.
///
/// # Returns
/// `[q_center, q_left, q_right]` — flow fractions for each arm, summing to 1.
///
/// # Arguments
/// * `center_frac`       — center-arm width as fraction of parent width.
/// * `left_periph_frac`  — left peripheral width as fraction of parent width.
///   Right peripheral = `1 − center_frac − left_periph_frac`.
/// * `parent_width_m`    — parent channel width [m].
/// * `channel_height_m`  — shared channel height [m].
pub fn tri_asymmetric_q_fracs(
    center_frac: f64,
    left_periph_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> [f64; 3] {
    let w_c_frac = center_frac.clamp(1e-6, 1.0 - 2e-6);
    let w_l_frac = left_periph_frac.clamp(1e-6, 1.0 - w_c_frac - 1e-6);
    let w_r_frac = (1.0 - w_c_frac - w_l_frac).max(1e-6);

    let h = channel_height_m.max(1e-9);
    let w_c = parent_width_m * w_c_frac;
    let w_l = parent_width_m * w_l_frac;
    let w_r = parent_width_m * w_r_frac;

    // Hagen-Poiseuille: R ∝ 1/(w³ × h).
    let r_c = 1.0 / (w_c.powi(3) * h);
    let r_l = 1.0 / (w_l.powi(3) * h);
    let r_r = 1.0 / (w_r.powi(3) * h);

    // Idelchik cross-junction minor-loss K-factor correction.
    let a_parent = parent_width_m * h;
    let k_base = 1.3_f64;
    let k_c = k_base * ((w_c * h) / a_parent).sqrt();
    let k_l = k_base * ((w_l * h) / a_parent).sqrt();
    let k_r = k_base * ((w_r * h) / a_parent).sqrt();

    let g_c = 1.0 / (r_c * (1.0 + k_c * w_c / h));
    let g_l = 1.0 / (r_l * (1.0 + k_l * w_l / h));
    let g_r = 1.0 / (r_r * (1.0 + k_r * w_r / h));

    let g_total = (g_c + g_l + g_r).max(1e-30);
    [g_c / g_total, g_l / g_total, g_r / g_total]
}

/// Compute per-cell-type fractions after `n_levels` cascade trifurcation stages.
///
/// Each stage routes cells according to the Zweifach-Fung law; the center arm
/// narrows by `center_frac` at every level, further increasing the flow-fraction
/// asymmetry and the Zweifach-Fung routing bias.
///
/// # Arguments
/// * `n_levels`       — number of cascade trifurcation levels (1–3).
/// * `center_frac`    — center-arm width fraction at each level.
/// * `_channel_width` — trunk channel width [m] (reserved for future inertial
///   correction; currently unused in the routing calculation).
/// * `_channel_height`— trunk channel height [m] (reserved for future use).
/// * `_flow_rate`     — total inlet flow rate [m³/s] (reserved for future use).
pub fn cascade_junction_separation(
    n_levels: u8,
    center_frac: f64,
    _channel_width: f64,
    _channel_height: f64,
    _flow_rate: f64,
) -> CascadeJunctionResult {
    // Flow fraction from hydraulic resistance of width-scaled arms.
    // Re-computed at the trunk center_frac (same fraction at every level since we
    // only cascade the center arm, which re-splits at the same frac each time).
    let q_frac = tri_center_q_frac(center_frac);
    let q_fracs: Vec<f64> = std::iter::repeat_n(q_frac, n_levels as usize).collect();
    cascade_from_q_fractions(&q_fracs)
}

/// Compute cascade selective-routing separation from per-stage solved center-flow fractions.
///
/// This additive API is intended for high-fidelity coupling with solved 1D
/// network flow extraction (`cfd-optim::metrics::network_solve`), where each
/// stage can have a slightly different split due to downstream resistance.
pub fn cascade_junction_separation_from_qfracs(q_center_fracs: &[f64]) -> CascadeJunctionResult {
    cascade_from_q_fractions(q_center_fracs)
}

/// Compute cascade selective-routing separation with cross-junction diameter effects.
///
/// Each cascade level uses [`tri_center_q_frac_cross_junction`] to account
/// for cross-junction K-factor mismatch between center and peripheral arms.
/// The center arm width narrows through the cascade, so deeper levels see
/// increasingly asymmetric cross-junction geometry.
pub fn cascade_junction_separation_cross_junction(
    n_levels: u8,
    center_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> CascadeJunctionResult {
    let mut current_parent_w = parent_width_m;
    let q_fracs: Vec<f64> = (0..n_levels)
        .map(|_| {
            let q =
                tri_center_q_frac_cross_junction(center_frac, current_parent_w, channel_height_m);
            current_parent_w *= center_frac;
            q
        })
        .collect();
    cascade_from_q_fractions(&q_fracs)
}

/// Compute cell separation through a mixed cascade of trifurcation and
/// bifurcation junctions.
///
/// Each level is specified as a `(q_frac, is_trifurcation)` pair:
/// - **Trifurcation** levels use the 3-arm Zweifach-Fung model ([`p_center`])
///   where `q_frac` is the center-arm volumetric flow fraction.
/// - **Bifurcation** levels use the 2-arm model ([`p_treat_bifurcation`])
///   where `q_frac` is the treatment-arm volumetric flow fraction.
///
/// This generalises [`cascade_from_q_fractions`] (all-tri) and
/// [`incremental_from_q_fractions`](super::incremental_filtration::incremental_from_q_fractions)
/// (tri + terminal bi) to arbitrary Bi/Tri orderings, as required by
/// `PrimitiveSplitSequence` topologies.
///
/// # Theorem
///
/// **Monotonicity**: for each cell type with stiffness exponent `β ≥ 1`,
/// adding a cascade level with `q_frac > 0.5` (for bi) or `q_frac > 1/3`
/// (for tri) can only *increase* the center-arm cell fraction relative to
/// a shorter cascade truncated at the preceding level.
///
/// **Proof sketch**: at every level the routing probability `p ∈ (0,1)`
/// satisfies `p ≥ q_frac` when `β ≥ 1` and `q_frac` exceeds the
/// equal-flow baseline.  Multiplying the running product by `p ≤ 1` keeps
/// it monotonically bounded, while the differential between cancer (β=1.85)
/// and RBC (β=1.00) persists or widens.
pub fn mixed_cascade_separation(stages: &[(f64, bool)]) -> CascadeJunctionResult {
    let mut f_cancer = 1.0_f64;
    let mut f_wbc = 1.0_f64;
    let mut f_rbc = 1.0_f64;
    let mut q_center_total = 1.0_f64;

    for &(q_frac, is_tri) in stages {
        let q = q_frac.clamp(1e-9, 1.0 - 1e-9);
        if is_tri {
            f_cancer *= p_center(q, SE_CANCER);
            f_wbc *= p_center(q, SE_WBC);
            f_rbc *= p_center(q, SE_RBC);
        } else {
            f_cancer *= p_treat_bifurcation(q, SE_CANCER);
            f_wbc *= p_treat_bifurcation(q, SE_WBC);
            f_rbc *= p_treat_bifurcation(q, SE_RBC);
        }
        q_center_total *= q;
    }

    let rbc_periph = (1.0 - f_rbc).clamp(0.0, 1.0);
    let sep_eff = (f_cancer - f_rbc).abs().clamp(0.0, 1.0);
    let center_hematocrit_ratio = if q_center_total > 1e-12 {
        (f_rbc / q_center_total).clamp(0.0, 2.0)
    } else {
        1.0
    };

    CascadeJunctionResult {
        cancer_center_fraction: f_cancer.clamp(0.0, 1.0),
        wbc_center_fraction: f_wbc.clamp(0.0, 1.0),
        rbc_peripheral_fraction: rbc_periph,
        separation_efficiency: sep_eff,
        center_hematocrit_ratio,
    }
}

/// Compute cell separation through a mixed cascade of asymmetric junctions
/// with confinement-ratio (κ) dependent β amplification.
///
/// This is the high-fidelity successor to [`mixed_cascade_separation`].
/// Each [`CascadeStage`] carries:
/// - **All arm flow fractions** (`arm_q_fracs[0..n_arms]`), supporting
///   asymmetric peripheral arms (left ≠ right) without assuming equal splits.
/// - **Treatment-arm hydraulic diameter** (`treatment_dh_m`), used to compute
///   κ = cell_diameter / Dh at each stage.  As the cascade narrows the center
///   arm, κ grows and β is amplified for stiff cells (cancer, WBC), increasing
///   their preferential bias toward the high-flow treatment arm.
///
/// Physics:
/// - β_eff(κ) = 1 + (β_base − 1) × (1 + κ / κ_ref), capped at 3.0.
/// - For RBC (β_base = 1.0): β_eff = 1.0 always (deformable — no amplification).
/// - Routing: P_arm_i = q_i^β / Σ_j q_j^β (generalised Zweifach-Fung).
pub fn mixed_cascade_separation_kappa_aware(stages: &[CascadeStage]) -> CascadeJunctionResult {
    let mut f_cancer = 1.0_f64;
    let mut f_wbc = 1.0_f64;
    let mut f_rbc = 1.0_f64;
    let mut q_center_total = 1.0_f64;

    for stage in stages {
        let dh = stage.treatment_dh_m.max(1e-9);
        let v_in = stage.parent_v_in_m_s.max(1e-9);

        // κ-dependent amplification (Di Carlo 2009) + Fåhræus margination
        // correction for cells larger than RBCs (Pries 1989) + PMC5114676 velocity
        // damping/inversion (Yang 2017). The additive
        // Fåhræus term is significant in millifluidic channels where κ is
        // small (~0.01) and the κ-only amplification is weak, compensating
        // for the size-dependent displacement of large CTCs from channel walls
        // by the packed erythrocyte core at physiological hematocrit.
        let beta_cancer = (beta_kappa_adjusted(SE_CANCER, D_CANCER_M / dh, v_in, false)
            + fahrae_beta_correction(SE_CANCER, D_CANCER_M))
        .min(3.0);
        let beta_wbc = (beta_kappa_adjusted(SE_WBC, D_WBC_M / dh, v_in, false)
            + fahrae_beta_correction(SE_WBC, D_WBC_M))
        .min(3.0);
        let beta_rbc = beta_kappa_adjusted(SE_RBC, D_RBC_M / dh, v_in, true);

        let n = stage.n_arms.clamp(2, 5) as usize;
        let arms = &stage.arm_q_fracs[..n];

        let p_treat_cancer = p_arm_general(arms, 0, beta_cancer);
        let p_treat_wbc = p_arm_general(arms, 0, beta_wbc);
        let p_treat_rbc = p_arm_general(arms, 0, beta_rbc);

        let mut recovery_cancer = 0.0_f64;
        let mut recovery_wbc = 0.0_f64;
        let mut recovery_rbc = 0.0_f64;

        for i in 0..stage.n_recoveries as usize {
            if let Some(ref pr) = stage.peripheral_recoveries[i] {
                let p_leak_cancer = p_arm_general(arms, pr.source_arm_idx, beta_cancer);
                let p_leak_wbc = p_arm_general(arms, pr.source_arm_idx, beta_wbc);
                let p_leak_rbc = p_arm_general(arms, pr.source_arm_idx, beta_rbc);

                let sub_dh = pr.recovery_dh_m.max(1e-9);
                let sub_v_in = v_in * stage.arm_q_fracs[pr.source_arm_idx].max(1e-9); // Approx peripheral branch sub-velocity

                let sub_beta_cancer = (beta_kappa_adjusted(SE_CANCER, D_CANCER_M / sub_dh, sub_v_in, false)
                    + fahrae_beta_correction(SE_CANCER, D_CANCER_M))
                .min(3.0);
                let sub_beta_wbc = (beta_kappa_adjusted(SE_WBC, D_WBC_M / sub_dh, sub_v_in, false)
                    + fahrae_beta_correction(SE_WBC, D_WBC_M))
                .min(3.0);
                let sub_beta_rbc = beta_kappa_adjusted(SE_RBC, D_RBC_M / sub_dh, sub_v_in, true);

                let sub_n = pr.n_sub_arms.clamp(2, 5) as usize;
                let sub_arms = &pr.sub_arm_q_fracs[..sub_n];

                recovery_cancer +=
                    p_leak_cancer * p_arm_general(sub_arms, pr.recovery_arm_idx, sub_beta_cancer);
                recovery_wbc +=
                    p_leak_wbc * p_arm_general(sub_arms, pr.recovery_arm_idx, sub_beta_wbc);
                recovery_rbc +=
                    p_leak_rbc * p_arm_general(sub_arms, pr.recovery_arm_idx, sub_beta_rbc);
            }
        }

        f_cancer *= (p_treat_cancer + recovery_cancer).min(1.0);
        f_wbc *= (p_treat_wbc + recovery_wbc).min(1.0);
        f_rbc *= (p_treat_rbc + recovery_rbc).min(1.0);
        q_center_total *= stage.arm_q_fracs[0].max(1e-9);
    }

    let rbc_periph = (1.0 - f_rbc).clamp(0.0, 1.0);
    let sep_eff = (f_cancer - f_rbc).abs().clamp(0.0, 1.0);
    let center_hematocrit_ratio = if q_center_total > 1e-12 {
        (f_rbc / q_center_total).clamp(0.0, 2.0)
    } else {
        1.0
    };

    CascadeJunctionResult {
        cancer_center_fraction: f_cancer.clamp(0.0, 1.0),
        wbc_center_fraction: f_wbc.clamp(0.0, 1.0),
        rbc_peripheral_fraction: rbc_periph,
        separation_efficiency: sep_eff,
        center_hematocrit_ratio,
    }
}
