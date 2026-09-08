//! Core Zweifach-Fung routing probability functions and cell constants.
//!
//! Contains the fundamental physics of flow-weighted cell routing at
//! microfluidic junction bifurcations and trifurcations, including:
//!
//! - Base routing probability (`p_center`, `p_treat_bifurcation`, `p_arm_general`)
//! - Confinement-ratio amplification (`beta_kappa_adjusted`)
//! - Fåhræus margination correction (`fahrae_beta_correction`)
//! - High-velocity inversion (`pmc5114676_velocity_inversion`)
//!
//! # References
//! - Fung, Y. C. (1969). Biorheology of soft tissues. *Biorheology*, 6, 409–419.
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Yang et al. (2017). PMC5114676.
//! - Pries et al. (1989). *Circ. Res.* 64, 1198–1207.

// ── Stiffness exponents ─────────────────────────────────────────────────────
//
// Empirical, per Fung 1969 + modern microfluidic literature.
//
// Cancer exponent updated to 1.85 based on MCF-7 CTC routing data: large
// (17.5 µm), stiff (DI = 0.15) CTCs exhibit enhanced Zweifach-Fung routing
// in millifluidic bifurcations, consistent with β ∈ [1.8, 2.2] (Hou et al.
// 2012, *Lab Chip* 12, 1952; Karabacak et al. 2014, *Nat. Protoc.* 9, 694).
pub(crate) const SE_CANCER: f64 = 1.85; // MCF-7 breast cancer cells (stiff, ~17.5 µm diameter)
pub(crate) const SE_WBC: f64 = 1.40; // WBCs (semi-rigid, ~10 µm diameter)
pub(crate) const SE_RBC: f64 = 1.00; // RBCs (deformable — distributes by flow fraction)

// Cell diameters for confinement-ratio (κ = a/Dh) corrections.
pub(crate) const D_CANCER_M: f64 = 17.5e-6; // MCF-7 breast cancer cell diameter [m]
pub(crate) const D_WBC_M: f64 = 10.0e-6; // WBC diameter [m]
pub(crate) const D_RBC_M: f64 = 7.0e-6; // RBC diameter [m]

// Reference confinement ratio where inertial effects onset (Di Carlo 2009).
pub(crate) const KAPPA_REF: f64 = 0.07;

/// Fåhræus margination size-enhancement coefficient.
///
/// At physiological hematocrit the packed erythrocyte core displaces cells
/// larger than RBCs further from channel walls.  At branch points this
/// pre-positioning amplifies Zweifach-Fung routing toward high-flow arms
/// for cells whose diameter exceeds the RBC reference (Fåhræus 1929;
/// Pries et al. 1989, *Circ. Res.* 64, 1198–1207).
///
/// Coefficient 0.12 fitted to millifluidic CTC isolation efficiencies
/// at 40% hematocrit (Hou et al. 2013, *Sci. Rep.* 3, 1259).
pub(crate) const FAHRAE_SIZE_ALPHA: f64 = 0.12;

/// Probability that a cell enters the center (highest-flow) arm at an asymmetric
/// trifurcation junction.
///
/// Uses the 3-arm extension of the Zweifach-Fung bifurcation law (1969):
///
/// ```text
/// P_center = r_c^β / (r_c^β + 2 · r_p^β)
/// ```
///
/// where `r_c = Q_center / Q_total` (center-arm flow fraction), and
/// `r_p = (1 − r_c) / 2` (each peripheral-arm flow fraction).  The denominator
/// sums over all three arms individually rather than collapsing both peripheral
/// arms into a single "(1 − r)" term.  This correctly gives `P = 1/3` for all
/// cell types when the trifurcation is perfectly symmetric (r_c = r_p = 1/3),
/// and biases stiff cells (higher β) toward the center arm when r_c > 1/3.
///
/// # Arguments
/// * `q_center_frac` — fraction of total inlet flow carried by the center arm (0–1).
/// * `stiffness_exp` — cell stiffness exponent β.  Larger values → stronger bias
///   toward the higher-flow arm.
pub(crate) fn p_center(q_center_frac: f64, stiffness_exp: f64) -> f64 {
    let r_c = q_center_frac.clamp(1e-9, 1.0 - 1e-9);
    let r_p = (1.0 - r_c) * 0.5; // each peripheral arm's flow fraction
    let r_beta = r_c.powf(stiffness_exp);
    let p_beta = r_p.powf(stiffness_exp);
    r_beta / (r_beta + 2.0 * p_beta)
}

/// Probability that a cell enters the designated treatment arm at an
/// asymmetric bifurcation.
///
/// Uses the two-arm extension of the same flow-weighted routing law:
///
/// ```text
/// P_treat = r_t^β / (r_t^β + r_b^β)
/// ```
///
/// where `r_t = Q_treat / Q_total` and `r_b = 1 − r_t`.
pub(crate) fn p_treat_bifurcation(q_treat_frac: f64, stiffness_exp: f64) -> f64 {
    let r_t = q_treat_frac.clamp(1e-9, 1.0 - 1e-9);
    let r_b = 1.0 - r_t;
    let t_beta = r_t.powf(stiffness_exp);
    let b_beta = r_b.powf(stiffness_exp);
    t_beta / (t_beta + b_beta)
}

/// Additive β correction for Fåhræus margination at branch points.
///
/// Cells larger than the RBC reference experience enhanced displacement
/// from channel walls by the packed erythrocyte core, amplifying their
/// bifurcation routing bias beyond the confinement-ratio (κ) amplification.
///
/// Returns an additive β increment (≥ 0) proportional to the excess cell
/// diameter and the base stiffness excess.  For RBCs the correction is zero
/// (both excess β and excess diameter are zero).
#[inline]
pub(crate) fn fahrae_beta_correction(beta_base: f64, cell_diameter_m: f64) -> f64 {
    let excess = (beta_base - 1.0).max(0.0);
    let size_ratio = (cell_diameter_m / D_RBC_M - 1.0).max(0.0);
    excess * FAHRAE_SIZE_ALPHA * size_ratio
}

/// PMC5114676 Velocity-Dependent Zweifach-Fung Inversion Coefficient.
///
/// According to Yang et al. (2017), the probability of RBCs entering the higher-flow
/// branch decreases as the total inflow velocity increases. Beyond a critical velocity
/// threshold `v_crit` (~0.05 m/s in millifluidic expansions), the Zweifach-Fung
/// effect inverts: RBCs are propelled disproportionally into the lower flow-rate
/// (peripheral) branches due to excessive inertial centering forces overcoming the
/// streamline partition.
///
/// This returns a velocity-dependent multiplier for the baseline stiffness exponent β.
/// - If `v_in < v_crit`, returns `1.0` (classic Zweifach-Fung holds).
/// - If `v_in >= v_crit`, returns `< 1.0` and can go slightly negative (inversion).
#[inline]
pub(crate) fn pmc5114676_velocity_inversion(beta_base: f64, v_in: f64, is_rbc: bool) -> f64 {
    // Velocity threshold for Zweifach-Fung inversion (Yang et al., 2017)
    let v_crit = 0.05_f64; // ~50 mm/s
    if !is_rbc || v_in < v_crit {
        return 1.0;
    }

    // Decay rate for the RBC inversion.
    // Higher velocities severely diminish the exponent, driving RBCs to bypass streams.
    // β_inverted = β_base * exp(-k * (v_in - v_crit))
    let decay_rate = 25.0_f64;
    beta_base * (-decay_rate * (v_in - v_crit)).exp()
}

/// Confinement-ratio-adjusted stiffness exponent β for Zweifach-Fung routing.
///
/// As the treatment arm narrows through cascade stages, the confinement ratio
/// κ = cell_diameter / Dh grows.  For stiff cells, a higher κ means the cell
/// sits proportionally closer to the wall, amplifying the cross-stream lift
/// force and strengthening the Zweifach-Fung bias toward the high-flow arm.
///
/// The amplification scales the *excess* β above 1.0 (deformable baseline):
///   β_eff = 1 + (β_base − 1) × (1 + κ / κ_ref)
///
/// Under PMC5114676, high inflow velocities cause RBCs (which lack excess β)
/// to experience a diminishing/inverting base β via the velocity multiplier.
pub(crate) fn beta_kappa_adjusted(beta_base: f64, kappa: f64, v_in: f64, is_rbc: bool) -> f64 {
    let excess = (beta_base - 1.0).max(0.0);
    let amplification = 1.0 + kappa.clamp(0.0, KAPPA_REF * 2.0) / KAPPA_REF;

    // Apply PMC5114676 velocity inversion to the base exponent (primarily affects RBCs)
    let vel_multiplier = pmc5114676_velocity_inversion(beta_base, v_in, is_rbc);
    let base_eff = if is_rbc { 1.0 * vel_multiplier } else { 1.0 };

    #[allow(clippy::manual_clamp)]
    (base_eff + excess * amplification).min(3.0).max(0.1) // Lower bound ensures numeric stability
}

/// Routing probability into a single arm of an N-arm junction.
///
/// Generalises the Zweifach-Fung law to arbitrary arm count and asymmetric
/// flow fractions.  For N = 2 this reduces to `p_treat_bifurcation`; for
/// N = 3 with equal peripherals it reduces to `p_center`.
///
/// ```text
/// P_arm_i = q_i^β / Σ_j q_j^β
/// ```
pub(crate) fn p_arm_general(arm_q_fracs: &[f64], target_arm: usize, beta: f64) -> f64 {
    let sum_beta: f64 = arm_q_fracs
        .iter()
        .map(|&q| q.max(1e-9).powf(beta))
        .sum::<f64>()
        .max(1e-30);
    arm_q_fracs[target_arm].max(1e-9).powf(beta) / sum_beta
}

#[cfg(test)]
mod proptests {
    use super::*;

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_cascade_routing_monotonicity(
            q_center in 0.334..0.99_f64,
            beta1 in 1.01..2.0_f64,
            beta2 in 2.01..3.0_f64
        ) {
            let p1 = p_center(q_center, beta1);
            let p2 = p_center(q_center, beta2);

            // Theorem 1: Cascade Separation Monotonicity
            // For any asymmetric flow partition (q_center > 1/3), separating a stiffer cell
            // (higher beta) strictly increases the collection probability in the high-flow arm.
            prop_assert!(p2 > p1, "Failed monotonic stiffness amplification: {} <= {} for beta {} vs {}", p2, p1, beta2, beta1);

            // Theorem 2: Zweifach-Fung Amplification Law
            // For any particle with beta > 1.0, the fraction of particles entering the high-flow
            // branch strictly exceeds the fraction of bulk fluid entering that branch.
            prop_assert!(p1 > q_center, "Zweifach-Fung effect failed: P_center ({}) <= Q_center ({}) for beta {}", p1, q_center, beta1);
        }
    }
}
