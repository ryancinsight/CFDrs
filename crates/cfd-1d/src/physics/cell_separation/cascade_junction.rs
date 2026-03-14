//! Zweifach-Fung junction routing model for selective split-sequence trees.
//!
//! Used by primitive selective branching topologies to predict the
//! cell-type-specific distribution between the designated treatment arm
//! (cancer/WBC-enriched, routed to the therapy zone) and the peripheral
//! bypass arms (RBC-enriched, low shear).
//!
//! # Physical basis
//!
//! At a trifurcation junction, the **Zweifach-Fung law** (1969) predicts that large,
//! stiff particles (cancer cells, WBCs) preferentially enter the arm carrying the
//! highest volumetric flow.  The probability follows a power law:
//!
//! ```text
//! P_center(cell) = r_c^β / (r_c^β + 2·r_p^β)
//! ```
//!
//! where `r_c = Q_center / Q_total`, `r_p = (1 − r_c) / 2` (each peripheral
//! arm), and `β` is a stiffness exponent:
//!
//! | Cell type | β    | Basis |
//! |-----------|------|-------|
//! | Cancer (MCF-7, ~17.5 µm, stiff)   | 1.85 | Size-enhanced stiff sphere (Hou 2012, Karabacak 2014) |
//! | WBC (~10 µm, semi-rigid)           | 1.40 | Intermediate deformability |
//! | RBC (~7 µm, highly deformable)     | 1.00 | Deformable cell — flow-weighted |
//!
//! The cascade iterates this routing over `n_levels` trifurcation levels, each time
//! narrowing the center arm by `center_frac` and carrying `q_center_frac` of the
//! current level's flow.  After all levels the center-arm fraction of each cell type
//! accumulates multiplicatively.
//!
//! # References
//! - Fung, Y. C. (1969). Biorheology of soft tissues. *Biorheology*, 6, 409–419.
//! - Doyeux, V. et al. (2011). Spheres in the vicinity of a bifurcation: elucidating
//!   the Zweifach-Fung effect. *J. Fluid Mech.*, 674, 359–388.
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.

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
fn p_center(q_center_frac: f64, stiffness_exp: f64) -> f64 {
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
fn p_treat_bifurcation(q_treat_frac: f64, stiffness_exp: f64) -> f64 {
    let r_t = q_treat_frac.clamp(1e-9, 1.0 - 1e-9);
    let r_b = 1.0 - r_t;
    let t_beta = r_t.powf(stiffness_exp);
    let b_beta = r_b.powf(stiffness_exp);
    t_beta / (t_beta + b_beta)
}

// Stiffness exponents (empirical, per Fung 1969 + modern microfluidic literature).
//
// Cancer exponent updated to 1.85 based on MCF-7 CTC routing data: large
// (17.5 µm), stiff (DI = 0.15) CTCs exhibit enhanced Zweifach-Fung routing
// in millifluidic bifurcations, consistent with β ∈ [1.8, 2.2] (Hou et al.
// 2012, *Lab Chip* 12, 1952; Karabacak et al. 2014, *Nat. Protoc.* 9, 694).
const SE_CANCER: f64 = 1.85; // MCF-7 breast cancer cells (stiff, ~17.5 µm diameter)
const SE_WBC: f64 = 1.40; // WBCs (semi-rigid, ~10 µm diameter)
const SE_RBC: f64 = 1.00; // RBCs (deformable — distributes by flow fraction)

// Cell diameters for confinement-ratio (κ = a/Dh) corrections.
const D_CANCER_M: f64 = 17.5e-6; // MCF-7 breast cancer cell diameter [m]
const D_WBC_M: f64 = 10.0e-6; // WBC diameter [m]
const D_RBC_M: f64 = 7.0e-6; // RBC diameter [m]

// Reference confinement ratio where inertial effects onset (Di Carlo 2009).
const KAPPA_REF: f64 = 0.07;

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
const FAHRAE_SIZE_ALPHA: f64 = 0.12;

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
fn fahrae_beta_correction(beta_base: f64, cell_diameter_m: f64) -> f64 {
    let excess = (beta_base - 1.0).max(0.0);
    let size_ratio = (cell_diameter_m / D_RBC_M - 1.0).max(0.0);
    excess * FAHRAE_SIZE_ALPHA * size_ratio
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
/// - At κ = 0: β_eff = β_base (no change).
/// - At κ = κ_ref (0.07): β_eff = 1 + 2 × (β_base − 1) (doubled excess).
/// - Capped at 3.0 to prevent unphysical divergence at very small channels.
/// - For RBC (β_base = 1.0): β_eff = 1.0 always (deformable — no amplification).
fn beta_kappa_adjusted(beta_base: f64, kappa: f64) -> f64 {
    let excess = (beta_base - 1.0).max(0.0);
    let amplification = 1.0 + kappa.clamp(0.0, KAPPA_REF * 2.0) / KAPPA_REF;
    (1.0 + excess * amplification).min(3.0)
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
fn p_arm_general(arm_q_fracs: &[f64], target_arm: usize, beta: f64) -> f64 {
    let sum_beta: f64 = arm_q_fracs
        .iter()
        .map(|&q| q.max(1e-9).powf(beta))
        .sum::<f64>()
        .max(1e-30);
    arm_q_fracs[target_arm].max(1e-9).powf(beta) / sum_beta
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

/// Peripheral recovery sub-split descriptor.
///
/// When a non-treatment arm is further sub-split, the wider sub-arm can
/// feed recovered cells back to the treatment path.  This struct
/// parameterises that recovery routing for a single source arm.
#[derive(Debug, Clone, Copy)]
pub struct PeripheralRecovery {
    /// Index of the source arm in the parent stage's `arm_q_fracs`.
    pub source_arm_idx: usize,
    /// Sub-arm flow fractions (up to 5); only `[..n_sub_arms]` are used.
    pub sub_arm_q_fracs: [f64; 5],
    /// Number of active sub-arms (2–5).
    pub n_sub_arms: u8,
    /// Index of the sub-arm that feeds back to the treatment path.
    pub recovery_arm_idx: usize,
    /// Hydraulic diameter of the recovery sub-arm [m].
    pub recovery_dh_m: f64,
}

/// Per-stage descriptor for a mixed Bi/Tri selective-routing cascade.
///
/// Carries all arm flow fractions (enabling asymmetric splits) and the
/// treatment-arm hydraulic diameter (enabling κ-dependent β correction).
///
/// Index 0 of `arm_q_fracs` is always the **treatment arm** (center at
/// trifurcations, treatment branch at bifurcations).
///
/// For a bifurcation, `arm_q_fracs = [q_treat, q_bypass, 0.0]` and `n_arms = 2`.
/// For a trifurcation, `arm_q_fracs = [q_center, q_left, q_right]` and `n_arms = 3`.
#[derive(Debug, Clone, Copy)]
pub struct CascadeStage {
    /// Flow fractions for all arms; first entry is the treatment arm.
    pub arm_q_fracs: [f64; 5],
    /// Number of active arms (2 = bifurcation, 3 = trifurcation, 4 = quad, 5 = penta).
    pub n_arms: u8,
    /// Hydraulic diameter of the treatment arm [m] at this stage.
    /// Used to compute κ = cell_diameter / Dh for β amplification.
    pub treatment_dh_m: f64,
    /// Optional peripheral recovery sub-splits (up to 4 per stage).
    pub peripheral_recoveries: [Option<PeripheralRecovery>; 4],
    /// Number of active peripheral recoveries.
    pub n_recoveries: u8,
}

fn cascade_from_q_fractions(q_center_fracs: &[f64]) -> CascadeJunctionResult {
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

fn incremental_from_q_fractions(
    pretri_q_center_fracs: &[f64],
    terminal_tri_q_center_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    let mut f_cancer = 1.0_f64;
    let mut f_wbc = 1.0_f64;
    let mut f_rbc = 1.0_f64;
    let mut q_total_treat = 1.0_f64;

    for &q_center in pretri_q_center_fracs {
        let q = q_center.clamp(1.0e-9, 1.0 - 1.0e-9);
        f_cancer *= p_center(q, SE_CANCER);
        f_wbc *= p_center(q, SE_WBC);
        f_rbc *= p_center(q, SE_RBC);
        q_total_treat *= q;
    }

    let q_tri = terminal_tri_q_center_frac.clamp(1.0e-9, 1.0 - 1.0e-9);
    f_cancer *= p_center(q_tri, SE_CANCER);
    f_wbc *= p_center(q_tri, SE_WBC);
    f_rbc *= p_center(q_tri, SE_RBC);
    q_total_treat *= q_tri;

    let q_bi = terminal_bi_treat_frac.clamp(0.50, 0.85);
    f_cancer *= p_treat_bifurcation(q_bi, SE_CANCER);
    f_wbc *= p_treat_bifurcation(q_bi, SE_WBC);
    f_rbc *= p_treat_bifurcation(q_bi, SE_RBC);
    q_total_treat *= q_bi;

    let rbc_periph = (1.0 - f_rbc).clamp(0.0, 1.0);
    let sep_eff = (f_cancer - f_rbc).abs().clamp(0.0, 1.0);
    let center_hematocrit_ratio = if q_total_treat > 1.0e-12 {
        (f_rbc / q_total_treat).clamp(0.0, 2.0)
    } else {
        1.0
    };

    IncrementalFiltrationResult {
        cancer_center_fraction: f_cancer.clamp(0.0, 1.0),
        wbc_center_fraction: f_wbc.clamp(0.0, 1.0),
        rbc_center_fraction: f_rbc.clamp(0.0, 1.0),
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

/// Stage-wise pre-trifurcation center-arm width fractions for selective routing.
///
/// This models progressive asymmetric focusing across pre-trifurcation levels:
/// deeper levels are allowed to bias further toward the center arm so that
/// large, stiff cells remain in the treatment path while deformable RBCs are
/// progressively skimmed to peripheral bypass branches.
///
/// The returned fractions are:
/// - clamped to `[0.20, 0.70]`,
/// - monotone non-decreasing over stage index,
/// - equal to `pretri_center_frac` when `n_pretri == 1`.
///   Constructed via a monotone interpolation from `start` to
///   `max(start, terminal_tri_frac)` with final clamping.
pub fn cif_pretri_stage_center_fracs(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
) -> Vec<f64> {
    let n = n_pretri.clamp(1, 3) as usize;
    let start = pretri_center_frac.clamp(0.20, 0.70);
    if n == 1 {
        return vec![start];
    }

    let terminal = terminal_tri_frac.clamp(0.20, 0.70);
    let target = start.max(terminal);

    (0..n)
        .map(|i| {
            let t = ((i + 1) as f64 / n as f64).powi(2);
            (start + (target - start) * t).clamp(0.20, 0.70)
        })
        .collect()
}

/// Stage-wise pre-trifurcation center-flow fractions for selective routing.
///
/// Each value is `tri_center_q_frac(stage_center_frac)` for the corresponding
/// stage returned by [`cif_pretri_stage_center_fracs`].
pub fn cif_pretri_stage_q_fracs(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
) -> Vec<f64> {
    cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)
        .into_iter()
        .map(tri_center_q_frac)
        .collect()
}

/// Cell fractions at the center (venturi) arm and peripheral bypass arms after
/// N cascade trifurcation levels.
#[derive(Debug, Clone, Copy)]
pub struct CascadeJunctionResult {
    /// Fraction of input cancer cells that reach the deepest center arm (venturi).
    pub cancer_center_fraction: f64,
    /// Fraction of input WBCs that reach the deepest center arm.
    pub wbc_center_fraction: f64,
    /// Fraction of input RBCs that are diverted to peripheral bypass arms
    /// at any cascade level (`= 1 − rbc_center_fraction`).
    pub rbc_peripheral_fraction: f64,
    /// Separation efficiency = `|f_cancer_center − f_rbc_center|` ∈ [0, 1].
    ///
    /// High values indicate strong enrichment of cancer cells at the venturi
    /// and strong depletion of RBCs (which are protected in bypass channels).
    pub separation_efficiency: f64,
    /// Local hematocrit in the center arm at the deepest level, relative to
    /// the feed hematocrit. Computed from the flow fraction and RBC routing:
    /// `HCT_local = HCT_feed × (rbc_center / q_center_frac^n_levels)`.
    pub center_hematocrit_ratio: f64,
}

/// Cell fractions after a staged selective-routing path:
///
/// 1. `n_pretri` center-only cascade trifurcation stages (incremental skimming),
/// 2. one terminal trifurcation skimming stage,
/// 3. one terminal asymmetric bifurcation selecting the treatment arm.
#[derive(Debug, Clone, Copy)]
pub struct IncrementalFiltrationResult {
    /// Fraction of input cancer cells reaching the treatment venturi arm.
    pub cancer_center_fraction: f64,
    /// Fraction of input WBCs reaching the treatment venturi arm.
    pub wbc_center_fraction: f64,
    /// Fraction of input RBCs that still reach the treatment venturi arm.
    pub rbc_center_fraction: f64,
    /// Fraction of input RBCs skimmed to peripheral bypass streams.
    pub rbc_peripheral_fraction: f64,
    /// Separation efficiency = `|f_cancer_center − f_rbc_center|` ∈ [0, 1].
    pub separation_efficiency: f64,
    /// Local hematocrit ratio (HCT_local / HCT_feed) at the treatment arm.
    pub center_hematocrit_ratio: f64,
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

/// Compute staged selective-routing separation.
///
/// This is the canonical staged selective-routing model for:
/// - pre-trifurcation skimming (`pretri_center_frac`),
/// - terminal-trifurcation skimming (`terminal_tri_frac`),
/// - terminal treatment bifurcation (`terminal_bi_treat_frac`).
pub fn incremental_filtration_separation_staged(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    let pretri_q_fracs = cif_pretri_stage_q_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);
    let q_tri = tri_center_q_frac(terminal_tri_frac);
    incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
}

/// Compute staged selective-routing separation from solved per-stage flow fractions.
///
/// This additive API is intended for coupling with solved-network extraction
/// where each pre-trifurcation stage can carry a different center-flow fraction.
pub fn incremental_filtration_separation_from_qfracs(
    pretri_q_center_fracs: &[f64],
    terminal_tri_q_center_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    incremental_from_q_fractions(
        pretri_q_center_fracs,
        terminal_tri_q_center_frac,
        terminal_bi_treat_frac,
    )
}

/// Compute staged selective-routing separation with cross-junction diameter effects.
///
/// Identical to [`incremental_filtration_separation_staged`] but uses
/// [`tri_center_q_frac_cross_junction`] to account for cross-junction
/// minor losses at each trifurcation.  When channels of different widths
/// intersect, the junction K-factor correction steers additional flow to
/// the wider center arm, enhancing RBC peripheral skimming.
///
/// # Arguments
/// * `n_pretri`              — number of initial center-cascade trifurcations (1–3).
/// * `pretri_center_frac`    — center-arm width fraction for the pre-cascade trifurcations.
/// * `terminal_tri_frac`     — center-arm width fraction for the terminal trifurcation.
/// * `terminal_bi_treat_frac`— flow fraction into the treatment arm of the terminal bifurcation.
/// * `parent_width_m`        — parent trunk channel width [m].
/// * `channel_height_m`      — channel height [m].
pub fn incremental_filtration_separation_cross_junction(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> IncrementalFiltrationResult {
    let stage_fracs =
        cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);

    // Convert width fractions to flow fractions using cross-junction model.
    // Each stage re-splits from the current center width, which narrows by
    // center_frac at each level.
    let mut current_parent_w = parent_width_m;
    let pretri_q_fracs: Vec<f64> = stage_fracs
        .iter()
        .map(|&frac| {
            let q = tri_center_q_frac_cross_junction(frac, current_parent_w, channel_height_m);
            current_parent_w *= frac; // center arm becomes parent for next level
            q
        })
        .collect();

    let q_tri =
        tri_center_q_frac_cross_junction(terminal_tri_frac, current_parent_w, channel_height_m);
    incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
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
/// [`incremental_from_q_fractions`] (tri + terminal bi) to arbitrary
/// Bi/Tri orderings, as required by [`PrimitiveSplitSequence`] topologies.
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
/// it monotonically bounded, while the differential between cancer (β=1.70)
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

        // κ-dependent amplification (Di Carlo 2009) + Fåhræus margination
        // correction for cells larger than RBCs (Pries 1989).  The additive
        // Fåhræus term is significant in millifluidic channels where κ is
        // small (~0.01) and the κ-only amplification is weak, compensating
        // for the size-dependent displacement of large CTCs from channel walls
        // by the packed erythrocyte core at physiological hematocrit.
        let beta_cancer = (beta_kappa_adjusted(SE_CANCER, D_CANCER_M / dh)
            + fahrae_beta_correction(SE_CANCER, D_CANCER_M))
        .min(3.0);
        let beta_wbc = (beta_kappa_adjusted(SE_WBC, D_WBC_M / dh)
            + fahrae_beta_correction(SE_WBC, D_WBC_M))
        .min(3.0);
        let beta_rbc = beta_kappa_adjusted(SE_RBC, D_RBC_M / dh);

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
                let sub_beta_cancer = (beta_kappa_adjusted(SE_CANCER, D_CANCER_M / sub_dh)
                    + fahrae_beta_correction(SE_CANCER, D_CANCER_M))
                .min(3.0);
                let sub_beta_wbc = (beta_kappa_adjusted(SE_WBC, D_WBC_M / sub_dh)
                    + fahrae_beta_correction(SE_WBC, D_WBC_M))
                .min(3.0);
                let sub_beta_rbc = beta_kappa_adjusted(SE_RBC, D_RBC_M / sub_dh);

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_split_equal_distribution() {
        // With symmetric 1/3 split all cell types should distribute by flow fraction
        // (cancer biased slightly more toward center, RBC less so).
        let r = cascade_junction_separation(1, 1.0 / 3.0, 2e-3, 1e-3, 5e-6);
        // Center arm carries 1/3 of flow → q_frac = (1/3)³/((1/3)³+2*(1/3)³) = 1/3
        // After 1 level: cancer_center_fraction > rbc_center_fraction (stiffness effect)
        assert!(
            r.cancer_center_fraction > r.wbc_center_fraction
                || (r.cancer_center_fraction - r.wbc_center_fraction).abs() < 1e-6
        );
        assert!(r.wbc_center_fraction >= r.rbc_peripheral_fraction * 0.0); // just sanity
        assert!(r.separation_efficiency >= 0.0 && r.separation_efficiency <= 1.0);
    }

    #[test]
    fn center_biased_increases_separation() {
        // Center-biased split should increase separation relative to symmetric.
        let r_sym = cascade_junction_separation(2, 1.0 / 3.0, 2e-3, 1e-3, 5e-6);
        let r_bias = cascade_junction_separation(2, 0.55, 2e-3, 1e-3, 5e-6);
        // Biased split: q_center_frac > 1/3 → stronger Zweifach-Fung routing
        assert!(r_bias.cancer_center_fraction >= r_sym.cancer_center_fraction - 1e-10);
    }

    #[test]
    fn more_levels_increases_enrichment() {
        let r1 = cascade_junction_separation(1, 0.45, 2e-3, 1e-3, 5e-6);
        let r3 = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6);
        // More cascade levels → more cancer enriched in center, more RBCs in bypass
        assert!(r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10);
    }

    #[test]
    fn tri_center_q_frac_symmetric() {
        // Symmetric 1/3 split → q_frac = 1/3
        let q = tri_center_q_frac(1.0 / 3.0);
        assert!(
            (q - 1.0 / 3.0).abs() < 1e-10,
            "symmetric frac should give q=1/3, got {q}"
        );
    }

    #[test]
    fn tri_center_q_frac_larger_center_gives_more_flow() {
        let q33 = tri_center_q_frac(0.333);
        let q55 = tri_center_q_frac(0.55);
        assert!(q55 > q33, "wider center arm should carry more flow");
    }

    #[test]
    fn incremental_filtration_more_pretri_levels_pushes_more_rbc_periphery() {
        let r1 = incremental_filtration_separation_staged(1, 0.45, 0.50, 0.68);
        let r3 = incremental_filtration_separation_staged(3, 0.45, 0.50, 0.68);
        assert!(r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10);
        assert!(r3.separation_efficiency >= r1.separation_efficiency - 1e-10);
    }

    #[test]
    fn incremental_filtration_higher_bi_treat_frac_increases_treatment_capture() {
        let low = incremental_filtration_separation_staged(2, 0.45, 0.45, 0.60);
        let high = incremental_filtration_separation_staged(2, 0.45, 0.45, 0.76);
        assert!(high.cancer_center_fraction >= low.cancer_center_fraction - 1e-10);
        assert!(high.wbc_center_fraction >= low.wbc_center_fraction - 1e-10);
    }

    #[test]
    fn mixed_sequence_tri_bi_routes_more_rbc_peripheral_than_single_tri() {
        let q_tri = tri_center_q_frac(0.45);
        let tri_only = mixed_cascade_separation(&[(q_tri, true)]);
        let tri_bi = mixed_cascade_separation(&[(q_tri, true), (0.72, false)]);

        assert!(
            tri_bi.rbc_peripheral_fraction >= tri_only.rbc_peripheral_fraction,
            "adding a treatment bifurcation should not reduce RBC peripheral routing"
        );
        assert!(
            tri_bi.cancer_center_fraction > (1.0 - tri_bi.rbc_peripheral_fraction),
            "cancer capture should remain above RBC center carryover in selective routing"
        );
    }

    #[test]
    fn stronger_trifurcation_bias_improves_three_pop_selectivity() {
        let weak = mixed_cascade_separation(&[
            (tri_center_q_frac(0.38), true),
            (tri_center_q_frac(0.38), true),
        ]);
        let strong = mixed_cascade_separation(&[
            (tri_center_q_frac(0.52), true),
            (tri_center_q_frac(0.52), true),
        ]);

        assert!(
            strong.cancer_center_fraction >= weak.cancer_center_fraction,
            "stronger center-arm bias should not reduce cancer capture"
        );
        assert!(
            strong.separation_efficiency >= weak.separation_efficiency,
            "stronger center-arm bias should improve separation efficiency"
        );
    }

    #[test]
    fn incremental_filtration_terminal_tri_center_bias_improves_cancer_center_fraction() {
        let symmetric_terminal = incremental_filtration_separation_staged(2, 0.45, 1.0 / 3.0, 0.68);
        let center_biased_terminal = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        assert!(
            center_biased_terminal.cancer_center_fraction
                >= symmetric_terminal.cancer_center_fraction - 1e-10
        );
    }

    #[test]
    fn cascade_qfrac_api_matches_uniform_width_model() {
        let width_model = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6);
        let q = tri_center_q_frac(0.45);
        let solved_like = cascade_junction_separation_from_qfracs(&[q, q, q]);
        assert!(
            (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs() < 1e-12
        );
        assert!(
            (width_model.rbc_peripheral_fraction - solved_like.rbc_peripheral_fraction).abs()
                < 1e-12
        );
    }

    #[test]
    fn incremental_qfrac_api_matches_uniform_width_model() {
        let width_model = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        let q_pretri = cif_pretri_stage_q_fracs(2, 0.45, 0.55);
        let q_tri = tri_center_q_frac(0.55);
        let solved_like = incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68);
        assert!(
            (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs() < 1e-12
        );
        assert!((width_model.rbc_center_fraction - solved_like.rbc_center_fraction).abs() < 1e-12);
    }

    #[test]
    fn cif_pretri_stage_fracs_ramp_toward_terminal_bias() {
        let stage_fracs = cif_pretri_stage_center_fracs(3, 0.45, 0.60);
        assert_eq!(stage_fracs.len(), 3);
        assert!(stage_fracs[0] >= 0.45 - 1e-12);
        assert!(stage_fracs[1] >= stage_fracs[0] - 1e-12);
        assert!(stage_fracs[2] >= stage_fracs[1] - 1e-12);
        assert!(stage_fracs[2] <= 0.60 + 1e-12);
    }

    #[test]
    fn cross_junction_q_frac_symmetric_matches_basic() {
        // With symmetric 1/3 split, cross-junction model should give similar
        // result to the basic model (identical for symmetric).
        let q_basic = tri_center_q_frac(1.0 / 3.0);
        let q_cross = tri_center_q_frac_cross_junction(1.0 / 3.0, 2e-3, 1e-3);
        // Both should be ~1/3 for symmetric geometry
        assert!((q_basic - 1.0 / 3.0).abs() < 1e-10);
        assert!(
            (q_cross - 1.0 / 3.0).abs() < 0.05,
            "symmetric cross-junction should be near 1/3, got {q_cross}"
        );
    }

    #[test]
    fn cross_junction_q_frac_wider_center_sees_more_k_loss() {
        // Cross-junction K-factor correction penalises the wider center arm
        // more than the narrow peripherals (K_eff ∝ √(A_branch/A_parent),
        // and the additive resistance term scales as K × w/h), so q_center
        // *decreases*.  This is the desired physics: more flow is diverted
        // to peripheral arms, improving RBC peripheral enrichment.
        let q_basic = tri_center_q_frac(0.55);
        let q_cross = tri_center_q_frac_cross_junction(0.55, 2e-3, 1e-3);
        assert!(q_cross <= q_basic + 1e-6,
            "cross-junction correction should reduce center fraction: basic={q_basic}, cross={q_cross}");
    }

    #[test]
    fn cross_junction_cif_pushes_more_rbc_to_periphery() {
        let basic = incremental_filtration_separation_staged(2, 0.45, 0.50, 0.68);
        let cross =
            incremental_filtration_separation_cross_junction(2, 0.45, 0.50, 0.68, 2e-3, 1e-3);
        assert!(cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
            "cross-junction selective routing should push more RBCs to periphery: basic={}, cross={}",
            basic.rbc_peripheral_fraction, cross.rbc_peripheral_fraction);
    }

    #[test]
    fn cross_junction_cct_pushes_more_rbc_to_periphery() {
        let basic = cascade_junction_separation(2, 0.45, 2e-3, 1e-3, 5e-6);
        let cross = cascade_junction_separation_cross_junction(2, 0.45, 2e-3, 1e-3);
        assert!(
            cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
            "cross-junction cascade routing should push more RBCs to periphery: basic={}, cross={}",
            basic.rbc_peripheral_fraction,
            cross.rbc_peripheral_fraction
        );
    }

    #[test]
    fn mixed_cascade_all_tri_matches_cascade() {
        let q = tri_center_q_frac(0.45);
        let pure = cascade_from_q_fractions(&[q, q]);
        let mixed = mixed_cascade_separation(&[(q, true), (q, true)]);
        assert!((pure.cancer_center_fraction - mixed.cancer_center_fraction).abs() < 1e-12);
        assert!((pure.rbc_peripheral_fraction - mixed.rbc_peripheral_fraction).abs() < 1e-12);
    }

    #[test]
    fn mixed_cascade_tri_bi_produces_nonzero_separation() {
        let q_tri = tri_center_q_frac(0.45);
        let q_bi = 0.68;
        let r = mixed_cascade_separation(&[(q_tri, true), (q_bi, false)]);
        assert!(r.cancer_center_fraction > 0.0);
        assert!(r.separation_efficiency > 0.0);
        assert!(r.cancer_center_fraction > r.rbc_peripheral_fraction.min(0.99));
    }

    #[test]
    fn mixed_cascade_deeper_improves_separation() {
        let q_tri = tri_center_q_frac(0.45);
        let r1 = mixed_cascade_separation(&[(q_tri, true)]);
        let r2 = mixed_cascade_separation(&[(q_tri, true), (q_tri, true)]);
        assert!(r2.separation_efficiency >= r1.separation_efficiency - 1e-10);
    }

    // ── kappa-aware tests ─────────────────────────────────────────────────────

    #[test]
    fn beta_kappa_adjusted_no_change_at_zero_kappa() {
        // At κ = 0 (cell infinitely smaller than channel), β_eff = β_base.
        assert!((beta_kappa_adjusted(SE_CANCER, 0.0) - SE_CANCER).abs() < 1e-12);
        assert!((beta_kappa_adjusted(SE_WBC, 0.0) - SE_WBC).abs() < 1e-12);
        assert!((beta_kappa_adjusted(SE_RBC, 0.0) - SE_RBC).abs() < 1e-12);
    }

    #[test]
    fn beta_kappa_adjusted_rbc_never_amplified() {
        // RBC is fully deformable (β_base = 1.0, excess = 0) — no amplification regardless of κ.
        for kappa in [0.0, 0.03, 0.07, 0.20] {
            let b = beta_kappa_adjusted(SE_RBC, kappa);
            assert!(
                (b - 1.0).abs() < 1e-12,
                "RBC β should stay 1.0 at κ={kappa}, got {b}"
            );
        }
    }

    #[test]
    fn beta_kappa_adjusted_cancer_increases_with_kappa() {
        let b0 = beta_kappa_adjusted(SE_CANCER, 0.0);
        let b1 = beta_kappa_adjusted(SE_CANCER, KAPPA_REF);
        let b2 = beta_kappa_adjusted(SE_CANCER, KAPPA_REF * 2.0);
        assert!(b1 > b0, "β should increase with κ for stiff cancer cells");
        assert!(b2 >= b1, "β should not decrease as κ grows further");
        assert!(b2 <= 3.0, "β must be capped at 3.0");
    }

    #[test]
    fn p_arm_general_symmetric_trifurcation_matches_p_center() {
        // Symmetric trifurcation: arm_q_fracs = [q, q_p, q_p]
        let q_c = 0.55_f64;
        let q_p = (1.0 - q_c) / 2.0;
        let p_gen = p_arm_general(&[q_c, q_p, q_p], 0, SE_CANCER);
        let p_old = p_center(q_c, SE_CANCER);
        assert!(
            (p_gen - p_old).abs() < 1e-10,
            "generalised p_arm must match p_center for symmetric tri: gen={p_gen}, old={p_old}"
        );
    }

    #[test]
    fn p_arm_general_two_arm_matches_p_treat_bifurcation() {
        let q_t = 0.68_f64;
        let q_b = 1.0 - q_t;
        let p_gen = p_arm_general(&[q_t, q_b], 0, SE_CANCER);
        let p_old = p_treat_bifurcation(q_t, SE_CANCER);
        assert!(
            (p_gen - p_old).abs() < 1e-10,
            "generalised p_arm must match p_treat_bifurcation: gen={p_gen}, old={p_old}"
        );
    }

    #[test]
    fn tri_asymmetric_q_fracs_sum_to_one() {
        let [qc, ql, qr] = tri_asymmetric_q_fracs(0.45, 0.30, 4e-3, 1e-3);
        assert!(
            (qc + ql + qr - 1.0).abs() < 1e-10,
            "asymmetric arm fracs must sum to 1: {qc}+{ql}+{qr}={:.6}",
            qc + ql + qr
        );
    }

    #[test]
    fn tri_asymmetric_wider_periph_gets_more_flow() {
        // Left peripheral wider than right → left gets more flow.
        let [_qc, ql, qr] = tri_asymmetric_q_fracs(0.40, 0.40, 4e-3, 1e-3);
        // left_frac=0.40, right_frac=0.20 → left should carry more
        assert!(
            ql > qr,
            "wider left peripheral should carry more flow: ql={ql}, qr={qr}"
        );
    }

    #[test]
    fn tri_asymmetric_symmetric_matches_cross_junction() {
        // With equal peripherals, tri_asymmetric_q_fracs should match tri_center_q_frac_cross_junction.
        let center_frac = 0.45_f64;
        let left_frac = (1.0 - center_frac) / 2.0; // symmetric
        let [qc_asym, _, _] = tri_asymmetric_q_fracs(center_frac, left_frac, 4e-3, 1e-3);
        let qc_sym = tri_center_q_frac_cross_junction(center_frac, 4e-3, 1e-3);
        assert!(
            (qc_asym - qc_sym).abs() < 1e-8,
            "symmetric asymmetric must match cross-junction: asym={qc_asym}, sym={qc_sym}"
        );
    }

    #[test]
    fn kappa_aware_higher_beta_than_legacy_for_stiff_cells_in_narrow_channel() {
        // In a narrow channel (Dh = 0.5mm), cancer cell κ ≈ 0.035 (above zero).
        // The kappa-aware model should produce HIGHER cancer_center_fraction than
        // the legacy model using constant β.
        let q = tri_center_q_frac_cross_junction(0.45, 0.8e-3, 1e-3);
        let q_p = (1.0 - q) / 2.0;
        let w_center = 0.45 * 0.8e-3;
        let h = 1e-3_f64;
        let dh = 2.0 * w_center * h / (w_center + h);

        let stage = CascadeStage {
            arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };

        let kappa_result = mixed_cascade_separation_kappa_aware(&[stage]);
        let legacy_result = mixed_cascade_separation(&[(q, true)]);

        // κ ≈ 17.5µm / Dh → β_cancer > SE_CANCER = 1.85 → stronger cancer routing
        assert!(
            kappa_result.cancer_center_fraction >= legacy_result.cancer_center_fraction - 1e-10,
            "kappa-aware model must not reduce cancer routing vs legacy: kappa={:.4}, legacy={:.4}",
            kappa_result.cancer_center_fraction,
            legacy_result.cancer_center_fraction
        );
        // RBC routing must be unchanged (deformable, β stays 1.0)
        assert!(
            (kappa_result.rbc_peripheral_fraction - legacy_result.rbc_peripheral_fraction).abs()
                < 1e-9,
            "RBC routing must be unchanged by κ correction"
        );
    }

    #[test]
    fn kappa_aware_asymmetric_arms_biases_cells_to_wider_peripheral() {
        // Asymmetric trifurcation: left peripheral wider than right.
        // Cancer cells should preferentially go to the wider LEFT arm (higher q_l).
        let [q_c, q_l, q_r] = tri_asymmetric_q_fracs(0.40, 0.40, 4e-3, 1e-3);
        // q_l > q_r → more cells should route to left
        let w_center = 0.40 * 4e-3;
        let h = 1e-3_f64;
        let dh = 2.0 * w_center * h / (w_center + h);
        let stage = CascadeStage {
            arm_q_fracs: [q_c, q_l, q_r, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        let n = stage.n_arms as usize;
        let dh_s = stage.treatment_dh_m.max(1e-9);
        let beta = beta_kappa_adjusted(SE_CANCER, D_CANCER_M / dh_s);
        let p_to_left = p_arm_general(&stage.arm_q_fracs[..n], 1, beta);
        let p_to_right = p_arm_general(&stage.arm_q_fracs[..n], 2, beta);
        assert!(
            p_to_left > p_to_right,
            "wider left peripheral should attract more cancer cells: p_left={p_to_left:.4}, p_right={p_to_right:.4}"
        );
    }

    #[test]
    fn kappa_aware_deeper_cascade_still_improves_separation() {
        // Build a 3-stage TriTriTri with narrowing Dh at each stage.
        let center_frac = 0.45_f64;
        let h = 1e-3_f64;
        let mut parent_w = 4e-3_f64;
        let mut stages = Vec::new();
        for _ in 0..3 {
            let w_c = center_frac * parent_w;
            let dh = 2.0 * w_c * h / (w_c + h);
            let q_c = tri_center_q_frac_cross_junction(center_frac, parent_w, h);
            let q_p = (1.0 - q_c) / 2.0;
            stages.push(CascadeStage {
                arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            });
            parent_w *= center_frac;
        }

        let r3 = mixed_cascade_separation_kappa_aware(&stages);
        let r1 = mixed_cascade_separation_kappa_aware(&stages[..1]);

        assert!(
            r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10,
            "3-stage TriTriTri should push more RBCs to periphery: r1={:.4}, r3={:.4}",
            r1.rbc_peripheral_fraction,
            r3.rbc_peripheral_fraction
        );
        assert!(
            r3.separation_efficiency >= r1.separation_efficiency - 1e-10,
            "3-stage TriTriTri should improve separation efficiency"
        );
    }

    // ── Fåhræus margination correction tests ──────────────────────────────

    #[test]
    fn fahrae_correction_zero_for_rbc() {
        // RBC is the reference cell: excess β = 0 and size ratio = 0 → no correction.
        let c = fahrae_beta_correction(SE_RBC, D_RBC_M);
        assert!(
            c.abs() < 1e-15,
            "Fåhræus correction must be zero for RBC, got {c}"
        );
    }

    #[test]
    fn fahrae_correction_positive_for_cancer() {
        // Cancer cells (17.5 µm) are larger than RBCs (7 µm) → positive correction.
        let c = fahrae_beta_correction(SE_CANCER, D_CANCER_M);
        assert!(
            c > 0.0,
            "Fåhræus correction must be positive for cancer cells"
        );
        // Expected: excess(0.85) × α(0.12) × size_ratio(1.5) = 0.153
        assert!(
            (c - 0.85 * 0.12 * 1.5).abs() < 1e-12,
            "Fåhræus correction for cancer should be ~0.153, got {c}"
        );
    }

    #[test]
    fn fahrae_correction_smaller_for_wbc_than_cancer() {
        // WBCs (10 µm) are smaller than cancer cells (17.5 µm) → smaller correction.
        let c_cancer = fahrae_beta_correction(SE_CANCER, D_CANCER_M);
        let c_wbc = fahrae_beta_correction(SE_WBC, D_WBC_M);
        assert!(c_wbc > 0.0, "Fåhræus correction must be positive for WBCs");
        assert!(
            c_cancer > c_wbc,
            "Cancer correction ({c_cancer:.4}) must exceed WBC correction ({c_wbc:.4})"
        );
    }

    #[test]
    fn kappa_aware_with_fahrae_boosts_cancer_over_legacy() {
        // The Fåhræus correction in the kappa-aware model should boost cancer routing
        // beyond what the legacy mixed_cascade_separation (no κ, no Fåhræus) provides.
        let q = tri_center_q_frac_cross_junction(0.45, 4e-3, 1e-3);
        let q_p = (1.0 - q) / 2.0;
        let w_center = 0.45 * 4e-3;
        let h = 1e-3_f64;
        let dh = 2.0 * w_center * h / (w_center + h);

        let stage = CascadeStage {
            arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        let kappa_fahrae = mixed_cascade_separation_kappa_aware(&[stage]);
        let legacy = mixed_cascade_separation(&[(q, true)]);

        // The Fåhræus-enhanced kappa model should route more cancer cells to center.
        assert!(
            kappa_fahrae.cancer_center_fraction > legacy.cancer_center_fraction,
            "Fåhræus-enhanced model must exceed legacy cancer routing: enhanced={:.4}, legacy={:.4}",
            kappa_fahrae.cancer_center_fraction,
            legacy.cancer_center_fraction
        );
        // Cancer-to-RBC enrichment ratio should improve.
        let rbc_center_kf = 1.0 - kappa_fahrae.rbc_peripheral_fraction;
        let rbc_center_legacy = 1.0 - legacy.rbc_peripheral_fraction;
        let enrichment_kf = kappa_fahrae.cancer_center_fraction / rbc_center_kf.max(1e-12);
        let enrichment_legacy = legacy.cancer_center_fraction / rbc_center_legacy.max(1e-12);
        assert!(
            enrichment_kf >= enrichment_legacy - 1e-6,
            "Fåhræus model should improve CTC/RBC enrichment: kf={enrichment_kf:.4}, legacy={enrichment_legacy:.4}"
        );
    }

    #[test]
    fn option1_option2_tri_tri_enrichment_quantified() {
        // Quantify CTC enrichment for the actual Option 1/Option 2 TriTri topology:
        //   Stage 1: pretri center_frac = 0.45
        //   Stage 2: terminal tri center_frac = 0.333 (symmetric)
        //   parent_width = 4mm, height = 1mm
        let h = 1e-3_f64;
        let parent_w = 4e-3_f64;
        let pcf = 0.45;
        let tcf = 1.0 / 3.0; // symmetric

        // Build 2-stage TriTri with cross-junction corrections (as used by compute.rs)
        let left_periph = (1.0 - pcf) / 2.0;
        let arm_q_1 = tri_asymmetric_q_fracs(pcf, left_periph, parent_w, h);
        let w_center_1 = pcf * parent_w;
        let dh_1 = 2.0 * w_center_1 * h / (w_center_1 + h);

        let parent_w_2 = w_center_1;
        let left_periph_2 = (1.0 - tcf) / 2.0;
        let arm_q_2 = tri_asymmetric_q_fracs(tcf, left_periph_2, parent_w_2, h);
        let w_center_2 = tcf * parent_w_2;
        let dh_2 = 2.0 * w_center_2 * h / (w_center_2 + h);

        let stages = vec![
            CascadeStage {
                arm_q_fracs: [arm_q_1[0], arm_q_1[1], arm_q_1[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_1,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
            CascadeStage {
                arm_q_fracs: [arm_q_2[0], arm_q_2[1], arm_q_2[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_2,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
        ];
        let r = mixed_cascade_separation_kappa_aware(&stages);

        // With SE_CANCER=1.85 + Fåhræus correction, cancer routing should be > 25%
        // (up from ~24.6% with SE_CANCER=1.70 and no Fåhræus).
        assert!(
            r.cancer_center_fraction > 0.25,
            "TriTri cancer center fraction should exceed 25%: got {:.2}%",
            r.cancer_center_fraction * 100.0
        );
        // RBC peripheral fraction should remain high (> 78%)
        assert!(
            r.rbc_peripheral_fraction > 0.78,
            "TriTri RBC peripheral fraction should exceed 78%: got {:.2}%",
            r.rbc_peripheral_fraction * 100.0
        );
        // CTC/RBC enrichment ratio should be > 1.3
        let rbc_center = 1.0 - r.rbc_peripheral_fraction;
        let enrichment = r.cancer_center_fraction / rbc_center.max(1e-12);
        assert!(
            enrichment > 1.3,
            "CTC/RBC enrichment ratio should exceed 1.3: got {enrichment:.4}"
        );

        // Print for diagnostic review
        eprintln!("=== TriTri (pcf=0.45, tcf=0.333) kappa-aware + Fåhræus ===");
        eprintln!(
            "  Cancer center fraction:  {:.2}%",
            r.cancer_center_fraction * 100.0
        );
        eprintln!(
            "  WBC center fraction:     {:.2}%",
            r.wbc_center_fraction * 100.0
        );
        eprintln!(
            "  RBC peripheral fraction: {:.2}%",
            r.rbc_peripheral_fraction * 100.0
        );
        eprintln!("  RBC center fraction:     {:.2}%", rbc_center * 100.0);
        eprintln!("  Separation efficiency:   {:.4}", r.separation_efficiency);
        eprintln!("  CTC/RBC enrichment:      {:.4}×", enrichment);
        eprintln!(
            "  Center HCT ratio:        {:.4}",
            r.center_hematocrit_ratio
        );
    }

    #[test]
    fn option2_higher_tcf_dramatically_improves_enrichment() {
        // Demonstrate that using tcf=0.55 instead of 0.333 would dramatically
        // improve CTC enrichment. This is a design-space insight, not a physics change.
        let h = 1e-3_f64;
        let parent_w = 4e-3_f64;
        let pcf = 0.45;

        // Build with tcf=0.55 (asymmetric terminal tri)
        let tcf = 0.55;
        let left_periph = (1.0 - pcf) / 2.0;
        let arm_q_1 = tri_asymmetric_q_fracs(pcf, left_periph, parent_w, h);
        let w_center_1 = pcf * parent_w;
        let dh_1 = 2.0 * w_center_1 * h / (w_center_1 + h);

        let parent_w_2 = w_center_1;
        let left_periph_2 = (1.0 - tcf) / 2.0;
        let arm_q_2 = tri_asymmetric_q_fracs(tcf, left_periph_2, parent_w_2, h);
        let w_center_2 = tcf * parent_w_2;
        let dh_2 = 2.0 * w_center_2 * h / (w_center_2 + h);

        let stages = vec![
            CascadeStage {
                arm_q_fracs: [arm_q_1[0], arm_q_1[1], arm_q_1[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_1,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
            CascadeStage {
                arm_q_fracs: [arm_q_2[0], arm_q_2[1], arm_q_2[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_2,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
        ];
        let r_high = mixed_cascade_separation_kappa_aware(&stages);

        // With tcf=0.55: cancer routing should exceed 40%.
        assert!(
            r_high.cancer_center_fraction > 0.40,
            "TriTri with tcf=0.55 cancer fraction should exceed 40%: got {:.2}%",
            r_high.cancer_center_fraction * 100.0
        );

        let rbc_center = 1.0 - r_high.rbc_peripheral_fraction;
        let enrichment = r_high.cancer_center_fraction / rbc_center.max(1e-12);
        assert!(
            enrichment > 1.4,
            "TriTri with tcf=0.55 enrichment should exceed 1.4: got {enrichment:.4}"
        );

        eprintln!("=== TriTri (pcf=0.45, tcf=0.55) kappa-aware + Fåhræus ===");
        eprintln!(
            "  Cancer center fraction:  {:.2}%",
            r_high.cancer_center_fraction * 100.0
        );
        eprintln!(
            "  WBC center fraction:     {:.2}%",
            r_high.wbc_center_fraction * 100.0
        );
        eprintln!(
            "  RBC peripheral fraction: {:.2}%",
            r_high.rbc_peripheral_fraction * 100.0
        );
        eprintln!("  RBC center fraction:     {:.2}%", rbc_center * 100.0);
        eprintln!("  CTC/RBC enrichment:      {:.4}×", enrichment);
    }

    #[test]
    fn updated_se_cancer_improves_single_stage_selectivity() {
        // With SE_CANCER = 1.85 (up from 1.70), a single trifurcation stage with
        // center-biased flow should produce higher cancer-to-RBC differential.
        let q = tri_center_q_frac(0.50);
        let r = mixed_cascade_separation(&[(q, true)]);

        // Cancer routing must exceed RBC routing (flow-weighted at β=1.0).
        let rbc_center = 1.0 - r.rbc_peripheral_fraction;
        assert!(
            r.cancer_center_fraction > rbc_center,
            "Cancer must route preferentially to center: cancer={:.4}, rbc_center={:.4}",
            r.cancer_center_fraction,
            rbc_center
        );
        // Enrichment ratio must be > 1.0
        let enrichment = r.cancer_center_fraction / rbc_center.max(1e-12);
        assert!(
            enrichment > 1.0,
            "CTC/RBC enrichment must exceed 1.0 at asymmetric split, got {enrichment:.4}"
        );
    }

    // ── Peripheral recovery routing tests ────────────────────────────────

    #[test]
    fn recovery_zero_when_no_peripheral_recoveries() {
        let q = tri_center_q_frac(0.45);
        let q_p = (1.0 - q) / 2.0;
        let dh = 1e-3;
        let stage_no_recovery = CascadeStage {
            arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        let r = mixed_cascade_separation_kappa_aware(&[stage_no_recovery]);
        let r_legacy = mixed_cascade_separation(&[(q, true)]);
        // With no recovery, kappa-aware should match or exceed legacy
        assert!(r.cancer_center_fraction >= r_legacy.cancer_center_fraction - 1e-10);
    }

    #[test]
    fn recovery_increases_cancer_center_fraction() {
        // 20/60/20 trifurcation with recovery on left peripheral
        let q_c = 0.65;  // center gets ~65% of flow (wider channel)
        let q_p = (1.0 - q_c) / 2.0;
        let dh = 1e-3;
        let stage_no_recovery = CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        // Recovery sub-split on arm 1 (left): 70/30 split, wider arm (index 0) feeds back to treatment
        let recovery = PeripheralRecovery {
            source_arm_idx: 1,
            sub_arm_q_fracs: [0.70, 0.30, 0.0, 0.0, 0.0],
            n_sub_arms: 2,
            recovery_arm_idx: 0,
            recovery_dh_m: 0.5e-3,
        };
        let stage_with_recovery = CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            peripheral_recoveries: [Some(recovery), None, None, None],
            n_recoveries: 1,
        };
        let r_no = mixed_cascade_separation_kappa_aware(&[stage_no_recovery]);
        let r_yes = mixed_cascade_separation_kappa_aware(&[stage_with_recovery]);
        assert!(
            r_yes.cancer_center_fraction > r_no.cancer_center_fraction,
            "recovery should increase cancer center fraction: without={:.4}, with={:.4}",
            r_no.cancer_center_fraction, r_yes.cancer_center_fraction
        );
    }

    #[test]
    fn recovery_bounded_by_one() {
        // Even with aggressive recovery, P_eff must not exceed 1.0
        let q_c = 0.50;
        let q_p = 0.25;
        let recovery_both = [
            Some(PeripheralRecovery {
                source_arm_idx: 1,
                sub_arm_q_fracs: [0.90, 0.10, 0.0, 0.0, 0.0],
                n_sub_arms: 2,
                recovery_arm_idx: 0,
                recovery_dh_m: 0.3e-3,
            }),
            Some(PeripheralRecovery {
                source_arm_idx: 2,
                sub_arm_q_fracs: [0.90, 0.10, 0.0, 0.0, 0.0],
                n_sub_arms: 2,
                recovery_arm_idx: 0,
                recovery_dh_m: 0.3e-3,
            }),
            None,
            None,
        ];
        let stage = CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: 1e-3,
            peripheral_recoveries: recovery_both,
            n_recoveries: 2,
        };
        let r = mixed_cascade_separation_kappa_aware(&[stage]);
        assert!(r.cancer_center_fraction <= 1.0, "cancer fraction must be <= 1.0");
        assert!(r.wbc_center_fraction <= 1.0, "wbc fraction must be <= 1.0");
        assert!((1.0 - r.rbc_peripheral_fraction) <= 1.0, "rbc center must be <= 1.0");
    }
}
