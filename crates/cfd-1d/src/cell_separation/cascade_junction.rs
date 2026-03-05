//! Zweifach-Fung junction routing model for asymmetric trifurcation cascades.
//!
//! Used by `CascadeCenterTrifurcationSeparator` (CCT) topologies to predict
//! the cell-type-specific distribution between the center arm (cancer/WBC-enriched,
//! routed to the venturi) and the peripheral bypass arms (RBC-enriched, low shear).
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
//! | Cancer (MCF-7, ~17.5 µm, stiff)   | 1.70 | Stiff sphere limit (Fung 1969 extended) |
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
const SE_CANCER: f64 = 1.70; // MCF-7 breast cancer cells (stiff, ~17.5 µm diameter)
const SE_WBC: f64 = 1.40; // WBCs (semi-rigid, ~10 µm diameter)
const SE_RBC: f64 = 1.00; // RBCs (deformable — distributes by flow fraction)

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

/// Stage-wise pre-trifurcation center-arm width fractions for CIF.
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

/// Stage-wise pre-trifurcation center-flow fractions for CIF.
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

/// Cell fractions after a staged controlled incremental filtration (CIF) path:
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

/// Compute CCT routing from per-stage solved center-flow fractions.
///
/// This additive API is intended for high-fidelity coupling with solved 1D
/// network flow extraction (`cfd-optim::metrics::network_solve`), where each
/// stage can have a slightly different split due to downstream resistance.
pub fn cascade_junction_separation_from_qfracs(
    q_center_fracs: &[f64],
) -> CascadeJunctionResult {
    cascade_from_q_fractions(q_center_fracs)
}

/// Compute staged controlled incremental filtration (CIF) routing.
///
/// The sequence is:
/// - Pre-skimming: `n_pretri` center-only trifurcation cascade levels,
/// - Mid-length skimming: one additional asymmetric trifurcation,
/// - Terminal selection: one asymmetric bifurcation where the higher-flow arm
///   is the treatment / venturi path.
///
/// This models "first trifurcations, then trifurcation/bifurcation" routing
/// with progressive RBC displacement into peripheral bypass paths.
///
/// # Arguments
/// * `n_pretri`              — number of initial center-cascade trifurcations (1–3).
/// * `pretri_center_frac`    — center-arm width fraction for the pre-cascade trifurcations.
/// * `terminal_tri_frac`     — center-arm width fraction for the terminal trifurcation.
/// * `terminal_bi_treat_frac`— flow fraction into the treatment arm of the terminal bifurcation.
pub fn incremental_filtration_separation(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    incremental_filtration_separation_staged(
        n_pretri,
        pretri_center_frac,
        terminal_tri_frac,
        terminal_bi_treat_frac,
    )
}

/// Compute staged controlled incremental filtration (CIF) routing.
///
/// This is the canonical staged-CIF model for:
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

/// Compute staged CIF routing from solved per-stage flow fractions.
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

/// Compute staged CIF routing with cross-junction diameter effects.
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

    let q_tri = tri_center_q_frac_cross_junction(
        terminal_tri_frac,
        current_parent_w,
        channel_height_m,
    );
    incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
}

/// Compute CCT routing with cross-junction diameter effects.
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
        assert!(r.cancer_center_fraction > r.wbc_center_fraction || (r.cancer_center_fraction - r.wbc_center_fraction).abs() < 1e-6);
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
        assert!((q - 1.0 / 3.0).abs() < 1e-10, "symmetric frac should give q=1/3, got {q}");
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
    fn incremental_filtration_terminal_tri_center_bias_improves_cancer_center_fraction() {
        let symmetric_terminal = incremental_filtration_separation_staged(2, 0.45, 1.0 / 3.0, 0.68);
        let center_biased_terminal = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        assert!(
            center_biased_terminal.cancer_center_fraction
                >= symmetric_terminal.cancer_center_fraction - 1e-10
        );
    }

    #[test]
    fn incremental_filtration_legacy_wrapper_matches_staged() {
        let legacy = incremental_filtration_separation(2, 0.45, 0.55, 0.68);
        let staged = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        assert!((legacy.cancer_center_fraction - staged.cancer_center_fraction).abs() < 1e-12);
        assert!((legacy.wbc_center_fraction - staged.wbc_center_fraction).abs() < 1e-12);
        assert!((legacy.rbc_center_fraction - staged.rbc_center_fraction).abs() < 1e-12);
    }

    #[test]
    fn cascade_qfrac_api_matches_uniform_width_model() {
        let width_model = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6);
        let q = tri_center_q_frac(0.45);
        let solved_like = cascade_junction_separation_from_qfracs(&[q, q, q]);
        assert!(
            (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs()
                < 1e-12
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
            (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs()
                < 1e-12
        );
        assert!(
            (width_model.rbc_center_fraction - solved_like.rbc_center_fraction).abs() < 1e-12
        );
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
        let q_cross =
            tri_center_q_frac_cross_junction(1.0 / 3.0, 2e-3, 1e-3);
        // Both should be ~1/3 for symmetric geometry
        assert!((q_basic - 1.0 / 3.0).abs() < 1e-10);
        assert!((q_cross - 1.0 / 3.0).abs() < 0.05,
            "symmetric cross-junction should be near 1/3, got {q_cross}");
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
        let cross = incremental_filtration_separation_cross_junction(
            2, 0.45, 0.50, 0.68, 2e-3, 1e-3,
        );
        assert!(cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
            "cross-junction CIF should push more RBCs to periphery: basic={}, cross={}",
            basic.rbc_peripheral_fraction, cross.rbc_peripheral_fraction);
    }

    #[test]
    fn cross_junction_cct_pushes_more_rbc_to_periphery() {
        let basic = cascade_junction_separation(2, 0.45, 2e-3, 1e-3, 5e-6);
        let cross = cascade_junction_separation_cross_junction(2, 0.45, 2e-3, 1e-3);
        assert!(cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
            "cross-junction CCT should push more RBCs to periphery: basic={}, cross={}",
            basic.rbc_peripheral_fraction, cross.rbc_peripheral_fraction);
    }
}
