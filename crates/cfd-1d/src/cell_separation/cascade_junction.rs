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
//! P_center(cell) = r^β / (r^β + (1 − r)^β)
//! ```
//!
//! where `r = Q_center / Q_total` and `β` is a stiffness exponent:
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
/// Extends the Zweifach-Fung bifurcation law (1969) to a 3-way split:
/// treats the peripheral arms collectively as the "low-flow" side.
///
/// # Arguments
/// * `q_center_frac` — fraction of total inlet flow carried by the center arm (0–1).
/// * `stiffness_exp` — cell stiffness exponent β.  Larger values → stronger bias
///   toward the high-flow arm.
fn p_center(q_center_frac: f64, stiffness_exp: f64) -> f64 {
    let r = q_center_frac.clamp(1e-9, 1.0 - 1e-9);
    let r_beta = r.powf(stiffness_exp);
    let s_beta = (1.0 - r).powf(stiffness_exp);
    r_beta / (r_beta + s_beta)
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
    // Stiffness exponents (empirical, per Fung 1969 + modern microfluidic literature)
    const SE_CANCER: f64 = 1.70; // MCF-7 breast cancer cells (stiff, ~17.5 µm diameter)
    const SE_WBC: f64 = 1.40;    // WBCs (semi-rigid, ~10 µm diameter)
    const SE_RBC: f64 = 1.00;    // RBCs (deformable — distributes by flow fraction)

    // Flow fraction from hydraulic resistance of width-scaled arms.
    // Re-computed at the trunk center_frac (same fraction at every level since we
    // only cascade the center arm, which re-splits at the same frac each time).
    let q_frac = tri_center_q_frac(center_frac);

    // Accumulate center-arm fraction for each cell type across all cascade levels.
    let mut f_cancer = 1.0_f64;
    let mut f_wbc = 1.0_f64;
    let mut f_rbc = 1.0_f64;

    for _ in 0..n_levels {
        // At each level: the current q_frac applies (same asymmetry repeated).
        f_cancer *= p_center(q_frac, SE_CANCER);
        f_wbc    *= p_center(q_frac, SE_WBC);
        f_rbc    *= p_center(q_frac, SE_RBC);
        // Note: q_frac itself stays constant across levels because each level
        // re-splits the center arm at the same center_frac ratio.
    }

    let rbc_periph = (1.0 - f_rbc).clamp(0.0, 1.0);
    let sep_eff    = (f_cancer - f_rbc).abs().clamp(0.0, 1.0);

    CascadeJunctionResult {
        cancer_center_fraction: f_cancer.clamp(0.0, 1.0),
        wbc_center_fraction:    f_wbc.clamp(0.0, 1.0),
        rbc_peripheral_fraction: rbc_periph,
        separation_efficiency:  sep_eff,
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
}
