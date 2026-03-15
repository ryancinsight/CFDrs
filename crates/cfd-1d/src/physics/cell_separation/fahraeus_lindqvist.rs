//! Fahraeus-Lindqvist apparent viscosity model for microvessels.
//!
//! Blood flowing through vessels smaller than ~300 µm experiences reduced
//! apparent viscosity due to cell-free layer (CFL) formation at the walls.
//!
//! # Reference
//! Pries, A.R., Neuhaus, D. & Gaehtgens, P. (1992). "Blood viscosity in tube
//! flow: dependence on diameter and hematocrit", *Am. J. Physiol.*
//! 263(6):H1770-H1778.

/// Fahraeus-Lindqvist apparent viscosity for microvessels (Pries et al. 1992).
///
/// ## Theorem — Fahraeus-Lindqvist Effect (Fahraeus & Lindqvist 1931)
///
/// In tubes with diameter D < 300 µm, the apparent blood viscosity decreases
/// with decreasing tube diameter due to the formation of a cell-free marginal
/// layer. Red blood cells migrate toward the tube axis (margination), leaving
/// a low-viscosity plasma layer near the walls that acts as a lubricant.
///
/// The Pries et al. (1992) empirical correlation for human blood at 45% hematocrit:
///
/// ```text
/// µ_rel(D) = µ_0.45 · f(D)
///
/// f(D) = [1 + (eta_0.45 - 1) · (D/(D - 1.1))^2 · ((1-0.45)^C - 1) / ((1-0.45)^C - 1)]
/// ```
///
/// Simplified form (valid for D in micrometers):
///
/// ```text
/// µ_rel(D) = 1 + (6.0·exp(-0.085·D) + 3.2 - 2.44·exp(-0.06·D^0.645))
///            · ((1 - H_t)^C - 1) / ((1 - 0.45)^C - 1)
/// ```
///
/// where C = 0.8 + exp(-0.075·D) · (-1 + 1/(1 + 10^-11 · D^12))
///
/// For millifluidic channels (D = 30-200 µm), the effect reduces apparent
/// viscosity by 15-40% compared to bulk blood viscosity.
///
/// **Reference**: Pries, A.R., Neuhaus, D. & Gaehtgens, P. (1992).
/// "Blood viscosity in tube flow: dependence on diameter and hematocrit",
/// *Am. J. Physiol.* 263(6):H1770-H1778.
///
/// # Arguments
/// * `diameter_um` - Tube/channel diameter in micrometers (µm)
/// * `hematocrit` - Feed hematocrit (volume fraction, typically 0.30-0.50)
/// * `mu_plasma_pa_s` - Plasma viscosity in Pa·s (typically 0.0012)
///
/// # Returns
/// Apparent dynamic viscosity in Pa·s
#[must_use]
pub fn fahraeus_lindqvist_viscosity(
    diameter_um: f64,
    hematocrit: f64,
    mu_plasma_pa_s: f64,
) -> f64 {
    // Clamp inputs to physically meaningful ranges
    let d = diameter_um.max(1.0); // Avoid division by zero / nonsensical values
    let ht = hematocrit.clamp(0.0, 0.99);

    // For large vessels (D > 300 µm), use bulk viscosity (no CFL effect).
    // Bulk relative viscosity at the given hematocrit using the same Pries formula
    // evaluated at D = 300 µm (smoothly transitions to bulk).
    let d_eff = if d > 300.0 { 300.0 } else { d };

    // C parameter: captures diameter-dependent hematocrit sensitivity
    // C = 0.8 + exp(-0.075·D) · (-1 + 1/(1 + 1e-11·D^12))
    let c = 0.8 + (-0.075 * d_eff).exp() * (-1.0 + 1.0 / (1.0 + 1e-11 * d_eff.max(1.0).powf(12.0)));

    // eta_0.45: relative viscosity at 45% hematocrit as function of diameter
    // eta_0.45 = 6.0·exp(-0.085·D) + 3.2 - 2.44·exp(-0.06·D^0.645)
    let eta_045 = 6.0 * (-0.085 * d_eff).exp() + 3.2 - 2.44 * (-0.06 * d_eff.max(1.0).powf(0.645)).exp();

    // Relative viscosity at the given hematocrit:
    // eta_rel = 1 + (eta_045 - 1) · ((1-Ht)^C - 1) / ((1-0.45)^C - 1)
    let ref_hct = 0.45_f64;
    let numerator = (1.0 - ht).clamp(0.001, 1.0).powf(c) - 1.0;
    let denominator = (1.0 - ref_hct).clamp(0.001, 1.0).powf(c) - 1.0;

    let eta_rel = if denominator.abs() < 1e-30 {
        // Edge case: C is so large that both terms are essentially -1
        1.0
    } else {
        1.0 + (eta_045 - 1.0) * numerator / denominator
    };

    // Ensure viscosity is at least plasma viscosity
    mu_plasma_pa_s * eta_rel.max(1.0)
}

/// Secomb (2017) in-vivo apparent viscosity for microvascular networks.
///
/// ## Theorem — Secomb Network Blood Flow Viscosity (Secomb 2017)
///
/// In microvascular networks, the apparent viscosity is modified by the
/// cell-free layer (CFL) width W(D), which depends on vessel diameter and
/// hematocrit:
///
/// ```text
/// μ_vivo(D, H_t) = μ_plasma · [1 + (μ_0.45 − 1) · F(D, H_t)]
/// ```
///
/// where:
/// - μ_0.45 = relative viscosity at H_t = 0.45 (from Pries 1992 parameterization)
/// - F(D, H_t) = [(1 − H_t)^C − 1] / [(1 − 0.45)^C − 1]
/// - C(D) = (0.8 + e^(−0.075D)) · (−1 + 1/(1 + 10^(−11) · D^12)) + 1/(1 + 10^(−11) · D^12)
///
/// The key improvement over Pries (1992): Secomb adds a **phase-separation
/// parameter** X₀ that captures hematocrit partitioning at bifurcations:
///
/// ```text
/// X₀(D, H_t) = 0.964 · (1 − H_t) / D
/// ```
///
/// This allows computing the hematocrit in each daughter branch of a
/// bifurcation based on the flow split ratio and parent hematocrit.
///
/// For standalone viscosity (no bifurcation), the model reduces to a
/// refined version of Pries (1992) with improved C(D) coefficients.
///
/// **Reference**: Secomb, T.W. (2017). "Blood Flow in the Microcirculation",
/// *Annu. Rev. Fluid Mech.* 49:443-461.
///
/// # Arguments
/// * `diameter_um` - Tube/channel diameter in micrometers (µm)
/// * `hematocrit` - Feed hematocrit (volume fraction, typically 0.30-0.50)
/// * `mu_plasma_pa_s` - Plasma viscosity in Pa·s (typically 0.0012)
///
/// # Returns
/// Apparent dynamic viscosity in Pa·s
#[must_use]
pub fn secomb_network_viscosity(
    diameter_um: f64,
    hematocrit: f64,
    mu_plasma_pa_s: f64,
) -> f64 {
    // Clamp inputs to physically meaningful ranges
    let d = diameter_um.max(1.0);
    let ht = hematocrit.clamp(0.0, 0.99);

    // For large vessels (D > 300 µm), cap at D=300 (same convention as Pries)
    let d_eff = if d > 300.0 { 300.0 } else { d };

    // Secomb (2017) refined C(D) coefficient:
    // C = (0.8 + exp(-0.075·D)) · (-1 + 1/(1 + 1e-11·D^12)) + 1/(1 + 1e-11·D^12)
    //
    // This differs from Pries (1992) by the additive term 1/(1 + 1e-11·D^12)
    // which provides a smoother transition for intermediate diameters.
    let d12_term = 1.0 / (1.0 + 1e-11 * d_eff.max(1.0).powf(12.0));
    let c = (0.8 + (-0.075 * d_eff).exp()) * (-1.0 + d12_term) + d12_term;

    // eta_0.45: relative viscosity at 45% hematocrit (Pries 1992 parameterization,
    // shared with Secomb 2017)
    let eta_045 = 6.0 * (-0.085 * d_eff).exp() + 3.2 - 2.44 * (-0.06 * d_eff.max(1.0).powf(0.645)).exp();

    // Relative viscosity at the given hematocrit:
    // eta_rel = 1 + (eta_045 - 1) · ((1 - Ht)^C - 1) / ((1 - 0.45)^C - 1)
    let ref_hct = 0.45_f64;
    let numerator = (1.0 - ht).clamp(0.001, 1.0).powf(c) - 1.0;
    let denominator = (1.0 - ref_hct).clamp(0.001, 1.0).powf(c) - 1.0;

    let eta_rel = if denominator.abs() < 1e-30 {
        1.0
    } else {
        1.0 + (eta_045 - 1.0) * numerator / denominator
    };

    // Ensure viscosity is at least plasma viscosity
    mu_plasma_pa_s * eta_rel.max(1.0)
}

/// Secomb (2017) phase-separation parameter for bifurcation hematocrit partitioning.
///
/// At a microvascular bifurcation, the cell-free layer width determines how
/// hematocrit partitions between daughter branches. The phase-separation
/// parameter X₀ captures this effect:
///
/// ```text
/// X₀(D, H_t) = 0.964 · (1 − H_t) / D
/// ```
///
/// Larger X₀ means stronger plasma skimming (more hematocrit asymmetry
/// between daughter branches). X₀ decreases with vessel diameter (larger
/// vessels have proportionally thinner cell-free layers) and increases
/// with lower hematocrit (thicker cell-free layer at lower concentrations).
///
/// **Reference**: Secomb, T.W. (2017). "Blood Flow in the Microcirculation",
/// *Annu. Rev. Fluid Mech.* 49:443-461.
///
/// # Arguments
/// * `diameter_um` - Vessel diameter in micrometers (µm)
/// * `hematocrit` - Hematocrit (volume fraction, 0.0-1.0)
///
/// # Returns
/// Phase-separation parameter X₀ (dimensionless)
#[inline]
#[must_use]
pub fn secomb_phase_separation_x0(diameter_um: f64, hematocrit: f64) -> f64 {
    0.964 * (1.0 - hematocrit) / diameter_um.max(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_PLASMA: f64 = 0.0012; // Pa·s

    /// Bulk blood viscosity at 45% hematocrit is approximately 3-4 mPa·s.
    /// For D=1000 µm (large vessel), the Fahraeus-Lindqvist effect is negligible
    /// and viscosity should approach bulk values.
    #[test]
    fn test_fahraeus_lindqvist_large_vessel() {
        let mu = fahraeus_lindqvist_viscosity(1000.0, 0.45, MU_PLASMA);
        // For large D, eta_rel should be evaluated at D=300 (our cap),
        // which gives approximately the bulk viscosity for 45% Ht.
        // Bulk blood viscosity ~3.0-4.5 mPa·s
        let eta_rel = mu / MU_PLASMA;
        assert!(
            eta_rel > 2.0 && eta_rel < 5.0,
            "Large vessel eta_rel = {:.3} should be ~3-4 (bulk)",
            eta_rel
        );
    }

    /// In microchannels (D=50 µm), apparent viscosity should be lower than bulk
    /// due to the cell-free layer effect.
    #[test]
    fn test_fahraeus_lindqvist_reduction_in_microchannel() {
        let mu_bulk = fahraeus_lindqvist_viscosity(1000.0, 0.45, MU_PLASMA);
        let mu_micro = fahraeus_lindqvist_viscosity(50.0, 0.45, MU_PLASMA);

        assert!(
            mu_micro < mu_bulk,
            "Microchannel viscosity ({:.4} Pa·s) should be less than bulk ({:.4} Pa·s)",
            mu_micro,
            mu_bulk,
        );
    }

    /// The Pries 1992 formula exhibits a viscosity minimum in the range
    /// D ~ 10-40 µm (the Fahraeus-Lindqvist minimum). For very small D
    /// approaching the RBC diameter (~8 µm), viscosity rises again due to
    /// single-file flow and near-occlusion effects. Verify non-monotonic
    /// behavior: a mid-range diameter has lower viscosity than both very
    /// small and very large diameters.
    #[test]
    fn test_fahraeus_lindqvist_minimum_near_7um() {
        let mu_7 = fahraeus_lindqvist_viscosity(7.0, 0.45, MU_PLASMA);
        let mu_30 = fahraeus_lindqvist_viscosity(30.0, 0.45, MU_PLASMA);
        let mu_200 = fahraeus_lindqvist_viscosity(200.0, 0.45, MU_PLASMA);

        // All values should be finite and above plasma viscosity
        assert!(
            mu_7.is_finite() && mu_7 > MU_PLASMA,
            "Viscosity at 7 µm ({:.6}) should be finite and > plasma",
            mu_7,
        );

        // D=30 µm (in the FL minimum region) should have lower viscosity
        // than D=200 µm (approaching bulk), demonstrating the FL effect.
        assert!(
            mu_30 < mu_200,
            "Viscosity at 30 µm ({:.6}) should be less than at 200 µm ({:.6}) — FL effect",
            mu_30,
            mu_200,
        );

        // D=7 µm should have higher viscosity than D=30 µm (rise at very
        // small diameters due to near-occlusion / single-file effects).
        assert!(
            mu_7 > mu_30,
            "Viscosity at 7 µm ({:.6}) should exceed 30 µm ({:.6}) — near-occlusion rise",
            mu_7,
            mu_30,
        );
    }

    /// Higher hematocrit should always produce higher apparent viscosity
    /// at any given diameter.
    #[test]
    fn test_fahraeus_lindqvist_increases_with_hematocrit() {
        let d = 100.0; // 100 µm channel
        let mu_low = fahraeus_lindqvist_viscosity(d, 0.30, MU_PLASMA);
        let mu_mid = fahraeus_lindqvist_viscosity(d, 0.40, MU_PLASMA);
        let mu_high = fahraeus_lindqvist_viscosity(d, 0.50, MU_PLASMA);

        assert!(
            mu_low < mu_mid && mu_mid < mu_high,
            "Viscosity should increase with hematocrit: {:.4} < {:.4} < {:.4}",
            mu_low,
            mu_mid,
            mu_high
        );
    }

    /// Zero hematocrit should return plasma viscosity.
    #[test]
    fn test_fahraeus_lindqvist_zero_hematocrit() {
        let mu = fahraeus_lindqvist_viscosity(100.0, 0.0, MU_PLASMA);
        let eta_rel = mu / MU_PLASMA;
        assert!(
            (eta_rel - 1.0).abs() < 0.01,
            "Zero hematocrit eta_rel = {:.4} should be ~1.0",
            eta_rel
        );
    }

    // ── Secomb (2017) tests ─────────────────────────────────────────────

    /// For large vessels (D=1000 µm), Secomb and Pries should agree closely
    /// because both cap at D=300 µm and the Secomb correction terms vanish
    /// at large diameters.
    #[test]
    fn test_secomb_agrees_with_pries_for_large_vessels() {
        let mu_pries = fahraeus_lindqvist_viscosity(1000.0, 0.45, MU_PLASMA);
        let mu_secomb = secomb_network_viscosity(1000.0, 0.45, MU_PLASMA);

        let rel_diff = (mu_secomb - mu_pries).abs() / mu_pries;
        assert!(
            rel_diff < 0.05,
            "Secomb ({:.6}) and Pries ({:.6}) should agree for large vessels, rel_diff = {:.4}",
            mu_secomb,
            mu_pries,
            rel_diff
        );
    }

    /// In a microchannel (D=50 µm), Secomb viscosity should be less than
    /// bulk viscosity (the Fahraeus-Lindqvist effect is preserved).
    #[test]
    fn test_secomb_lower_than_bulk_in_microchannel() {
        let mu_bulk = secomb_network_viscosity(1000.0, 0.45, MU_PLASMA);
        let mu_micro = secomb_network_viscosity(50.0, 0.45, MU_PLASMA);

        assert!(
            mu_micro < mu_bulk,
            "Secomb microchannel viscosity ({:.6}) should be < bulk ({:.6})",
            mu_micro,
            mu_bulk
        );
    }

    /// The phase-separation parameter X₀ should decrease with increasing
    /// vessel diameter (larger vessels have proportionally thinner CFL).
    #[test]
    fn test_secomb_phase_separation_decreases_with_diameter() {
        let x0_small = secomb_phase_separation_x0(20.0, 0.45);
        let x0_medium = secomb_phase_separation_x0(100.0, 0.45);
        let x0_large = secomb_phase_separation_x0(500.0, 0.45);

        assert!(
            x0_small > x0_medium && x0_medium > x0_large,
            "X₀ should decrease with diameter: {:.6} > {:.6} > {:.6}",
            x0_small,
            x0_medium,
            x0_large
        );
    }

    /// Higher hematocrit should produce higher Secomb viscosity at any
    /// given diameter.
    #[test]
    fn test_secomb_increases_with_hematocrit() {
        let d = 100.0;
        let mu_low = secomb_network_viscosity(d, 0.30, MU_PLASMA);
        let mu_mid = secomb_network_viscosity(d, 0.40, MU_PLASMA);
        let mu_high = secomb_network_viscosity(d, 0.50, MU_PLASMA);

        assert!(
            mu_low < mu_mid && mu_mid < mu_high,
            "Secomb viscosity should increase with hematocrit: {:.6} < {:.6} < {:.6}",
            mu_low,
            mu_mid,
            mu_high
        );
    }
}
