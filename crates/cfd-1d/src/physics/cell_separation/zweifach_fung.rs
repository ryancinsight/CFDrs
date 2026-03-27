//! Zweifach-Fung bifurcation law for RBC partitioning at microvascular bifurcations.
//!
//! # Theorem — Zweifach-Fung Effect (Svanes & Zweifach 1968, Fung 1973)
//!
//! At a microvascular bifurcation, red blood cells preferentially enter the
//! daughter branch with the higher flow rate.  Below a critical flow ratio
//! $Q_{high}/Q_{low} \approx 2.5$, both branches receive RBCs in proportion
//! to their flow.  Above this critical ratio, the low-flow branch experiences
//! near-complete plasma separation (the "all-or-nothing" phenomenon).
//!
//! ## Governing equation
//!
//! The fractional RBC flux entering the high-flow daughter is modelled as a
//! sigmoid transition centred on the critical flow ratio:
//!
//! ```text
//! f_RBC(Q_r) = 1 / (1 + exp(-k · (Q_r − Q_r,crit)))
//! ```
//!
//! where:
//! - $Q_r = Q_{high} / Q_{total}$ is the fractional flow to the high-flow daughter
//! - $Q_{r,crit}$ is the critical fractional flow for all-or-nothing capture
//! - $k$ controls the sharpness of the transition (depends on confinement ratio)
//!
//! ## Confinement dependence
//!
//! The critical flow ratio and transition sharpness depend on the confinement
//! ratio $\kappa = d_{RBC} / D_{channel}$:
//!
//! ```text
//! Q_r,crit = 0.5 + 0.14 · κ^{−0.5}       (Doyeux et al. 2011)
//! k        = 10 + 40 · κ²                   (higher confinement = sharper transition)
//! ```
//!
//! At strong confinement ($\kappa > 0.5$), the transition is nearly step-like
//! (single-file regime); at weak confinement ($\kappa < 0.1$), the transition
//! is gradual and the Pries logit model dominates.
//!
//! ## Relationship to plasma skimming
//!
//! The Zweifach-Fung model complements the Pries (1990) plasma skimming logit
//! model in `plasma_skimming.rs`:
//! - **Pries logit**: continuous hematocrit partitioning, valid for all flow ratios
//! - **Zweifach-Fung**: captures the critical-ratio threshold and all-or-nothing
//!   behaviour at high confinement
//!
//! # References
//! - Svanes, K. & Zweifach, B.W. (1968). Variations in small blood vessel
//!   hematocrits produced in hypothermic rats by micro-occlusion.
//!   *Microvasc. Res.* 1:210-220.
//! - Fung, Y.C. (1973). Stochastic flow in capillary blood vessels.
//!   *Microvasc. Res.* 5:34-48.
//! - Doyeux, V. et al. (2011). Spheres in the vicinity of a bifurcation:
//!   elucidating the Zweifach-Fung effect. *J. Fluid Mech.* 674:359-388.

/// Mean human RBC diameter [µm].
const RBC_DIAMETER_UM: f64 = 8.0;

/// Confinement ratio: $\kappa = d_{RBC} / D_{channel}$.
///
/// Higher confinement means RBCs are large relative to the channel, increasing
/// the strength of the Zweifach-Fung effect.
///
/// # Arguments
/// * `channel_diameter_um` — Channel hydraulic diameter [µm]
///
/// # Returns
/// Confinement ratio $\kappa \in (0, 1]$, clamped to a minimum channel
/// diameter of 1 µm to avoid singularity.
#[inline]
#[must_use]
pub fn confinement_ratio(channel_diameter_um: f64) -> f64 {
    RBC_DIAMETER_UM / channel_diameter_um.max(RBC_DIAMETER_UM)
}

/// Critical fractional flow for the Zweifach-Fung all-or-nothing transition.
///
/// # Theorem
///
/// The critical fractional flow $Q_{r,crit}$ above which the high-flow
/// daughter branch captures nearly all RBCs is:
///
/// ```text
/// Q_r,crit = 0.5 + 0.14 · κ^{-0.5}
/// ```
///
/// where $\kappa = d_{RBC}/D$ is the confinement ratio.
///
/// At $\kappa \to 0$ (large channels), $Q_{r,crit} \to \infty$ (no ZF effect;
/// clamped to 0.99).  At $\kappa = 1$ (single-file), $Q_{r,crit} \approx 0.64$.
///
/// **Derivation**: Doyeux et al. (2011) show via boundary-element simulations
/// that the critical ratio transitions from ~2.5:1 at $\kappa = 0.3$ to
/// ~1.8:1 at $\kappa = 0.7$, consistent with this parameterisation.
///
/// # Arguments
/// * `channel_diameter_um` — Channel hydraulic diameter [µm]
///
/// # Returns
/// Critical fractional flow $Q_{r,crit} \in [0.5, 0.99]$.
#[inline]
#[must_use]
pub fn critical_fractional_flow(channel_diameter_um: f64) -> f64 {
    let kappa = confinement_ratio(channel_diameter_um);
    // For kappa → 0 (very large channels), the 1/sqrt(kappa) term diverges,
    // so we clamp output to [0.5, 0.99].
    (0.5 + 0.14 * kappa.powf(-0.5)).clamp(0.5, 0.99)
}

/// Transition sharpness parameter for the sigmoid model.
///
/// Higher confinement produces a sharper (more step-like) transition
/// between the "shared" and "all-or-nothing" regimes.
///
/// ```text
/// k = 10 + 40 · κ²
/// ```
///
/// At $\kappa = 1$: $k = 50$ (sharp step).
/// At $\kappa = 0.1$: $k = 10.4$ (gradual).
#[inline]
#[must_use]
fn transition_sharpness(channel_diameter_um: f64) -> f64 {
    let kappa = confinement_ratio(channel_diameter_um);
    10.0 + 40.0 * kappa * kappa
}

/// Compute the fractional RBC flux entering the high-flow daughter branch
/// at a bifurcation, using the Zweifach-Fung sigmoid transition model.
///
/// # Theorem — Zweifach-Fung RBC Partitioning
///
/// The fraction of RBCs entering the higher-flow daughter is:
///
/// ```text
/// f_RBC(Q_r) = 1 / (1 + exp(-k · (Q_r − Q_r,crit)))
/// ```
///
/// **Properties**:
/// - $f_{RBC}(Q_{r,crit}) = 0.5$ (equal partition at the critical point)
/// - $f_{RBC} \to 1$ as $Q_r \to 1$ (all RBCs to high-flow branch)
/// - $f_{RBC} \to 0$ as $Q_r \to 0$ (no RBCs to this branch)
///
/// **Proof of boundedness**: The sigmoid $1/(1+e^{-x})$ maps $\mathbb{R} \to (0,1)$,
/// so $f_{RBC} \in (0, 1)$ for all finite $Q_r$.  ∎
///
/// # Arguments
/// * `flow_fraction` — Fractional volumetric flow to this daughter branch,
///   $Q_{daughter}/Q_{total} \in [0, 1]$
/// * `channel_diameter_um` — Channel hydraulic diameter [µm]
///
/// # Returns
/// Fractional RBC flux to this daughter branch $\in (0, 1)$.
#[inline]
#[must_use]
pub fn zweifach_fung_rbc_fraction(flow_fraction: f64, channel_diameter_um: f64) -> f64 {
    let qr = flow_fraction.clamp(0.0, 1.0);
    let kappa = confinement_ratio(channel_diameter_um);
    let qr_crit = critical_fractional_flow(channel_diameter_um);
    let k = transition_sharpness(channel_diameter_um);

    // Sigmoid RBC fraction (ZF-dominated)
    let f_sigmoid = 1.0 / (1.0 + (-k * (qr - qr_crit)).exp());

    // Confinement weight: controls blending between sigmoid and linear.
    // At κ=1 (single-file), weight=1 → pure sigmoid.
    // At κ→0 (large channel), weight→0 → pure linear (proportional).
    // Smooth transition via κ² to avoid abrupt cutoff.
    let zf_weight = (kappa * kappa).clamp(0.0, 1.0);

    // Blend: linear (proportional) at low confinement, sigmoid at high
    zf_weight * f_sigmoid + (1.0 - zf_weight) * qr
}

/// Compute both daughter hematocrits from Zweifach-Fung partitioning.
///
/// Given a parent vessel with feed hematocrit $H_F$ splitting into two
/// daughters with fractional flows $f_1$ and $f_2 = 1 - f_1$, this function
/// returns $(H_1, H_2)$ such that RBC mass is conserved:
///
/// ```text
/// H_F · Q_total = H_1 · Q_1 + H_2 · Q_2
/// ```
///
/// # Arguments
/// * `feed_hematocrit` — Parent vessel hematocrit $\in [0, 1]$
/// * `flow_fraction_1` — Fractional flow to daughter 1, $Q_1/Q_{total} \in [0, 1]$
/// * `channel_diameter_um` — Channel hydraulic diameter [µm]
///
/// # Returns
/// `(h1, h2)` — Hematocrits of daughter branches 1 and 2, both $\in [0, 1]$.
#[must_use]
pub fn zweifach_fung_daughter_hematocrits(
    feed_hematocrit: f64,
    flow_fraction_1: f64,
    channel_diameter_um: f64,
) -> (f64, f64) {
    let ht = feed_hematocrit.clamp(0.0, 1.0);
    let f1 = flow_fraction_1.clamp(0.0, 1.0);
    let f2 = 1.0 - f1;

    if ht < 1e-15 {
        return (0.0, 0.0);
    }
    if f1 < 1e-15 {
        return (0.0, ht);
    }
    if f2 < 1e-15 {
        return (ht, 0.0);
    }

    // Fractional RBC flux to daughter 1
    let rbc_frac_1 = zweifach_fung_rbc_fraction(f1, channel_diameter_um);
    // RBC mass balance: H_1 · f1 = H_F · rbc_frac_1  →  H_1 = H_F · rbc_frac_1 / f1
    let h1 = (ht * rbc_frac_1 / f1).clamp(0.0, 1.0);
    // Conservation: H_F = H_1·f1 + H_2·f2  →  H_2 = (H_F - H_1·f1) / f2
    let h2 = ((ht - h1 * f1) / f2).clamp(0.0, 1.0);

    (h1, h2)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Symmetric split (50/50) in large channel: both daughters should receive
    /// approximately equal hematocrit (no Zweifach-Fung bias).
    #[test]
    fn symmetric_split_equal_hematocrit() {
        let (h1, h2) = zweifach_fung_daughter_hematocrits(0.45, 0.5, 200.0);
        assert!(
            (h1 - h2).abs() < 0.05,
            "Symmetric split: h1={h1:.4}, h2={h2:.4} should be approximately equal"
        );
        // Mass conservation
        let mass_err = (0.45 - (h1 * 0.5 + h2 * 0.5)).abs();
        assert!(mass_err < 0.01, "Mass conservation error: {mass_err:.6}");
    }

    /// Asymmetric split (90/10) in strongly confined channel: high-flow branch
    /// should receive disproportionately more RBCs.
    #[test]
    fn asymmetric_split_high_flow_gets_more_rbcs() {
        // D=15 µm → κ=0.53, strong ZF effect
        let (h1, h2) = zweifach_fung_daughter_hematocrits(0.45, 0.9, 15.0);
        assert!(
            h1 > h2,
            "High-flow daughter h1={h1:.4} should exceed low-flow h2={h2:.4}"
        );
        // The RBC fraction to the high-flow branch should exceed the flow fraction
        let rbc_frac = zweifach_fung_rbc_fraction(0.9, 15.0);
        assert!(
            rbc_frac > 0.9,
            "RBC fraction {rbc_frac:.4} should exceed flow fraction 0.9"
        );
    }

    /// At extreme flow ratio (99/1) in highly confined channel: near-complete
    /// capture (all-or-nothing).
    #[test]
    fn extreme_ratio_all_or_nothing() {
        // D=10 µm → κ=0.8, strong ZF effect
        let rbc_frac = zweifach_fung_rbc_fraction(0.99, 10.0);
        assert!(
            rbc_frac > 0.99,
            "Extreme flow ratio should give near-complete RBC capture, got {rbc_frac:.6}"
        );
    }

    /// Zero flow fraction means zero RBC flux.
    #[test]
    fn zero_flow_zero_rbcs() {
        let (h1, _h2) = zweifach_fung_daughter_hematocrits(0.45, 0.0, 50.0);
        assert!(h1.abs() < 1e-10, "Zero flow → zero hematocrit, got {h1:.10}");
    }

    /// Full flow fraction means all RBCs go to this branch.
    #[test]
    fn full_flow_all_rbcs() {
        let (h1, h2) = zweifach_fung_daughter_hematocrits(0.45, 1.0, 50.0);
        assert!(
            (h1 - 0.45).abs() < 0.01,
            "Full flow should give feed hematocrit, got {h1:.4}"
        );
        assert!(h2.abs() < 1e-10, "No flow → zero hematocrit, got {h2:.10}");
    }

    /// Confinement ratio is bounded [0, 1].
    #[test]
    fn confinement_ratio_bounded() {
        assert!((0.0..=1.0).contains(&confinement_ratio(1000.0)));
        assert!((0.0..=1.0).contains(&confinement_ratio(8.0)));
        assert!((0.0..=1.0).contains(&confinement_ratio(5.0)));
    }

    /// Critical fractional flow decreases with confinement (smaller channels
    /// have lower critical ratio → easier to trigger all-or-nothing).
    #[test]
    fn critical_flow_decreases_with_confinement() {
        let qr_large = critical_fractional_flow(200.0);
        let qr_small = critical_fractional_flow(20.0);
        assert!(
            qr_small < qr_large,
            "Critical flow in small channel ({qr_small:.4}) should be < large ({qr_large:.4})"
        );
    }

    /// RBC fraction is monotonically increasing with flow fraction.
    #[test]
    fn rbc_fraction_monotonic() {
        let d = 50.0;
        let mut prev = zweifach_fung_rbc_fraction(0.01, d);
        for fq in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99] {
            let current = zweifach_fung_rbc_fraction(fq, d);
            assert!(
                current >= prev,
                "RBC fraction should increase: f({fq})={current:.6} < f(prev)={prev:.6}"
            );
            prev = current;
        }
    }

    /// Daughter hematocrits are always in [0, 1].
    #[test]
    fn daughter_hematocrits_bounded() {
        let cases = [
            (0.45, 0.1, 20.0),
            (0.45, 0.5, 100.0),
            (0.45, 0.9, 50.0),
            (0.80, 0.3, 30.0),
            (0.10, 0.7, 200.0),
        ];
        for &(ht, fq, d) in &cases {
            let (h1, h2) = zweifach_fung_daughter_hematocrits(ht, fq, d);
            assert!(
                (0.0..=1.0).contains(&h1) && (0.0..=1.0).contains(&h2),
                "Ht={ht}, fq={fq}, d={d}: h1={h1:.4}, h2={h2:.4} out of [0,1]"
            );
        }
    }
}
