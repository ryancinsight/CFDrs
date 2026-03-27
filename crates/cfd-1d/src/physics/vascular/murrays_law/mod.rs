//! Murray's Law for optimal vascular bifurcation geometry.
//!
//! ## Theorem: Murray's Law of Minimal Work
//!
//! **Theorem**: The optimal radius r of a vascular segment that minimizes the total
//! biological work W_total required to maintain steady blood flow Q scales as r ∝ Q^(1/3).
//!
//! **Proof Outline**: The total work is the sum of viscous dissipation power (Poiseuille
//! resistance) and the metabolic cost of maintaining the blood volume:
//! ```text
//! W_total = Q² · (8μL / πr⁴) + λ · πr²L
//! ```
//! Setting ∂W/∂r = 0 yields Q ∝ r³.
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`law`] | `MurraysLaw<T>` — diameter calculations, deviation, validation |
//! | [`bifurcation`] | `OptimalBifurcation<T>` — symmetric/asymmetric geometry |
//!
//! ## References
//!
//! - Murray, C.D. (1926) *J. Gen. Physiol.* 9, 835-841.
//! - Sherman, T.F. (1981) *J. Gen. Physiol.* 78, 431-453.
//! - Zamir, M. (1978) *J. Theor. Biol.* 74, 227-250.
//! - Revellin, R. et al. (2009) *Theor. Biol. Med. Model.* 6:7.

pub mod bifurcation;
pub mod law;

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use bifurcation::OptimalBifurcation;
pub use law::MurraysLaw;

// ── Non-Newtonian Flow-Split Exponent ───────────────────────────────────────

/// Non-Newtonian flow-split exponent for power-law fluids.
///
/// ## Theorem (Revellin et al. 2009)
///
/// For a power-law fluid with constitutive relation τ = K·γ̇ⁿ, the volumetric
/// flow rate through a circular tube of radius R under pressure gradient ΔP/L is:
///
/// ```text
/// Q = (nπ/(3n+1)) · (ΔP/(2KL))^(1/n) · R^((3n+1)/n)
/// ```
///
/// At a bifurcation with equal pressure drop across both daughter branches,
/// the flow-split ratio is:
///
/// ```text
/// Q₁/Q₂ = (D₁/D₂)^((3n+1)/n)
/// ```
///
/// The exponent m = (3n+1)/n determines how strongly vessel geometry affects
/// the flow distribution in non-Newtonian fluids:
///
/// | n   | Fluid type                | m = (3n+1)/n |
/// |-----|---------------------------|--------------|
/// | 1.0 | Newtonian                 | 4.0          |
/// | 0.9 | Mildly shear-thinning     | 4.11         |
/// | 0.8 | Moderately shear-thinning | 4.25         |
/// | 0.5 | Strongly shear-thinning   | 5.0          |
///
/// **Reference**: Revellin et al. (2009). *Theor. Biol. Med. Model.* 6:7.
///
/// # Panics
/// Panics if `power_law_index_n` is not positive.
pub fn non_newtonian_flow_split_exponent(power_law_index_n: f64) -> f64 {
    assert!(
        power_law_index_n > 0.0,
        "Power-law index must be positive, got {power_law_index_n}"
    );
    (3.0 * power_law_index_n + 1.0) / power_law_index_n
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_symmetric_daughter_k3() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d_daughter = murray.symmetric_daughter_diameter(d0);

        let expected = d0 / 2.0_f64.powf(1.0 / 3.0);
        assert_relative_eq!(d_daughter, expected, epsilon = 1e-10);
        assert_relative_eq!(d_daughter, 7.937, epsilon = 0.001);
    }

    #[test]
    fn test_murray_deviation_perfect() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d1 = murray.symmetric_daughter_diameter(d0);
        let d2 = d1;
        let deviation = murray.deviation(d0, d1, d2);
        assert!(deviation < 1e-10, "Perfect bifurcation should have zero deviation");
    }

    #[test]
    fn test_murray_deviation_imperfect() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0;
        let d1 = 7.5;
        let d2 = 7.5;
        let deviation = murray.deviation(d0, d1, d2);
        assert!(
            deviation > 0.0 && deviation < 0.2,
            "Deviation {} should be positive but small",
            deviation
        );
    }

    #[test]
    fn test_asymmetric_daughter() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = 10.0_f64;
        let d1 = 9.0_f64;
        let d2 = murray.asymmetric_daughter_diameter(d0, d1).unwrap();
        let lhs = d0.powi(3);
        let rhs = d1.powi(3) + d2.powi(3);
        assert_relative_eq!(lhs, rhs, epsilon = 1e-10);
    }

    #[test]
    fn test_asymmetric_daughter_impossible() {
        let murray = MurraysLaw::<f64>::new();
        let result = murray.asymmetric_daughter_diameter(10.0, 12.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_parent_diameter_reconstruction() {
        let murray = MurraysLaw::<f64>::new();
        let d0 = murray.parent_diameter(7.937, 7.937);
        assert_relative_eq!(d0, 10.0, epsilon = 0.01);
    }

    #[test]
    fn test_ideal_area_ratio_k3() {
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.ideal_area_ratio();
        assert_relative_eq!(ratio, 2.0_f64.powf(1.0 / 3.0), epsilon = 1e-10);
        assert_relative_eq!(ratio, 1.26, epsilon = 0.01);
    }

    #[test]
    fn test_area_preserving() {
        let murray = MurraysLaw::<f64>::area_preserving();
        let ratio = murray.ideal_area_ratio();
        assert_relative_eq!(ratio, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_symmetric_bifurcation() {
        let bif = OptimalBifurcation::<f64>::symmetric(0.01, 1e-6);
        assert!(bif.is_murray_compliant(0.001));
        assert!(bif.mass_conservation_error() < 1e-10);
        assert_relative_eq!(bif.daughter1_diameter, bif.daughter2_diameter, epsilon = 1e-10);
    }

    #[test]
    fn test_asymmetric_bifurcation() {
        let bif = OptimalBifurcation::<f64>::asymmetric(0.01, 1e-6, 0.7);
        assert!(bif.mass_conservation_error() < 1e-10);
        assert_relative_eq!(bif.daughter1_flow / bif.parent_flow, 0.7, epsilon = 1e-10);
        assert!(bif.daughter1_diameter > bif.daughter2_diameter);
    }

    #[test]
    fn test_area_ratio_symmetric() {
        let bif = OptimalBifurcation::<f64>::symmetric(0.01, 1e-6);
        let ratio = bif.area_ratio();
        let murray = MurraysLaw::<f64>::new();
        assert_relative_eq!(ratio, murray.ideal_area_ratio(), epsilon = 0.01);
    }

    #[test]
    fn test_trifurcation_extension() {
        let d0 = 10.0_f64;
        let d1 = 6.0_f64;
        let d2 = 5.0_f64;
        let d3_cubed = d0.powi(3) - d1.powi(3) - d2.powi(3);
        let d3 = d3_cubed.powf(1.0 / 3.0);
        let sum = d1.powi(3) + d2.powi(3) + d3.powi(3);
        assert_relative_eq!(d0.powi(3), sum, epsilon = 1e-10);
    }

    #[test]
    fn test_pressure_drop() {
        let bif = OptimalBifurcation::<f64>::symmetric(0.01, 1e-6);
        let dp = bif.pressure_drop_daughter1(0.0035, 0.1);
        assert!(dp > 0.0 && dp.is_finite());
    }

    // ── Non-Newtonian flow-split exponent tests ─────────────────────────

    #[test]
    fn test_non_newtonian_exponent_newtonian_limit() {
        let m = non_newtonian_flow_split_exponent(1.0);
        assert_relative_eq!(m, 4.0, epsilon = 1e-12);
    }

    #[test]
    fn test_non_newtonian_exponent_shear_thinning() {
        let m = non_newtonian_flow_split_exponent(0.5);
        assert_relative_eq!(m, 5.0, epsilon = 1e-12);
    }

    #[test]
    fn test_non_newtonian_exponent_mildly_shear_thinning() {
        let m = non_newtonian_flow_split_exponent(0.9);
        assert_relative_eq!(m, 3.0 + 1.0 / 0.9, epsilon = 1e-12);
    }

    #[test]
    fn test_non_newtonian_exponent_positive() {
        for &n in &[0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0] {
            let m = non_newtonian_flow_split_exponent(n);
            assert!(m > 0.0, "Exponent must be positive for n={}, got m={}", n, m);
        }
    }

    #[test]
    #[should_panic(expected = "Power-law index must be positive")]
    fn test_non_newtonian_exponent_zero_panics() {
        non_newtonian_flow_split_exponent(0.0);
    }

    #[test]
    fn test_flow_split_newtonian_matches_cubic() {
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.flow_split_ratio(2.0, 1.0, None);
        assert_relative_eq!(ratio, 8.0, epsilon = 1e-10);
    }

    #[test]
    fn test_flow_split_newtonian_via_power_law() {
        let murray = MurraysLaw::<f64>::new();
        let ratio = murray.flow_split_ratio(2.0, 1.0, Some(1.0));
        assert_relative_eq!(ratio, 16.0, epsilon = 1e-10);
    }

    #[test]
    fn test_flow_split_shear_thinning_stronger() {
        let murray = MurraysLaw::<f64>::new();
        let newtonian = murray.flow_split_ratio(2.0, 1.0, Some(1.0));
        let shear_thin = murray.flow_split_ratio(2.0, 1.0, Some(0.5));

        assert_relative_eq!(newtonian, 16.0, epsilon = 1e-10);
        assert_relative_eq!(shear_thin, 32.0, epsilon = 1e-10);
        assert!(shear_thin > newtonian);
    }

    #[test]
    fn test_non_newtonian_constructor_preserves_k3() {
        let murray = MurraysLaw::<f64>::non_newtonian(0.5);
        assert_relative_eq!(murray.exponent, 3.0, epsilon = 1e-12);
        let d0 = 10.0;
        let d_nn = murray.symmetric_daughter_diameter(d0);
        let d_standard = MurraysLaw::<f64>::new().symmetric_daughter_diameter(d0);
        assert_relative_eq!(d_nn, d_standard, epsilon = 1e-12);
    }
}
