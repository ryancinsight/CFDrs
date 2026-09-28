use super::*;
use aequitas::systems::si::quantities::{DynamicViscosity, Length, VolumetricFlowRate};
use eunomia::assert_relative_eq;

#[test]
fn test_symmetric_daughter_k3() {
    let murray = MurraysLaw::<f64>::new();
    let d0 = 10.0;
    let d_daughter = murray.symmetric_daughter_diameter(d0);

    let expected = d0 / 2.0_f64.cbrt();
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
    assert!(
        deviation < 1e-10,
        "Perfect bifurcation should have zero deviation"
    );
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
        "Deviation {deviation} should be positive but small"
    );
}

#[test]
fn test_asymmetric_daughter() {
    let murray = MurraysLaw::<f64>::new();
    let d0 = 10.0_f64;
    let d1 = 9.0_f64;
    let d2 = murray
        .asymmetric_daughter_diameter(d0, d1)
        .expect("expected value");
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
    assert_relative_eq!(ratio, 2.0_f64.cbrt(), epsilon = 1e-10);
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
    let bif = OptimalBifurcation::<f64>::symmetric(
        Length::from_base(0.01),
        VolumetricFlowRate::from_base(1e-6),
    );
    assert!(bif.is_murray_compliant(0.001));
    assert!(bif.mass_conservation_error().into_base() < 1e-10);
    assert_relative_eq!(
        bif.daughter1_diameter.into_base(),
        bif.daughter2_diameter.into_base(),
        epsilon = 1e-10
    );
}

#[test]
fn test_asymmetric_bifurcation() {
    let bif = OptimalBifurcation::<f64>::asymmetric(
        Length::from_base(0.01),
        VolumetricFlowRate::from_base(1e-6),
        0.7,
    );
    assert!(bif.mass_conservation_error().into_base() < 1e-10);
    assert_relative_eq!(
        bif.daughter1_flow.into_base() / bif.parent_flow.into_base(),
        0.7,
        epsilon = 1e-10
    );
    assert!(bif.daughter1_diameter.into_base() > bif.daughter2_diameter.into_base());
}

#[test]
fn test_area_ratio_symmetric() {
    let bif = OptimalBifurcation::<f64>::symmetric(
        Length::from_base(0.01),
        VolumetricFlowRate::from_base(1e-6),
    );
    let ratio = bif.area_ratio().into_base();
    let murray = MurraysLaw::<f64>::new();
    assert_relative_eq!(ratio, murray.ideal_area_ratio(), epsilon = 0.01);
}

#[test]
fn test_trifurcation_extension() {
    let d0 = 10.0_f64;
    let d1 = 6.0_f64;
    let d2 = 5.0_f64;
    let d3_cubed = d0.powi(3) - d1.powi(3) - d2.powi(3);
    let d3 = d3_cubed.cbrt();
    let sum = d1.powi(3) + d2.powi(3) + d3.powi(3);
    assert_relative_eq!(d0.powi(3), sum, epsilon = 1e-10);
}

#[test]
fn test_pressure_drop() {
    let bif = OptimalBifurcation::<f64>::symmetric(
        Length::from_base(0.01),
        VolumetricFlowRate::from_base(1e-6),
    );
    let dp =
        bif.pressure_drop_daughter1(DynamicViscosity::from_base(0.0035), Length::from_base(0.1));
    assert!(dp.into_base() > 0.0 && dp.into_base().is_finite());
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
        assert!(m > 0.0, "Exponent must be positive for n={n}, got m={m}");
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
