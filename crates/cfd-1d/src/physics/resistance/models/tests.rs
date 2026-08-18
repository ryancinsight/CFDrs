use super::traits::{FlowConditions, ResistanceModel};
use super::{darcy_weisbach::DarcyWeisbachModel, hagen_poiseuille::HagenPoiseuilleModel};
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, SpecificHeatCapacity, ThermalConductivity, Velocity,
};
use cfd_core::physics::fluid::ConstantPropertyFluid;

#[test]
fn test_hagen_poiseuille_resistance_matches_formula() {
    // Parameters: μ=1e-3 Pa·s (water), L=1 m, D=0.01 m
    let fluid = ConstantPropertyFluid::new(
        "Test Water".to_string(),
        MassDensity::from_base(1000.0),
        DynamicViscosity::from_base(1.0e-3),
        SpecificHeatCapacity::from_base(4186.0),
        ThermalConductivity::from_base(0.6),
        Velocity::from_base(1500.0),
    );
    let model = HagenPoiseuilleModel::new(0.01f64, 1.0f64);
    let cond = FlowConditions::new(0.0f64);

    let r = model
        .calculate_resistance(&fluid, &cond)
        .expect("expected value");
    // R = 128 μ L / (π D^4)
    let expected =
        128.0f64 * fluid.viscosity.into_base() * 1.0f64 / (std::f64::consts::PI * 0.01f64.powi(4));
    let rel_err = (r - expected).abs() / expected.abs();
    assert!(
        rel_err < 1e-12,
        "Hagen-Poiseuille resistance mismatch: r={r} expected={expected} rel_err={rel_err}"
    );
}

#[test]
fn test_darcy_weisbach_laminar_limit_friction_factor() {
    // Laminar regime: Re=1000 < 2300, expect f ≈ 64/Re
    let fluid = ConstantPropertyFluid::new(
        "Test Water".to_string(),
        MassDensity::from_base(1000.0),
        DynamicViscosity::from_base(1.0e-3),
        SpecificHeatCapacity::from_base(4186.0),
        ThermalConductivity::from_base(0.6),
        Velocity::from_base(1500.0),
    );
    let model = DarcyWeisbachModel::circular(0.01f64, 1.0f64, 0.0f64);
    let mut cond = FlowConditions::new(0.0f64);
    cond.reynolds_number = Some(1000.0f64);

    let r = model
        .calculate_resistance(&fluid, &cond)
        .expect("expected value");
    let area = std::f64::consts::PI * (0.01f64.powi(2)) / 4.0f64;
    // Darcy-Weisbach with f=64/Re and Re=rho*V*D/mu reduces to
    // R = 32*mu*L/(A*D^2), the linear laminar resistance.
    let expected = 32.0f64 * fluid.viscosity.into_base() * 1.0f64 / (area * 0.01f64.powi(2));
    let rel_err = (r - expected).abs() / expected.abs();
    assert!(
        rel_err < 1e-3,
        "Darcy-Weisbach laminar limit mismatch: r={r} expected={expected} rel_err={rel_err}"
    );
}
