use super::traits::{FlowConditions, ResistanceModel};
use super::{hagen_poiseuille::HagenPoiseuilleModel, darcy_weisbach::DarcyWeisbachModel};
use cfd_core::fluid::Fluid;

#[test]
fn test_hagen_poiseuille_resistance_matches_formula() {
    // Parameters: μ=1e-3 Pa·s (water), L=1 m, D=0.01 m
    let fluid = Fluid { density: 1000.0f64, viscosity: 1.0e-3f64 };
    let model = HagenPoiseuilleModel::new(0.01f64, 1.0f64);
    let cond = FlowConditions::new(0.0f64);

    let r = model.calculate_resistance(&fluid, &cond).unwrap();
    // R = 128 μ L / (π D^4)
    let expected = 128.0f64 * fluid.viscosity * 1.0f64 / (std::f64::consts::PI * 0.01f64.powi(4));
    let rel_err = (r - expected).abs() / expected.abs();
    assert!(rel_err < 1e-12, "Hagen-Poiseuille resistance mismatch: r={} expected={} rel_err={}", r, expected, rel_err);
}

#[test]
fn test_darcy_weisbach_laminar_limit_friction_factor() {
    // Laminar regime: Re=1000 < 2000, expect f ≈ 64/Re
    let fluid = Fluid { density: 1000.0f64, viscosity: 1.0e-3f64 };
    let model = DarcyWeisbachModel::new(0.01f64, 1.0f64, 0.0f64);
    let mut cond = FlowConditions::new(0.0f64);
    cond.reynolds_number = Some(1000.0f64);

    let r = model.calculate_resistance(&fluid, &cond).unwrap();
    // Compute equivalent resistance using f=64/Re
    let re = 1000.0f64;
    let f = 64.0f64 / re;
    let area = std::f64::consts::PI * (0.01f64.powi(2)) / 4.0f64;
    let expected = f * 1.0f64 * fluid.density / (2.0f64 * area * 0.01f64.powi(2));
    let rel_err = (r - expected).abs() / expected.abs();
    assert!(rel_err < 1e-3, "Darcy-Weisbach laminar limit mismatch: r={} expected={} rel_err={}", r, expected, rel_err);
}
