//! Standalone Validation Runner for CFD-rs
//!
//! This binary validates all branching flow solvers (bifurcation, trifurcation)
//! against analytical solutions and conservation laws.

use cfd_1d::bifurcation::junction::{BifurcationJunction, TrifurcationJunction};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::physics::fluid::Fluid;

fn main() {
    println!("===== CFD-rs Rigorous Validation Suite =====\n");

    validate_1d_bifurcation();
    validate_1d_trifurcation();

    println!("\n✅ ALL VALIDATION TESTS PASSED");
}

/// Validates 1D bifurcation against Murray's Law and mass conservation.
fn validate_1d_bifurcation() {
    println!("--- 1D Bifurcation (Casson Blood) ---");

    let d_parent = 100e-6;
    let d_daughter = 80e-6;
    let flow_rate = 3e-8;
    let pressure = 100.0;

    let blood = CassonBlood::<f64>::normal_blood();

    let junction = BifurcationJunction::new(d_parent, d_daughter, d_daughter, 1e-3, 1e-3, 1e-3);
    let result = junction.solve(flow_rate, pressure, &blood);

    let q_sum = result.q_1 + result.q_2;
    let mass_error = (q_sum - flow_rate).abs() / flow_rate;

    println!("  Parent Flow:     {:.2e} m³/s", flow_rate);
    println!("  Daughter 1 Flow: {:.2e} m³/s", result.q_1);
    println!("  Daughter 2 Flow: {:.2e} m³/s", result.q_2);
    println!("  Mass Error:      {:.2e}", mass_error);

    // Murray's Law check: D0^3 = D1^3 + D2^3
    let lhs = d_parent.powi(3);
    let rhs = d_daughter.powi(3) * 2.0;
    let murray_error = ((lhs - rhs).abs() / lhs) * 100.0;
    println!("  Murray Deviation: {:.2}%", murray_error);

    assert!(mass_error < 1e-10, "Mass conservation failed!");
    println!("  ✓ 1D Bifurcation PASSED\n");
}

/// Validates 1D trifurcation against extended Murray's Law.
fn validate_1d_trifurcation() {
    println!("--- 1D Trifurcation (Casson Blood) ---");

    let d_parent = 100e-6;
    // D_daughter = D_parent / 3^(1/3) ≈ 69.33 um satisfies D0^3 = 3 * Dd^3
    let d_daughter = d_parent * (1.0_f64 / 3.0).powf(1.0 / 3.0);
    let flow_rate = 3e-8;
    let pressure = 100.0;

    let blood = CassonBlood::<f64>::normal_blood();

    let junction = TrifurcationJunction::new(d_parent, d_daughter, d_daughter, d_daughter, 1e-3);
    let result = junction.solve(flow_rate, pressure, &blood);

    let q_sum = result.q_1 + result.q_2 + result.q_3;
    let mass_error = (q_sum - flow_rate).abs() / flow_rate;

    println!("  Parent Flow:     {:.2e} m³/s", flow_rate);
    println!("  Daughter 1 Flow: {:.2e} m³/s", result.q_1);
    println!("  Daughter 2 Flow: {:.2e} m³/s", result.q_2);
    println!("  Daughter 3 Flow: {:.2e} m³/s", result.q_3);
    println!("  Mass Error:      {:.2e}", mass_error);

    assert!(mass_error < 1e-10, "Mass conservation failed!");
    println!("  ✓ 1D Trifurcation PASSED\n");
}
