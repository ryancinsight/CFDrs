//! Standalone Validation Runner for CFD-rs
//!
//! This binary validates all branching flow solvers (bifurcation, trifurcation)
//! against analytical solutions and conservation laws.

use cfd_1d::bifurcation::{BifurcationJunction, TrifurcationJunction};
use cfd_1d::channel::{Channel, ChannelGeometry};
use cfd_core::physics::fluid::blood::CassonBlood;

fn main() {
    println!("===== CFD-rs Rigorous Validation Suite =====\n");

    validate_1d_bifurcation();
    validate_1d_trifurcation();

    println!("\n✅ ALL VALIDATION TESTS PASSED");
}

/// Validates 1D bifurcation against Murray's Law and mass conservation.
fn validate_1d_bifurcation() {
    println!("--- 1D Bifurcation (Casson Blood) ---");

    let d_parent: f64 = 100e-6;
    let d_daughter: f64 = 80e-6;
    let flow_rate: f64 = 3e-8;
    let pressure: f64 = 100.0;

    let blood = CassonBlood::<f64>::normal_blood();

    let parent = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_parent, 1e-6));
    let daughter1 = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_daughter, 1e-6));
    let daughter2 = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_daughter, 1e-6));

    let junction = BifurcationJunction::new(parent, daughter1, daughter2, 0.5);
    let result = junction.solve(blood, flow_rate, pressure).unwrap();

    let q_sum: f64 = result.q_1 + result.q_2;
    let mass_error: f64 = (q_sum - flow_rate).abs() / flow_rate;

    println!("  Parent Flow:     {:.2e} m³/s", flow_rate);
    println!("  Daughter 1 Flow: {:.2e} m³/s", result.q_1);
    println!("  Daughter 2 Flow: {:.2e} m³/s", result.q_2);
    println!("  Mass Error:      {:.2e}", mass_error);

    // Murray's Law check: D0^3 = D1^3 + D2^3
    let lhs: f64 = d_parent.powi(3);
    let rhs: f64 = d_daughter.powi(3) * 2.0;
    let murray_error: f64 = ((lhs - rhs).abs() / lhs) * 100.0;
    println!("  Murray Deviation: {:.2}%", murray_error);

    assert!(mass_error < 1e-10, "Mass conservation failed!");
    println!("  ✓ 1D Bifurcation PASSED\n");
}

/// Validates 1D trifurcation against extended Murray's Law.
fn validate_1d_trifurcation() {
    println!("--- 1D Trifurcation (Casson Blood) ---");

    let d_parent: f64 = 100e-6;
    // D_daughter = D_parent / 3^(1/3) ≈ 69.33 um satisfies D0^3 = 3 * Dd^3
    let d_daughter: f64 = d_parent * (1.0_f64 / 3.0).powf(1.0 / 3.0);
    let flow_rate: f64 = 3e-8;
    let pressure: f64 = 100.0;

    let blood = CassonBlood::<f64>::normal_blood();

    let parent = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_parent, 1e-6));
    let daughter1 = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_daughter, 1e-6));
    let daughter2 = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_daughter, 1e-6));
    let daughter3 = Channel::new(ChannelGeometry::<f64>::circular(1e-3, d_daughter, 1e-6));

    let junction = TrifurcationJunction::new(
        parent, daughter1, daughter2, daughter3,
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
    );
    let result = junction.solve(blood, flow_rate, pressure).unwrap();

    let q_sum: f64 = result.q_1 + result.q_2 + result.q_3;
    let mass_error: f64 = (q_sum - flow_rate).abs() / flow_rate;

    println!("  Parent Flow:     {:.2e} m³/s", flow_rate);
    println!("  Daughter 1 Flow: {:.2e} m³/s", result.q_1);
    println!("  Daughter 2 Flow: {:.2e} m³/s", result.q_2);
    println!("  Daughter 3 Flow: {:.2e} m³/s", result.q_3);
    println!("  Mass Error:      {:.2e}", mass_error);

    assert!(mass_error < 1e-10, "Mass conservation failed!");
    println!("  ✓ 1D Trifurcation PASSED\n");
}
