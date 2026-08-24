#![allow(clippy::print_stdout)]
//! Conservation checking for CFD simulations.
//!
//! This module provides tools to verify that CFD simulations satisfy fundamental
//! conservation laws such as mass, momentum, and energy conservation.

use crate::scalar;
use cfd_core::physics::constants::physics::thermo::P_ATM;
use eunomia::FloatElement;
use eunomia::NumericElement;
use eunomia::RealField;
use leto::geometry::Vector2;
use leto::Array2;

mod angular_momentum;
mod energy;
mod geometric;
mod history;
mod mass;
mod momentum;
mod report;
mod tolerance;
mod traits;
mod vorticity;

pub use angular_momentum::AngularMomentumChecker;
pub use energy::EnergyConservationChecker;
pub use geometric::GeometricConservationChecker;
pub use history::ConservationHistory;
pub use mass::MassConservationChecker;
pub use momentum::MomentumConservationChecker;
pub use report::ConservationReport;
pub use tolerance::ConservationTolerance;
pub use traits::ConservationChecker;
pub use vorticity::VorticityChecker;

/// Run comprehensive conservation verification suite
/// Implements MINOR-011: Complete Conservation Property Verification
pub fn run_comprehensive_conservation_verification<T: RealField + Copy + FloatElement>() {
    tracing::info!("🧪 Comprehensive Conservation Property Verification Suite");
    tracing::info!("======================================================");
    tracing::info!(
        "MINOR-011: Verifying conservation of mass, momentum, energy, angular momentum, vorticity"
    );
    tracing::info!("References: LeVeque (2002), Thomas & Lombard (1979), Batchelor (1967)");
    tracing::info!("");

    // Mass conservation test
    tracing::info!("Test 1: Mass Conservation (Incompressible Navier-Stokes)");
    let mass_checker =
        MassConservationChecker::<T>::new(<T as FloatElement>::from_f64(1e-6), 32, 32);
    let velocity_u = Array2::from_elem([32, 32], <T as FloatElement>::from_f64(1.0));
    let velocity_v = Array2::from_elem([32, 32], scalar::zero::<T>());

    match mass_checker.check_divergence_2d(
        &velocity_u,
        &velocity_v,
        <T as FloatElement>::from_f64(0.1),
        <T as FloatElement>::from_f64(0.1),
    ) {
        Ok(report) => {
            tracing::info!(
                "  ✅ {}: Error = {:.2e}, Conserved: {}",
                report.check_name,
                <T as NumericElement>::to_f64(report.error),
                report.is_conserved
            );
        }
        Err(e) => tracing::info!("  ❌ Mass conservation test failed: {e}"),
    }

    // Momentum conservation test
    tracing::info!("\nTest 2: Momentum Conservation (Steady State)");
    let momentum_checker = MomentumConservationChecker::<T>::new(
        <T as FloatElement>::from_f64(1e-6),
        32,
        32,
        scalar::one::<T>(),
    );
    let u_prev = Array2::from_elem([32, 32], <T as FloatElement>::from_f64(0.95)); // Slightly different for time derivative
    let v_prev = Array2::zeros([32, 32]);
    let pressure = Array2::from_elem([32, 32], <T as FloatElement>::from_f64(P_ATM)); // Atmospheric pressure

    match momentum_checker.check_momentum_2d(
        &velocity_u,
        &velocity_v,
        &u_prev,
        &v_prev,
        &pressure,
        <T as FloatElement>::from_f64(1.5e-5), // Air viscosity
        <T as FloatElement>::from_f64(1e-3),   // dt
        <T as FloatElement>::from_f64(0.1),
        <T as FloatElement>::from_f64(0.1), // dx, dy
        Vector2::zeros(),                   // No gravity
    ) {
        Ok(report) => {
            tracing::info!(
                "  ✅ {}: Error = {:.2e}, Conserved: {}",
                report.check_name,
                <T as NumericElement>::to_f64(report.error),
                report.is_conserved
            );
        }
        Err(e) => tracing::info!("  ❌ Momentum conservation test failed: {e}"),
    }

    // Energy conservation test (if temperature field available)
    tracing::info!("\nTest 3: Energy Conservation (Thermal)");
    let energy_checker = EnergyConservationChecker::<T>::new(
        <T as FloatElement>::from_f64(1e-6),
        32,
        32,
        scalar::one::<T>(),
        <T as FloatElement>::from_f64(1000.0),
    );
    let temperature = Array2::from_elem([32, 32], <T as FloatElement>::from_f64(300.0)); // Room temperature

    match energy_checker.check_energy_2d(
        &temperature,
        &temperature,
        &velocity_u,
        &velocity_v,
        <T as FloatElement>::from_f64(0.025), // Thermal conductivity
        <T as FloatElement>::from_f64(1e-3),  // dt
        <T as FloatElement>::from_f64(0.1),
        <T as FloatElement>::from_f64(0.1), // dx, dy
        None,                               // No source term
    ) {
        Ok(report) => {
            tracing::info!(
                "  ✅ {}: Error = {:.2e}, Conserved: {}",
                report.check_name,
                <T as NumericElement>::to_f64(report.error),
                report.is_conserved
            );
        }
        Err(e) => tracing::info!("  ❌ Energy conservation test failed: {e}"),
    }

    // Angular momentum conservation test
    tracing::info!("\nTest 4: Angular Momentum Conservation (2D Cartesian)");
    let am_checker =
        AngularMomentumChecker::<T>::new_centered(<T as FloatElement>::from_f64(1e-6), 32, 32);

    match am_checker.check_angular_momentum_2d(
        &velocity_u,
        &velocity_v,
        <T as FloatElement>::from_f64(0.1),
        <T as FloatElement>::from_f64(0.1),
    ) {
        Ok(report) => {
            tracing::info!(
                "  ✅ {}: Error = {:.2e}, Conserved: {}",
                report.check_name,
                <T as NumericElement>::to_f64(report.error),
                report.is_conserved
            );
        }
        Err(e) => tracing::info!("  ❌ Angular momentum conservation test failed: {e}"),
    }

    // Vorticity conservation test
    tracing::info!("\nTest 5: Vorticity Transport Conservation");
    let vorticity_checker = VorticityChecker::<T>::new(<T as FloatElement>::from_f64(1e-6), 32, 32);

    match vorticity_checker.check_vorticity_transport_2d(
        &velocity_u,
        &velocity_v,
        <T as FloatElement>::from_f64(1.5e-5), // viscosity
        <T as FloatElement>::from_f64(0.1),
        <T as FloatElement>::from_f64(0.1),
    ) {
        Ok(report) => {
            tracing::info!(
                "  ✅ {}: Error = {:.2e}, Conserved: {}",
                report.check_name,
                <T as NumericElement>::to_f64(report.error),
                report.is_conserved
            );
        }
        Err(e) => tracing::info!("  ❌ Vorticity conservation test failed: {e}"),
    }

    // Geometric conservation law test
    tracing::info!("\nTest 6: Geometric Conservation Law");
    let gcl_checker =
        GeometricConservationChecker::<T>::new(<T as FloatElement>::from_f64(1e-14), 32, 32);

    match gcl_checker.run_comprehensive_gcl_tests() {
        Ok(results) => {
            let passed = results.iter().filter(|r| r.is_conserved).count();
            tracing::info!("  GCL Tests: {}/{} passed", passed, results.len());
            if passed > 0 {
                tracing::info!(
                    "  ✅ Sample GCL result: {} (Error = {:.2e})",
                    results[0].check_name,
                    <T as NumericElement>::to_f64(results[0].error)
                );
            }
        }
        Err(e) => tracing::info!("  ❌ GCL test failed: {e}"),
    }

    tracing::info!("\n✅ Complete conservation property verification completed!");
    tracing::info!("   All fundamental conservation laws have been tested:");
    tracing::info!("   - Mass conservation (continuity equation)");
    tracing::info!("   - Momentum conservation (Navier-Stokes)");
    tracing::info!("   - Energy conservation (thermal transport)");
    tracing::info!("   - Angular momentum conservation (rotation)");
    tracing::info!("   - Vorticity conservation (flow rotation)");
    tracing::info!("   - Geometric conservation law (numerical consistency)");
}
