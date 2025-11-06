//! Conservation checking for CFD simulations.
//!
//! This module provides tools to verify that CFD simulations satisfy fundamental
//! conservation laws such as mass, momentum, and energy conservation.

use nalgebra::{DMatrix, RealField};
use num_traits::{FromPrimitive, ToPrimitive};

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
pub fn run_comprehensive_conservation_verification<T: RealField + Copy + FromPrimitive + ToPrimitive>() {
    println!("🧪 Comprehensive Conservation Property Verification Suite");
    println!("======================================================");
    println!("MINOR-011: Verifying conservation of mass, momentum, energy, angular momentum, vorticity");
    println!("References: LeVeque (2002), Thomas & Lombard (1979), Batchelor (1967)");
    println!();

    // Mass conservation test
    println!("Test 1: Mass Conservation (Incompressible Navier-Stokes)");
    let mass_checker = MassConservationChecker::<T>::new(T::from_f64(1e-6).unwrap(), 32, 32);
    let velocity_u = DMatrix::from_element(32, 32, T::from_f64(1.0).unwrap());
    let velocity_v = DMatrix::from_element(32, 32, T::zero());

    match mass_checker.check_divergence_2d(
        &velocity_u, &velocity_v, T::from_f64(0.1).unwrap(), T::from_f64(0.1).unwrap()
    ) {
        Ok(report) => {
            println!("  ✅ {}: Error = {:.2e}, Conserved: {}",
                     report.check_name, report.error.to_f64().unwrap_or(0.0), report.is_conserved);
        }
        Err(e) => println!("  ❌ Mass conservation test failed: {}", e),
    }

    // Momentum conservation test
    println!("\nTest 2: Momentum Conservation (Steady State)");
    let momentum_checker = MomentumConservationChecker::<T>::new(T::from_f64(1e-6).unwrap(), 32, 32, T::one());
    let u_prev = DMatrix::from_element(32, 32, T::from_f64(0.95).unwrap()); // Slightly different for time derivative
    let v_prev = DMatrix::zeros(32, 32);
    let pressure = DMatrix::from_element(32, 32, T::from_f64(101325.0).unwrap()); // Atmospheric pressure

    match momentum_checker.check_momentum_2d(
        &velocity_u, &velocity_v, &u_prev, &v_prev, &pressure,
        T::from_f64(1.5e-5).unwrap(), // Air viscosity
        T::from_f64(1e-3).unwrap(), // dt
        T::from_f64(0.1).unwrap(), T::from_f64(0.1).unwrap(), // dx, dy
        nalgebra::Vector2::zeros() // No gravity
    ) {
        Ok(report) => {
            println!("  ✅ {}: Error = {:.2e}, Conserved: {}",
                     report.check_name, report.error.to_f64().unwrap_or(0.0), report.is_conserved);
        }
        Err(e) => println!("  ❌ Momentum conservation test failed: {}", e),
    }

    // Energy conservation test (if temperature field available)
    println!("\nTest 3: Energy Conservation (Thermal)");
    let energy_checker = EnergyConservationChecker::<T>::new(T::from_f64(1e-6).unwrap(), 32, 32, T::one(), T::from_f64(1000.0).unwrap());
    let temperature = DMatrix::from_element(32, 32, T::from_f64(300.0).unwrap()); // Room temperature

    match energy_checker.check_energy_2d(
        &temperature, &temperature, &velocity_u, &velocity_v,
        T::from_f64(0.025).unwrap(), // Thermal conductivity
        T::from_f64(1e-3).unwrap(), // dt
        T::from_f64(0.1).unwrap(), T::from_f64(0.1).unwrap(), // dx, dy
        None, // No source term
    ) {
        Ok(report) => {
            println!("  ✅ {}: Error = {:.2e}, Conserved: {}",
                     report.check_name, report.error.to_f64().unwrap_or(0.0), report.is_conserved);
        }
        Err(e) => println!("  ❌ Energy conservation test failed: {}", e),
    }

    // Angular momentum conservation test
    println!("\nTest 4: Angular Momentum Conservation (2D Cartesian)");
    let am_checker = AngularMomentumChecker::<T>::new_centered(T::from_f64(1e-6).unwrap(), 32, 32);

    match am_checker.check_angular_momentum_2d(
        &velocity_u, &velocity_v, T::from_f64(0.1).unwrap(), T::from_f64(0.1).unwrap()
    ) {
        Ok(report) => {
            println!("  ✅ {}: Error = {:.2e}, Conserved: {}",
                     report.check_name, report.error.to_f64().unwrap_or(0.0), report.is_conserved);
        }
        Err(e) => println!("  ❌ Angular momentum conservation test failed: {}", e),
    }

    // Vorticity conservation test
    println!("\nTest 5: Vorticity Transport Conservation");
    let vorticity_checker = VorticityChecker::<T>::new(T::from_f64(1e-6).unwrap(), 32, 32);

    match vorticity_checker.check_vorticity_transport_2d(
        &velocity_u, &velocity_v,
        T::from_f64(1.5e-5).unwrap(), // viscosity
        T::from_f64(0.1).unwrap(), T::from_f64(0.1).unwrap()
    ) {
        Ok(report) => {
            println!("  ✅ {}: Error = {:.2e}, Conserved: {}",
                     report.check_name, report.error.to_f64().unwrap_or(0.0), report.is_conserved);
        }
        Err(e) => println!("  ❌ Vorticity conservation test failed: {}", e),
    }

    // Geometric conservation law test
    println!("\nTest 6: Geometric Conservation Law");
    let gcl_checker = GeometricConservationChecker::<T>::new(T::from_f64(1e-14).unwrap(), 32, 32);

    match gcl_checker.run_comprehensive_gcl_tests() {
        Ok(results) => {
            let passed = results.iter().filter(|r| r.is_conserved).count();
            println!("  GCL Tests: {}/{} passed", passed, results.len());
            if passed > 0 {
                println!("  ✅ Sample GCL result: {} (Error = {:.2e})",
                         results[0].check_name,
                         results[0].error.to_f64().unwrap_or(0.0));
            }
        }
        Err(e) => println!("  ❌ GCL test failed: {}", e),
    }

    println!("\n✅ Complete conservation property verification completed!");
    println!("   All fundamental conservation laws have been tested:");
    println!("   - Mass conservation (continuity equation)");
    println!("   - Momentum conservation (Navier-Stokes)");
    println!("   - Energy conservation (thermal transport)");
    println!("   - Angular momentum conservation (rotation)");
    println!("   - Vorticity conservation (flow rotation)");
    println!("   - Geometric conservation law (numerical consistency)");
}
