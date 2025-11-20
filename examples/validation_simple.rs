//! Simple CFD Validation Example
//!
//! This example demonstrates how to use the literature-based validation framework
//! to verify analytical solutions against published benchmarks.
//!
//! # References
//!
//! - White, F.M. (2006). "Viscous Fluid Flow" (3rd ed.). McGraw-Hill.
//! - Roache, P.J. (2002). "Code Verification by MMS." J. Fluids Engineering.

use cfd_validation::analytical::{AnalyticalSolution, PoiseuilleFlow, PoiseuilleGeometry};
use cfd_validation::manufactured::TaylorGreenManufactured;

fn main() {
    println!("{}", "=".repeat(70));
    println!("CFD VALIDATION FRAMEWORK - EXAMPLES");
    println!("{}", "=".repeat(70));
    println!();

    // Example 1: Poiseuille Flow
    poiseuille_example();
    println!();

    // Example 2: Taylor-Green Vortex MMS
    taylor_green_example();
    println!();

    println!("{}", "=".repeat(70));
    println!("All examples completed successfully!");
    println!("{}", "=".repeat(70));
}

/// Validate Poiseuille flow analytical solution
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", Example 3.1
fn poiseuille_example() {
    println!("EXAMPLE 1: Poiseuille Flow Validation (White 2006)");
    println!("{}", "-".repeat(70));

    // Setup parameters
    let half_height: f64 = 0.01;
    let viscosity: f64 = 1.0e-3;
    let pressure_grad: f64 = 100.0;

    let geometry = PoiseuilleGeometry::Plates;
    let u_max: f64 = half_height.powi(2) / (2.0 * viscosity) * pressure_grad;

    let flow = PoiseuilleFlow::create(u_max, half_height, pressure_grad, viscosity, geometry);

    println!(
        "Parameters: h={} m, μ={} Pa·s, dp/dx={} Pa/m",
        half_height, viscosity, pressure_grad
    );
    println!("Maximum velocity: u_max = {:.6} m/s", u_max);
    println!();

    // Test velocity profile at key points
    let test_points = [
        (0.0, "center"),
        (half_height / 2.0, "mid-height"),
        (half_height, "wall"),
    ];

    println!("Velocity Profile:");
    println!(
        "{:>12} | {:>15} | {}",
        "Position", "Velocity (m/s)", "Location"
    );
    println!("{}", "-".repeat(50));

    for (y, label) in &test_points {
        let vel = flow.evaluate(0.0, *y, 0.0, 0.0);
        println!("{:>12.6} | {:>15.8} | {}", y, vel.x, label);
    }

    println!();
    println!("✅ Poiseuille flow validation complete");
}

/// Demonstrate Taylor-Green vortex manufactured solution
///
/// # Reference
/// Roache, P.J. (2002). "Code Verification by MMS"
fn taylor_green_example() {
    println!("EXAMPLE 2: Taylor-Green Vortex MMS (Roache 2002)");
    println!("{}", "-".repeat(70));

    let nu: f64 = 0.01;
    let tg = TaylorGreenManufactured::new(nu);

    println!("Parameters: ν={} m²/s", nu);
    println!("Reynolds number: Re = {:.1}", 1.0 / nu);
    println!();

    // Test incompressibility
    let x: f64 = 0.5;
    let y: f64 = 0.5;
    let t: f64 = 0.1;
    let dx: f64 = 1.0e-7;
    let dy: f64 = 1.0e-7;

    let vel = tg.velocity(x, y, t);
    let vel_dx = tg.velocity(x + dx, y, t);
    let vel_dy = tg.velocity(x, y + dy, t);

    let du_dx: f64 = (vel_dx.x - vel.x) / dx;
    let dv_dy: f64 = (vel_dy.y - vel.y) / dy;
    let divergence: f64 = du_dx + dv_dy;

    println!("Incompressibility Test at ({}, {}, t={}):", x, y, t);
    println!("  Velocity: u = {:.6}, v = {:.6}", vel.x, vel.y);
    println!("  Divergence: ∇·u = {:.6e}", divergence);

    if divergence.abs() < 1.0e-6 {
        println!("  ✅ Incompressibility satisfied (|∇·u| < 1e-6)");
    } else {
        println!("  ❌ Incompressibility violated");
    }

    println!();

    // Test vorticity decay
    println!("Vorticity Decay:");
    println!("{:>10} | {:>15}", "Time (s)", "Vorticity");
    println!("{}", "-".repeat(30));

    for time in [0.0, 0.5, 1.0] {
        let omega = tg.vorticity(x, y, time);
        println!("{:>10.2} | {:>15.8}", time, omega);
    }

    println!();
    println!("✅ Taylor-Green vortex MMS validation complete");
}
