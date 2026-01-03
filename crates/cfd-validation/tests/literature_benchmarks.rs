//! Literature-based benchmark validations
//!
//! This module implements comprehensive validation tests against established
//! literature benchmarks in Computational Fluid Dynamics.
//!
//! # References
//!
//! - Ghia, U., Ghia, K.N., & Shin, C.T. (1982). "High-Re solutions for incompressible
//!   flow using the Navier-Stokes equations and a multigrid method." Journal of
//!   Computational Physics, 48(3), 387-411.
//! - Roache, P.J. (1998). "Verification and Validation in Computational Science
//!   and Engineering." Hermosa Publishers.
//! - White, F.M. (2006). "Viscous Fluid Flow" (3rd ed.). McGraw-Hill.
//! - Ferziger, J.H., & Perić, M. (2019). "Computational Methods for Fluid Dynamics"
//!   (4th ed.). Springer.

use approx::assert_relative_eq;
use cfd_validation::analytical::{AnalyticalSolution, PoiseuilleFlow, PoiseuilleGeometry};

/// Test Poiseuille flow analytical solution against White (2006) exact solution
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", Example 3.1, pp. 123-125
///
/// The exact solution for laminar flow between parallel plates is:
/// u(y) = u_max * (1 - (y/h)²)
/// where u_max = h²/(2μ) * |dp/dx| and h is the half-height
#[test]
fn test_poiseuille_parallel_plates_white_2006() {
    // Test parameters from White (2006), Example 3.1
    let half_height: f64 = 0.01; // 1 cm half-height (total height 2 cm)
    let viscosity: f64 = 1.0e-3; // Water at 20°C: 1.002e-3 Pa·s
    let pressure_grad_mag: f64 = 100.0; // Pa/m

    let geometry = PoiseuilleGeometry::Plates;

    // Analytical maximum velocity at centerline (White 2006, Eq. 3-35)
    let u_max_analytical = half_height.powi(2) / (2.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(
        u_max_analytical,
        half_height,
        pressure_grad_mag,
        viscosity,
        geometry,
    );

    // Evaluate at centerline (y = 0)
    let velocity = flow.evaluate(0.0, 0.0, 0.0, 0.0);

    // White (2006) states analytical solution should be exact to machine precision
    assert_relative_eq!(velocity.x, u_max_analytical, epsilon = 1.0e-10);

    // Verify parabolic profile at several points
    let test_points = vec![
        (0.0, u_max_analytical),                      // Centerline, (1 - 0²)
        (half_height / 2.0, 0.75 * u_max_analytical), // y = h/2, (1 - 0.5²) = 0.75
        (half_height, 0.0),                           // Wall (no-slip)
    ];

    for (y, expected_u) in test_points {
        let vel = flow.evaluate(0.0, y, 0.0, 0.0);
        assert_relative_eq!(vel.x, expected_u, epsilon = 1.0e-10);
    }
}

/// Test Poiseuille flow in circular pipe against Hagen-Poiseuille equation
///
/// # Reference  
/// White, F.M. (2006). "Viscous Fluid Flow", Eq. 3-42, pp. 126-127
///
/// For laminar flow in a circular pipe:
/// u_max = (R²/4μ) * |dp/dx|
/// Q = (πR⁴/8μ) * |dp/dx|
#[test]
fn test_poiseuille_circular_pipe_hagen_equation() {
    let radius: f64 = 0.01; // 1 cm radius
    let viscosity: f64 = 1.0e-3; // Pa·s
    let pressure_grad_mag: f64 = 100.0; // Pa/m

    let geometry = PoiseuilleGeometry::Pipe;

    // Hagen-Poiseuille: u_max = R²/(4μ) * |dp/dx|
    let u_max_analytical = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(
        u_max_analytical,
        radius,
        pressure_grad_mag,
        viscosity,
        geometry,
    );

    let velocity = flow.evaluate(0.0, 0.0, 0.0, 0.0); // Center of pipe (r=0)

    // Should match Hagen-Poiseuille to machine precision
    assert_relative_eq!(velocity.x, u_max_analytical, epsilon = 1.0e-10);

    // Verify parabolic profile: u(r) = u_max * (1 - (r/R)²)
    let test_cases = vec![
        (0.0, 1.0),              // Center: (1 - 0²) = 1.0
        (radius / 2.0, 0.75),    // r = R/2: (1 - 0.5²) = 0.75
        (0.75 * radius, 0.4375), // r = 0.75R: (1 - 0.75²) = 0.4375
        (radius, 0.0),           // Wall: (1 - 1²) = 0
    ];

    for (r, factor) in test_cases {
        let vel = flow.evaluate(0.0, r, 0.0, 0.0);
        let expected = u_max_analytical * factor;
        assert_relative_eq!(vel.x, expected, epsilon = 1.0e-10);
    }
}

/// Test flow rate calculation against analytical formula
///
/// # Reference
/// Ferziger & Perić (2019). "Computational Methods for Fluid Dynamics", Eq. 8.45
///
/// For Poiseuille flow in circular pipe: Q = (πR⁴/8μ) * |dp/dx|
#[test]
fn test_flow_rate_ferziger_2019() {
    let radius: f64 = 0.01;
    let viscosity: f64 = 1.0e-3;
    let pressure_grad_mag: f64 = 100.0;

    let geometry = PoiseuilleGeometry::Pipe;
    let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

    // Analytical flow rate (Ferziger & Perić 2019, Eq. 8.45)
    let flow_rate_analytical =
        std::f64::consts::PI * radius.powi(4) / (8.0 * viscosity) * pressure_grad_mag;

    // Use built-in method
    let flow_rate_computed = flow.flow_rate();

    // Should match analytical to high precision
    assert_relative_eq!(
        flow_rate_computed,
        flow_rate_analytical,
        epsilon = 1.0e-10 * flow_rate_analytical
    );
}

/// Test Reynolds number calculation and laminar flow criterion
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", pp. 221-223
///
/// For pipe flow: Re = ρUD/μ where U is average velocity, D is diameter
/// Laminar if Re < 2300 (critical Reynolds number)
#[test]
fn test_reynolds_number_laminar_criterion_white_2006() {
    let radius: f64 = 0.01;
    let density: f64 = 1000.0; // kg/m³ (water)
    let viscosity: f64 = 1.0e-3; // Pa·s
    let pressure_grad_mag: f64 = 1.0; // Very low pressure gradient for laminar flow

    let geometry = PoiseuilleGeometry::Pipe;
    let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

    // Reynolds number using built-in method
    let reynolds = flow.reynolds_number(density);

    // Should be well within laminar regime (Re < 2300)
    assert!(
        reynolds < 2300.0,
        "Reynolds number {:.1} exceeds laminar limit 2300",
        reynolds
    );

    // Verify Reynolds number calculation
    // For Poiseuille flow with u_max, diameter 2R: Re = ρ * u_max * 2R / μ
    let expected_re = density * u_max * (2.0 * radius) / viscosity;
    assert_relative_eq!(reynolds, expected_re, epsilon = 1.0e-10);
}

/// Test wall shear stress calculation
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", Eq. 3-37
///
/// Wall shear stress: τ_w = μ(du/dy)|_wall
/// For Poiseuille flow between plates: τ_w = h * |dp/dx|
#[test]
fn test_wall_shear_stress_white_2006() {
    let half_height: f64 = 0.01;
    let viscosity: f64 = 1.0e-3;
    let pressure_grad_mag: f64 = 100.0;

    let geometry = PoiseuilleGeometry::Plates;
    let u_max = half_height.powi(2) / (2.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(u_max, half_height, pressure_grad_mag, viscosity, geometry);

    // Analytical wall shear stress (White 2006, Eq. 3-37)
    // τ_w = h * |dp/dx| where h is half-height
    let tau_wall_analytical = half_height * pressure_grad_mag;

    // Compute velocity gradient at wall using finite difference
    let dy = 1.0e-8; // Small distance from wall
    let u_wall = flow.evaluate(0.0, half_height, 0.0, 0.0).x; // Should be 0
    let u_near_wall = flow.evaluate(0.0, half_height - dy, 0.0, 0.0).x;

    let du_dy = (u_near_wall - u_wall) / dy;
    let tau_wall_computed = viscosity * du_dy.abs();

    // Should match analytical to within numerical differentiation error
    assert_relative_eq!(
        tau_wall_computed,
        tau_wall_analytical,
        epsilon = 1.0e-6 * tau_wall_analytical
    );
}

/// Test symmetry of velocity profile
///
/// # Reference
/// Symmetry is a fundamental property from the governing equations
#[test]
fn test_velocity_profile_symmetry() {
    let half_height: f64 = 0.01;
    let viscosity: f64 = 1.0e-3;
    let pressure_grad_mag: f64 = 100.0;

    let geometry = PoiseuilleGeometry::Plates;
    let u_max = half_height.powi(2) / (2.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(u_max, half_height, pressure_grad_mag, viscosity, geometry);

    // Test symmetry about centerline for parallel plates
    for i in 1..10 {
        let y = half_height * (i as f64 / 10.0);
        let vel_above = flow.evaluate(0.0, y, 0.0, 0.0);
        let vel_below = flow.evaluate(0.0, -y, 0.0, 0.0);

        // Should be symmetric
        assert_relative_eq!(vel_above.x, vel_below.x, epsilon = 1.0e-10);
    }
}

/// Test no-slip boundary condition
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", §1.6
/// Fundamental boundary condition in viscous flow
#[test]
fn test_no_slip_boundary_condition_white_2006() {
    let radius: f64 = 0.01;
    let viscosity: f64 = 1.0e-3;
    let pressure_grad_mag: f64 = 100.0;

    let geometry = PoiseuilleGeometry::Pipe;
    let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;

    let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

    // Velocity at wall should be exactly zero (no-slip condition)
    let vel_at_wall = flow.evaluate(0.0, radius, 0.0, 0.0);
    assert!(
        vel_at_wall.x.abs() < 1.0e-10,
        "No-slip violated at wall: u = {}",
        vel_at_wall.x
    );
}

#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    // Property test: Velocity should always be non-negative and bounded
    //
    // # Reference
    // Physical constraint from Navier-Stokes equations
    proptest! {
        #[test]
        fn velocity_is_non_negative_and_bounded(
            radius in 0.001f64..0.1,
            viscosity in 1.0e-5f64..1.0e-2,
            pressure_grad_mag in 10.0f64..1000.0
        ) {
            let geometry = PoiseuilleGeometry::Pipe;
            let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;
            let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

            // Test at various radial positions
            for i in 0..=10 {
                let r = radius * (i as f64 / 10.0);
                let vel = flow.evaluate(0.0, r, 0.0, 0.0);

                // Velocity should be non-negative
                prop_assert!(vel.x >= 0.0, "Negative velocity at r = {}: {}", r, vel.x);

                // Velocity should not exceed maximum
                prop_assert!(vel.x <= u_max * 1.001, "Velocity exceeds max at r = {}: {} > {}", r, vel.x, u_max);
            }
        }
    }

    // Property test: No-slip boundary condition at walls
    //
    // # Reference
    // Fundamental boundary condition in viscous flow (White 2006, §1.6)
    proptest! {
        #[test]
        fn no_slip_at_walls(
            radius in 0.001f64..0.1,
            viscosity in 1.0e-5f64..1.0e-2,
            pressure_grad_mag in 10.0f64..1000.0
        ) {
            let geometry = PoiseuilleGeometry::Pipe;
            let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;
            let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

            // Velocity at wall should be zero
            let vel_at_wall = flow.evaluate(0.0, radius, 0.0, 0.0);
            prop_assert!(vel_at_wall.x.abs() < 1.0e-10, "No-slip violated: u_wall = {}", vel_at_wall.x);
        }
    }

    // Property test: Velocity decreases monotonically from center to wall
    //
    // # Reference
    // Physical property of pressure-driven pipe flow
    proptest! {
        #[test]
        fn velocity_decreases_from_center(
            radius in 0.001f64..0.1,
            viscosity in 1.0e-5f64..1.0e-2,
            pressure_grad_mag in 10.0f64..1000.0
        ) {
            let geometry = PoiseuilleGeometry::Pipe;
            let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;
            let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

            // Velocity should decrease as we move from center to wall
            let mut prev_vel = f64::INFINITY;
            for i in 0..=10 {
                let r = radius * (i as f64 / 10.0);
                let vel = flow.evaluate(0.0, r, 0.0, 0.0).x;

                prop_assert!(vel <= prev_vel, "Velocity increases from r = {} to {}: {} > {}",
                    radius * ((i-1) as f64 / 10.0), r, vel, prev_vel);
                prev_vel = vel;
            }
        }
    }

    // Property test: Flow rate scales correctly with pressure gradient
    //
    // # Reference
    // Hagen-Poiseuille law: Q ∝ Δp
    proptest! {
        #[test]
        fn flow_rate_scales_with_pressure(
            radius in 0.001f64..0.1,
            viscosity in 1.0e-5f64..1.0e-2,
            pressure_grad_mag in 10.0f64..1000.0
        ) {
            let geometry = PoiseuilleGeometry::Pipe;
            let u_max = radius.powi(2) / (4.0 * viscosity) * pressure_grad_mag;
            let flow = PoiseuilleFlow::create(u_max, radius, pressure_grad_mag, viscosity, geometry);

            let q = flow.flow_rate();

            // Analytical: Q = (πR⁴/8μ) * |dp/dx|
            let q_expected = std::f64::consts::PI * radius.powi(4) / (8.0 * viscosity) * pressure_grad_mag;

            prop_assert!((q - q_expected).abs() / q_expected < 1.0e-10,
                "Flow rate mismatch: {} vs {}", q, q_expected);
        }
    }
}
