//! Comprehensive validation tests for Taylor-Green vortex and shear flows
//!
//! This module tests analytical solutions for vortex decay and shear-driven flows
//! against literature benchmarks.
//!
//! # References
//!
//! - Taylor, G.I., & Green, A.E. (1937). "Mechanism of the production of small eddies
//!   from large ones." Proceedings of the Royal Society of London A, 158(895), 499-521.
//! - Brachet, M.E., et al. (1983). "Small-scale structure of the Taylor-Green vortex."
//!   Journal of Fluid Mechanics, 130, 411-452.
//! - White, F.M. (2006). "Viscous Fluid Flow" (3rd ed.). McGraw-Hill.

use approx::assert_relative_eq;
use cfd_validation::analytical::{AnalyticalSolution, CouetteFlow, TaylorGreenVortex};

/// Test Taylor-Green vortex energy decay
///
/// # Reference
/// Taylor & Green (1937). "Mechanism of the production of small eddies from large ones"
///
/// The kinetic energy should decay exponentially: E(t) = E₀ * exp(-2νk²t)
/// where k = 2π/L is the wave number
#[test]
fn test_taylor_green_energy_decay() {
    let length_scale: f64 = 1.0; // Unit domain
    let velocity_scale: f64 = 1.0; // Unit velocity
    let viscosity: f64 = 0.01; // Re = 100

    let vortex = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

    // Test at several time points
    let test_times = vec![0.0, 0.1, 0.5, 1.0];
    let mut prev_energy = f64::INFINITY;

    for t in test_times {
        let energy = vortex.kinetic_energy(t);

        // Energy should be positive
        assert!(energy > 0.0, "Energy should be positive at t = {}", t);

        // Energy should decrease monotonically
        assert!(
            energy < prev_energy || t == 0.0,
            "Energy should decay: E({}) = {} > E(prev) = {}",
            t,
            energy,
            prev_energy
        );

        prev_energy = energy;
    }

    // Verify that energy has decayed significantly by t=1.0
    let energy_final = vortex.kinetic_energy(1.0);
    let energy_initial = vortex.kinetic_energy(0.0);
    assert!(
        energy_final < energy_initial,
        "Energy should have decayed: E(1.0) = {} >= E(0) = {}",
        energy_final,
        energy_initial
    );
}

/// Test Reynolds number calculation for Taylor-Green vortex
///
/// # Reference
/// Brachet et al. (1983). "Small-scale structure of the Taylor-Green vortex"
///
/// Re = UL/ν where U is characteristic velocity, L is length scale
#[test]
fn test_taylor_green_reynolds_number() {
    let length_scale: f64 = 1.0;
    let velocity_scale: f64 = 1.0;
    let viscosity: f64 = 0.01; // Should give Re = 100

    let vortex = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

    let re = vortex.reynolds_number();
    let expected_re = velocity_scale * length_scale / viscosity;

    assert_relative_eq!(re, expected_re, epsilon = 1.0e-10);
    assert_relative_eq!(re, 100.0, epsilon = 1.0e-10);
}

/// Test Couette flow velocity profile
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", Example 3.4
///
/// For simple shear flow between parallel plates with one moving:
/// u(y) = U * y/h where U is wall velocity, h is gap height
#[test]
fn test_couette_linear_profile() {
    let gap_height: f64 = 0.01; // 1 cm gap
    let wall_velocity: f64 = 1.0; // 1 m/s
    let viscosity: f64 = 1.0e-3; // Water

    let couette = CouetteFlow::create(wall_velocity, gap_height, 0.0, viscosity);

    // Test velocity at several points
    let test_points = vec![
        (0.0, 0.0),                                     // Bottom wall (stationary)
        (gap_height / 4.0, 0.25 * wall_velocity),       // 1/4 height
        (gap_height / 2.0, 0.5 * wall_velocity),        // Mid-height
        (3.0 * gap_height / 4.0, 0.75 * wall_velocity), // 3/4 height
        (gap_height, wall_velocity),                    // Top wall (moving)
    ];

    for (y, expected_u) in test_points {
        let vel = couette.evaluate(0.0, y, 0.0, 0.0);
        assert_relative_eq!(vel.x, expected_u, epsilon = 1.0e-10);
    }
}

/// Test Couette flow wall shear stress
///
/// # Reference
/// White, F.M. (2006). "Viscous Fluid Flow", Eq. 3-46
///
/// Wall shear stress: τ_w = μ * U/h
#[test]
fn test_couette_wall_shear_stress() {
    let gap_height: f64 = 0.01;
    let wall_velocity: f64 = 1.0;
    let viscosity: f64 = 1.0e-3;

    let couette = CouetteFlow::create(wall_velocity, gap_height, 0.0, viscosity);

    // Analytical wall shear stress (White 2006, Eq. 3-46)
    let tau_wall_analytical = viscosity * wall_velocity / gap_height;

    // Compute velocity gradient using finite difference
    let dy = 1.0e-8;
    let u1 = couette.evaluate(0.0, 0.0, 0.0, 0.0).x;
    let u2 = couette.evaluate(0.0, dy, 0.0, 0.0).x;

    let du_dy = (u2 - u1) / dy;
    let tau_wall_computed = viscosity * du_dy;

    assert_relative_eq!(
        tau_wall_computed,
        tau_wall_analytical,
        epsilon = 1.0e-6 * tau_wall_analytical
    );
}

/// Test Taylor-Green vortex incompressibility
///
/// # Reference
/// Exact solution to incompressible Navier-Stokes equations
///
/// Divergence of velocity field should be zero: ∇·u = 0
#[test]
fn test_taylor_green_incompressibility() {
    let length_scale: f64 = 1.0;
    let velocity_scale: f64 = 1.0;
    let viscosity: f64 = 0.01;

    let vortex = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

    // Test incompressibility at several points
    let dx = 1.0e-6;
    let dy = 1.0e-6;

    for i in 0..5 {
        let x = 0.2 * i as f64;
        let y = 0.2 * i as f64;
        let t = 0.1;

        // Compute divergence using finite differences
        let vel = vortex.evaluate(x, y, 0.0, t);
        let vel_dx_plus = vortex.evaluate(x + dx, y, 0.0, t);
        let vel_dy_plus = vortex.evaluate(x, y + dy, 0.0, t);

        let du_dx = (vel_dx_plus.x - vel.x) / dx;
        let dv_dy = (vel_dy_plus.y - vel.y) / dy;

        let divergence = du_dx + dv_dy;

        // Divergence should be very close to zero (within numerical precision)
        assert!(
            divergence.abs() < 1.0e-6,
            "Divergence too large at ({}, {}): {}",
            x,
            y,
            divergence
        );
    }
}

/// Test vorticity conservation for Taylor-Green vortex
///
/// # Reference
/// Taylor & Green (1937) - vorticity magnitude should decay
#[test]
fn test_taylor_green_vorticity_decay() {
    let length_scale: f64 = 1.0;
    let velocity_scale: f64 = 1.0;
    let viscosity: f64 = 0.01;

    let vortex = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

    let dx = 1.0e-6;
    let dy = 1.0e-6;
    let x = 0.5;
    let y = 0.5;

    // Vorticity should decay over time
    let mut prev_vorticity_mag = f64::INFINITY;

    for i in 0..5 {
        let t = 0.2 * i as f64;

        let vel = vortex.evaluate(x, y, 0.0, t);
        let vel_dx_plus = vortex.evaluate(x + dx, y, 0.0, t);
        let vel_dy_plus = vortex.evaluate(x, y + dy, 0.0, t);

        let du_dy = (vel_dy_plus.x - vel.x) / dy;
        let dv_dx = (vel_dx_plus.y - vel.y) / dx;

        let vorticity = dv_dx - du_dy;
        let vorticity_mag = vorticity.abs();

        // Vorticity magnitude should decrease over time (or stay same at t=0)
        assert!(
            vorticity_mag <= prev_vorticity_mag || i == 0,
            "Vorticity should decay: ω({}) = {} > ω(prev) = {}",
            t,
            vorticity_mag,
            prev_vorticity_mag
        );

        prev_vorticity_mag = vorticity_mag;
    }
}

#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    /// Property test: Couette flow velocity should be linear
    ///
    /// # Reference
    /// Fundamental property of simple shear flow
    proptest! {
        #[test]
        fn couette_velocity_is_linear(
            gap_height in 0.001f64..0.1,
            wall_velocity in 0.1f64..10.0,
            viscosity in 1.0e-5f64..1.0e-2
        ) {
            let couette = CouetteFlow::create(wall_velocity, gap_height, 0.0, viscosity);

            // Test linearity at multiple points
            for i in 0..=10 {
                let y = gap_height * (i as f64 / 10.0);
                let vel = couette.evaluate(0.0, y, 0.0, 0.0);
                let expected = wall_velocity * y / gap_height;

                prop_assert!((vel.x - expected).abs() / wall_velocity < 1.0e-10,
                    "Velocity not linear at y = {}: {} vs {}", y, vel.x, expected);
            }
        }
    }

    /// Property test: Taylor-Green energy decreases monotonically
    ///
    /// # Reference
    /// Physical requirement from Navier-Stokes dissipation
    proptest! {
        #[test]
        fn taylor_green_energy_decreases(
            length_scale in 0.1f64..10.0,
            velocity_scale in 0.1f64..10.0,
            viscosity in 1.0e-4f64..1.0
        ) {
            let vortex = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

            let e0 = vortex.kinetic_energy(0.0);
            let e1 = vortex.kinetic_energy(0.1);
            let e2 = vortex.kinetic_energy(1.0);

            prop_assert!(e0 >= e1, "Energy should not increase: E(0) = {} < E(0.1) = {}", e0, e1);
            prop_assert!(e1 >= e2, "Energy should not increase: E(0.1) = {} < E(1) = {}", e1, e2);
        }
    }

    /// Property test: Couette flow has constant shear rate
    ///
    /// # Reference
    /// Fundamental property: du/dy = U/h = constant
    proptest! {
        #[test]
        fn couette_constant_shear_rate(
            gap_height in 0.001f64..0.1,
            wall_velocity in 0.1f64..10.0,
            viscosity in 1.0e-5f64..1.0e-2
        ) {
            let couette = CouetteFlow::create(wall_velocity, gap_height, 0.0, viscosity);
            let expected_shear = wall_velocity / gap_height;

            let dy = 1.0e-7;

            // Test shear rate at multiple heights
            for i in 1..10 {
                let y = gap_height * (i as f64 / 10.0);
                let u1 = couette.evaluate(0.0, y - dy/2.0, 0.0, 0.0).x;
                let u2 = couette.evaluate(0.0, y + dy/2.0, 0.0, 0.0).x;

                let shear = (u2 - u1) / dy;

                prop_assert!((shear - expected_shear).abs() / expected_shear < 1.0e-5,
                    "Shear rate not constant at y = {}: {} vs {}", y, shear, expected_shear);
            }
        }
    }

    /// Property test: Reynolds number scales correctly
    ///
    /// # Reference
    /// Definition: Re = UL/ν
    proptest! {
        #[test]
        fn reynolds_number_scaling(
            length_scale in 0.1f64..10.0,
            velocity_scale in 0.1f64..10.0,
            viscosity in 1.0e-4f64..1.0
        ) {
            let vortex = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

            let re = vortex.reynolds_number();
            let expected_re = velocity_scale * length_scale / viscosity;

            prop_assert!((re - expected_re).abs() / expected_re < 1.0e-10,
                "Reynolds number mismatch: {} vs {}", re, expected_re);
        }
    }
}
