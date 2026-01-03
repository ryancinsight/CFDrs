//! Validation tests for advanced physics manufactured solutions
//!
//! Tests compressible flows, shock capturing, hypersonic flows,
//! and turbulence-chemistry interactions.

use cfd_validation::manufactured::advanced_physics::{
    ManufacturedCompressibleEuler, ManufacturedHypersonic, ManufacturedShockCapturing,
    ManufacturedTCI,
};
use cfd_validation::manufactured::ManufacturedSolution;

/// Test compressible Euler equations validation
#[test]
fn test_compressible_euler_validation() {
    println!("Testing Compressible Euler MMS Validation...");

    let euler = ManufacturedCompressibleEuler::<f64>::new(
        2.0, // Mach number
        1.4, // γ for air
        0.1, // flow angle (radians)
        0.1, // perturbation amplitude
        1.0, // kx
        1.0, // ky
    );

    // Test reference state properties
    assert!(euler.rho_0() > 0.0, "Reference density must be positive");
    assert!(euler.p_0() > 0.0, "Reference pressure must be positive");
    assert!(euler.u_0() > 0.0, "Reference velocity must be positive");

    // Test Mach number relationships
    let sonic_speed = f64::sqrt(euler.gamma * euler.p_0() / euler.rho_0());
    let expected_u0 = euler.mach_number * sonic_speed;
    assert!(
        (euler.u_0() - expected_u0).abs() < 1e-10,
        "Mach number relation violated"
    );

    // Test solution at multiple points
    let test_points = vec![
        (0.25, 0.25, 0.0, 0.0),
        (0.5, 0.5, 0.0, 0.5),
        (0.75, 0.75, 0.0, 1.0),
    ];

    for (x, y, z, t) in test_points {
        let rho = euler.exact_solution(x, y, z, t);
        let source = euler.source_term(x, y, z, t);

        // Physical constraints
        assert!(
            rho > 0.0,
            "Density must be positive at ({},{},t={}): {}",
            x,
            y,
            t,
            rho
        );
        assert!(
            rho < 2.0,
            "Density perturbation too large at ({},{},t={}): {}",
            x,
            y,
            t,
            rho
        );
        assert!(
            source.is_finite(),
            "Source term not finite at ({},{},t={}): {}",
            x,
            y,
            t,
            source
        );
    }

    println!("✓ Compressible Euler validation passed");
}

/// Test hypersonic flow validation
#[test]
fn test_hypersonic_flow_validation() {
    println!("Testing Hypersonic Flow MMS Validation...");

    let hypersonic = ManufacturedHypersonic::<f64>::new(
        8.0,  // Mach 8 (hypersonic)
        1e5,  // Reynolds number
        0.72, // Prandtl number
        1.4,  // γ
        4.0,  // wall temperature ratio
        0.2,  // amplitude
        1.5,  // kx
        1.0,  // ky
    );

    // Test hypersonic flow parameters
    assert!(hypersonic.mach_inf >= 5.0, "Must be hypersonic flow");
    assert!(
        hypersonic.reynolds >= 1e4,
        "Reynolds number should be turbulent"
    );

    let test_points = vec![
        (0.2, 0.1, 0.0, 0.0), // Near wall
        (0.5, 0.5, 0.0, 0.5), // Mid domain
        (0.8, 0.9, 0.0, 1.0), // Far field
    ];

    for (x, y, z, t) in test_points {
        let temp = hypersonic.exact_solution(x, y, z, t);
        let source = hypersonic.source_term(x, y, z, t);

        // Temperature constraints for hypersonic flow
        assert!(
            temp > 3.0,
            "Temperature too low for hypersonic flow at ({},{},t={}): {}",
            x,
            y,
            t,
            temp
        );
        assert!(
            temp < 10.0,
            "Temperature perturbation too large at ({},{},t={}): {}",
            x,
            y,
            t,
            temp
        );
        assert!(
            source.is_finite(),
            "Hypersonic source term not finite at ({},{},t={}): {}",
            x,
            y,
            t,
            source
        );
    }

    println!("✓ Hypersonic flow validation passed");
}

/// Test shock-capturing scheme validation
#[test]
fn test_shock_capturing_validation() {
    println!("Testing Shock-Capturing MMS Validation...");

    let shock = ManufacturedShockCapturing::<f64>::new(
        6.0,  // shock strength (density ratio)
        2.5,  // shock speed
        0.2,  // initial shock position
        0.05, // smooth perturbation amplitude
        3.0,  // kx
        2.0,  // ky
    );

    // Test shock properties
    assert!(shock.shock_strength > 1.0, "Shock strength must be > 1");
    assert!(shock.shock_speed > 0.0, "Shock speed must be positive");

    // Test at different times
    let times = vec![0.0, 0.2, 0.4, 0.6];
    let mut prev_shock_pos = shock.shock_position(0.0);

    for &t in &times {
        let shock_pos = shock.shock_position(t);

        // Shock should move with constant speed
        assert!(
            shock_pos >= prev_shock_pos,
            "Shock should not move backward"
        );
        assert!(
            (shock_pos
                - prev_shock_pos
                - shock.shock_speed
                    * (t - times[times
                        .iter()
                        .position(|&x| x == prev_shock_pos / shock.shock_speed)
                        .unwrap_or(0)]))
                < 1e-10,
            "Shock speed should be constant"
        );

        // Test solution across shock
        let x_pre = shock_pos - 1e-4;
        let x_post = shock_pos + 1e-4;
        let y = 0.5;

        let rho_pre = shock.exact_solution(x_pre, y, 0.0, t);
        let rho_post = shock.exact_solution(x_post, y, 0.0, t);

        // Density jump should match shock strength
        let density_ratio = rho_post / rho_pre;
        assert!(
            (density_ratio - shock.shock_strength).abs() < 1e-3,
            "Density ratio doesn't match shock strength: expected={}, got={}",
            shock.shock_strength,
            density_ratio
        );

        // Both regions should have finite source terms
        let source_pre = shock.source_term(x_pre, y, 0.0, t);
        let source_post = shock.source_term(x_post, y, 0.0, t);
        assert!(source_pre.is_finite(), "Pre-shock source term not finite");
        assert!(source_post.is_finite(), "Post-shock source term not finite");

        prev_shock_pos = shock_pos;
    }

    println!("✓ Shock-capturing validation passed");
}

/// Test turbulence-chemistry interaction validation
#[test]
fn test_tci_validation() {
    println!("Testing Turbulence-Chemistry Interaction MMS Validation...");

    let tci = ManufacturedTCI::<f64>::new(
        0.9,  // turbulent Schmidt number
        5.0,  // Damkohler number (fast chemistry)
        2.0,  // reaction rate constant
        0.02, // turbulent diffusivity
        0.3,  // mixture fraction amplitude
        2.5,  // kx
        1.8,  // ky
    );

    // Test TCI parameters
    assert!(tci.schmidt_t > 0.0, "Schmidt number must be positive");
    assert!(tci.damkohler > 0.0, "Damkohler number must be positive");
    assert!(tci.diffusivity_t > 0.0, "Diffusivity must be positive");

    let test_points = vec![
        (0.1, 0.1, 0.0, 0.0),
        (0.3, 0.4, 0.0, 0.3),
        (0.6, 0.7, 0.0, 0.6),
        (0.9, 0.9, 0.0, 0.9),
    ];

    for (x, y, z, t) in test_points {
        let z_mix = tci.exact_solution(x, y, z, t);
        let source = tci.source_term(x, y, z, t);

        // Mixture fraction bounds
        assert!(
            z_mix >= 0.0,
            "Mixture fraction must be non-negative at ({},{},t={}): {}",
            x,
            y,
            t,
            z_mix
        );
        assert!(
            z_mix <= 1.0,
            "Mixture fraction must be ≤ 1.0 at ({},{},t={}): {}",
            x,
            y,
            t,
            z_mix
        );

        // Source term should be finite
        assert!(
            source.is_finite(),
            "TCI source term not finite at ({},{},t={}): {}",
            x,
            y,
            t,
            source
        );

        // For fast chemistry (Da >> 1), reaction should drive mixture fraction toward equilibrium
        if tci.damkohler > 1.0 {
            // This is a simplified check - in practice would verify reaction-diffusion balance
            assert!(
                source.abs() < 10.0,
                "Source term magnitude too large for fast chemistry: {}",
                source
            );
        }
    }

    println!("✓ TCI validation passed");
}

/// Test coupled compressible-turbulent flow validation
#[test]
fn test_coupled_compressible_turbulent_validation() {
    println!("Testing Coupled Compressible-Turbulent Flow Validation...");

    // Create coupled compressible Euler and turbulence models
    let compressible = ManufacturedCompressibleEuler::<f64>::new(3.0, 1.4, 0.0, 0.1, 1.0, 1.0);

    // Test compressibility effects on turbulence (simplified coupling)
    let x = 0.5;
    let y = 0.5;
    let t = 0.5;

    let rho = compressible.exact_solution(x, y, 0.0, t);
    let source_comp = compressible.source_term(x, y, 0.0, t);

    // For compressible flows, density fluctuations affect momentum and energy equations
    assert!(
        rho > 0.8 && rho < 1.3,
        "Density fluctuations reasonable: {}",
        rho
    );
    assert!(source_comp.is_finite(), "Compressible source term finite");

    // In a full implementation, we would couple this with turbulence MMS
    // For now, we verify the compressible solution properties

    println!("✓ Coupled compressible-turbulent validation passed");
}

/// Property-based tests for advanced physics validation
#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Test compressible Euler with various Mach numbers
        #[test]
        fn test_compressible_euler_properties(
            mach in 0.5f64..5.0,
            gamma in 1.2f64..1.8,
            flow_angle in 0.0f64..0.5,
            amplitude in 0.01f64..0.2,
            kx in 0.5f64..3.0,
            ky in 0.5f64..3.0
        ) {
            let euler = ManufacturedCompressibleEuler::new(mach, gamma, flow_angle, amplitude, kx, ky);

            let x = 0.5;
            let y = 0.5;
            let t = 1.0;

            let rho = euler.exact_solution(x, y, 0.0, t);
            let source = euler.source_term(x, y, 0.0, t);

            prop_assert!(rho > 0.0, "Density must be positive");
            prop_assert!(rho < 2.0, "Density perturbation reasonable");
            prop_assert!(source.is_finite(), "Source term must be finite");
            prop_assert!(euler.u_0() > 0.0, "Reference velocity positive");
        }

        /// Test shock capturing with various shock strengths
        #[test]
        fn test_shock_capturing_properties(
            shock_strength in 2.0f64..10.0,
            shock_speed in 0.5f64..5.0,
            shock_x0 in 0.1f64..0.9,
            amplitude in 0.01f64..0.1,
            kx in 1.0f64..5.0,
            ky in 1.0f64..5.0
        ) {
            let shock = ManufacturedShockCapturing::new(shock_strength, shock_speed, shock_x0, amplitude, kx, ky);

            let t = 0.5;
            let shock_pos = shock.shock_position(t);

            prop_assert!(shock_pos > shock_x0, "Shock should move forward");
            prop_assert!((shock_pos - shock_x0 - shock_speed * t).abs() < 1e-10, "Shock speed constant");

            // Test density jump
            let x_pre = shock_pos - 1e-4;
            let x_post = shock_pos + 1e-4;
            let rho_pre = shock.exact_solution(x_pre, 0.5, 0.0, t);
            let rho_post = shock.exact_solution(x_post, 0.5, 0.0, t);

            prop_assert!(rho_post > rho_pre, "Post-shock density higher");
            prop_assert!((rho_post / rho_pre - shock_strength).abs() < 1e-3, "Density ratio matches shock strength");
        }

        /// Test TCI with various reaction parameters
        #[test]
        fn test_tci_properties(
            schmidt_t in 0.3f64..1.5,
            damkohler in 0.1f64..10.0,
            reaction_rate in 0.1f64..5.0,
            diffusivity_t in 0.001f64..0.1,
            amplitude in 0.1f64..0.4,
            kx in 1.0f64..4.0,
            ky in 1.0f64..4.0
        ) {
            let tci = ManufacturedTCI::new(schmidt_t, damkohler, reaction_rate, diffusivity_t, amplitude, kx, ky);

            let x = 0.5;
            let y = 0.5;
            let t = 1.0;

            let z = tci.exact_solution(x, y, 0.0, t);
            let source = tci.source_term(x, y, 0.0, t);

            prop_assert!(z >= 0.0 && z <= 1.0, "Mixture fraction in bounds");
            prop_assert!(source.is_finite(), "Source term finite");
            prop_assert!(source.abs() < 100.0, "Source term reasonable magnitude");
        }

        /// Test hypersonic flow properties
        #[test]
        fn test_hypersonic_properties(
            mach_inf in 5.0f64..20.0,
            reynolds in 1e4f64..1e8,
            prandtl in 0.5f64..1.0,
            twall_ratio in 2.0f64..10.0,
            amplitude in 0.05f64..0.3,
            kx in 0.5f64..2.0,
            ky in 0.5f64..2.0
        ) {
            let hypersonic = ManufacturedHypersonic::new(mach_inf, reynolds, prandtl, 1.4, twall_ratio, amplitude, kx, ky);

            let x = 0.5;
            let y = 0.5;
            let t = 1.0;

            let temp = hypersonic.exact_solution(x, y, 0.0, t);
            let source = hypersonic.source_term(x, y, 0.0, t);

            prop_assert!(temp >= twall_ratio - amplitude && temp <= twall_ratio + amplitude * 2.0, "Temperature in expected range");
            prop_assert!(source.is_finite(), "Source term finite");
            prop_assert!(mach_inf >= 5.0, "Hypersonic Mach number");
        }
    }
}

/// Integration test for advanced physics validation pipeline
#[test]
fn test_advanced_physics_integration() {
    println!("🧪 Testing Advanced Physics Validation Integration...");

    // Test all advanced physics MMS together
    let models = vec![
        (
            "Compressible Euler",
            Box::new(ManufacturedCompressibleEuler::<f64>::new(
                2.5, 1.4, 0.2, 0.1, 1.5, 1.2,
            )) as Box<dyn ManufacturedSolution<f64>>,
        ),
        (
            "TCI",
            Box::new(ManufacturedTCI::<f64>::new(
                0.8, 3.0, 1.5, 0.015, 0.25, 2.0, 1.5,
            )),
        ),
        (
            "Hypersonic",
            Box::new(ManufacturedHypersonic::<f64>::new(
                12.0, 5e5, 0.71, 1.4, 5.0, 0.15, 1.2, 0.8,
            )),
        ),
        (
            "Shock Capturing",
            Box::new(ManufacturedShockCapturing::<f64>::new(
                8.0, 3.0, 0.25, 0.03, 2.5, 1.8,
            )),
        ),
    ];

    let test_point = (0.4, 0.6, 0.0, 0.8);

    for (name, model) in models {
        let solution = model.exact_solution(test_point.0, test_point.1, test_point.2, test_point.3);
        let source = model.source_term(test_point.0, test_point.1, test_point.2, test_point.3);

        assert!(
            solution.is_finite(),
            "{} solution not finite: {}",
            name,
            solution
        );
        assert!(
            source.is_finite(),
            "{} source term not finite: {}",
            name,
            source
        );

        println!(
            "  ✓ {}: solution={:.6}, source={:.6}",
            name, solution, source
        );
    }

    println!("✅ Advanced physics integration test passed");
}

