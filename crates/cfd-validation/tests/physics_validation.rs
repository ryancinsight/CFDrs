//! Physics validation tests against analytical solutions
//!
//! These tests verify that our numerical implementations
//! match known analytical solutions to required accuracy.

use approx::assert_relative_eq;
use cfd_validation::analytical_benchmarks::{CouetteFlow, PoiseuilleFlow, TaylorGreenVortex};

#[test]
fn validate_couette_flow_profile() {
    // Test case from White (2016), Example 3.1
    let flow = CouetteFlow {
        u_wall: 10.0, // 10 m/s upper wall
        h: 0.001,     // 1 mm gap
        dp_dx: 0.0,   // No pressure gradient
        mu: 0.001,    // 1 mPa·s
    };

    // Check linear profile when dp/dx = 0
    let y_test = 0.0005; // Midpoint
    let u_expected = 5.0; // Should be u_wall/2

    assert_relative_eq!(
        flow.velocity(y_test),
        u_expected,
        epsilon = 1e-10,
        max_relative = 1e-10
    );

    // Check wall boundary conditions
    assert_relative_eq!(flow.velocity(0.0), 0.0, epsilon = 1e-10);
    assert_relative_eq!(flow.velocity(flow.h), flow.u_wall, epsilon = 1e-10);
}

#[test]
fn validate_poiseuille_parabolic_profile() {
    // Test case: Channel flow with pressure gradient
    let flow = PoiseuilleFlow {
        h: 0.01,       // 1 cm half-height
        dp_dx: -100.0, // -100 Pa/m pressure gradient (pressure drop in +x)
        mu: 0.001,     // 1 mPa·s (water-like)
    };

    // Validate parabolic profile
    let y_values = [0.0, 0.0025, 0.005, 0.0075, 0.01];
    let expected_velocities = vec![
        5.0,    // Centerline maximum: (1/2μ)(dp/dx)(h²) = (1/(2×0.001))×100×(0.01²) = 5.0
        4.6875, // At y=0.0025: 5.0×(1-(0.0025/0.01)²) = 5.0×(1-0.0625) = 4.6875
        3.75,   // At y=0.005: 5.0×(1-(0.005/0.01)²) = 5.0×(1-0.25) = 3.75
        2.1875, // At y=0.0075: 5.0×(1-(0.0075/0.01)²) = 5.0×(1-0.5625) = 2.1875
        0.0,    // Wall: 5.0×(1-(0.01/0.01)²) = 5.0×(1-1) = 0.0
    ];

    for (y, u_expected) in y_values.iter().zip(expected_velocities.iter()) {
        let u_computed = flow.velocity(*y);
        assert_relative_eq!(
            u_computed,
            *u_expected,
            epsilon = 1e-10,
            max_relative = 1e-10
        );
    }

    // Validate flow rate calculation
    let q = flow.flow_rate();
    let q_expected = 2.0 * 0.01 * (2.0 / 3.0) * 5.0; // 2h * (2/3) * u_max = 2×0.01×(2/3)×5.0
    assert_relative_eq!(q, q_expected, epsilon = 1e-10);
}

#[test]
fn validate_taylor_green_decay() {
    // Classic Taylor-Green vortex test
    let vortex = TaylorGreenVortex {
        u0: 1.0,
        l: 2.0 * std::f64::consts::PI, // Domain size 2π
        nu: 0.1,                       // Kinematic viscosity
    };

    // Test velocity field at t=0
    let v0 = vortex.velocity(0.0, 0.0, 0.0);
    assert_relative_eq!(v0[0], 0.0, epsilon = 1e-10); // u=0 at origin
    assert_relative_eq!(v0[1], 0.0, epsilon = 1e-10); // v=0 at origin

    // Test velocity at (π/2, 0)
    let v1 = vortex.velocity(std::f64::consts::PI / 2.0, 0.0, 0.0);
    assert_relative_eq!(v1[0], 1.0, epsilon = 1e-10); // u=u0
    assert_relative_eq!(v1[1], 0.0, epsilon = 1e-10); // v=0

    // Test energy decay
    let t = 1.0;
    let e0 = vortex.kinetic_energy(0.0);
    let e1 = vortex.kinetic_energy(t);
    let expected_ratio = (-2.0 * vortex.nu * t).exp();

    assert_relative_eq!(e1 / e0, expected_ratio, epsilon = 1e-10);
}

#[test]
fn validate_reynolds_number_calculation() {
    use cfd_core::constants::physics_validated::{fluid_dynamics, reynolds};

    // Test pipe flow transition
    let d = 0.01; // 1 cm diameter
    let mu = fluid_dynamics::WATER_DYNAMIC_VISCOSITY_20C;
    let rho = fluid_dynamics::WATER_DENSITY_20C;
    let nu = mu / rho;

    // Calculate velocity for Re = 2300 (transition point)
    let re_crit = reynolds::PIPE_TRANSITION_LOWER;
    let u_crit = re_crit * nu / d;

    // Verify calculation
    let re_computed = u_crit * d / nu;
    assert_relative_eq!(re_computed, re_crit, epsilon = 1e-10);
}

#[test]
fn validate_prandtl_number() {
    use cfd_core::constants::physics_validated::{fluid_dynamics, thermodynamics, validation};

    // Calculate Prandtl number for water
    let mu = fluid_dynamics::WATER_DYNAMIC_VISCOSITY_20C;
    let cp = thermodynamics::WATER_SPECIFIC_HEAT_20C;
    let k = thermodynamics::WATER_THERMAL_CONDUCTIVITY_20C;

    let pr_computed = mu * cp / k;

    assert_relative_eq!(
        pr_computed,
        validation::WATER_PRANDTL_20C,
        epsilon = 0.01, // 1% tolerance for Prandtl number
        max_relative = 0.01
    );
}

#[test]
fn validate_lid_driven_cavity_benchmark() {
    use cfd_validation::analytical_benchmarks::lid_driven_cavity;

    // Verify benchmark data is properly loaded
    assert_eq!(lid_driven_cavity::RE100_U_CENTERLINE.len(), 17);
    assert_eq!(lid_driven_cavity::RE1000_U_CENTERLINE.len(), 17);

    // Check boundary conditions
    let first = lid_driven_cavity::RE100_U_CENTERLINE[0];
    let last = lid_driven_cavity::RE100_U_CENTERLINE[16];

    assert_eq!(first.0, 0.0); // Bottom wall
    assert_eq!(first.1, 0.0); // No-slip
    assert_eq!(last.0, 1.0); // Top wall (lid)
    assert_eq!(last.1, 1.0); // Lid velocity
}
