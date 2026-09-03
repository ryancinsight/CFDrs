#![allow(missing_docs)]
//! Physics validation tests against analytical solutions
//!
//! These tests verify that our numerical implementations
//! match known analytical solutions to required accuracy.

use aequitas::systems::si::quantities::{
    DynamicViscosity, KinematicViscosity, Length, PressureGradient, Time, Velocity,
};
use cfd_validation::analytical::poiseuille::{PoiseuilleFlowRate, PoiseuilleGeometry};
use cfd_validation::analytical::taylor_green::TaylorGreenKineticEnergy;
use cfd_validation::analytical::{
    AnalyticalSolution, CouetteFlow, PoiseuilleFlow, TaylorGreenVortex,
};
use cfd_validation::benchmarks::cavity::LidDrivenCavity;
use eunomia::assert_relative_eq;

#[test]
fn validate_couette_flow_profile() {
    // Test case from White (2016), Example 3.1
    let flow = CouetteFlow::create(
        Velocity::from_base(10.0),
        Length::from_base(0.001),
        PressureGradient::from_base(0.0),
        DynamicViscosity::from_base(0.001),
    );

    // Check linear profile when dp/dx = 0
    let y_test = 0.0005; // Midpoint
    let u_expected = 5.0; // Should be u_wall/2

    assert_relative_eq!(
        flow.evaluate(0.0, y_test, 0.0, 0.0)[0],
        u_expected,
        epsilon = 1e-10,
        max_relative = 1e-10
    );

    // Check wall boundary conditions
    assert_relative_eq!(flow.evaluate(0.0, 0.0, 0.0, 0.0)[0], 0.0, epsilon = 1e-10);
    assert_relative_eq!(
        flow.evaluate(0.0, flow.gap_height.into_base(), 0.0, 0.0)[0],
        flow.wall_velocity.into_base(),
        epsilon = 1e-10
    );
}

#[test]
fn validate_poiseuille_parabolic_profile() {
    // Test case: Channel flow with pressure gradient
    let flow = PoiseuilleFlow::create(
        Velocity::from_base(5.0),
        Length::from_base(0.01),
        PressureGradient::from_base(100.0),
        DynamicViscosity::from_base(0.001),
        PoiseuilleGeometry::Plates,
    );

    // Validate parabolic profile
    let y_values = [0.0, 0.0025, 0.005, 0.0075, 0.01];
    let expected_velocities = [
        5.0,    // Centerline maximum: (1/2μ)(dp/dx)(h²) = (1/(2×0.001))×100×(0.01²) = 5.0
        4.6875, // At y=0.0025: 5.0×(1-(0.0025/0.01)²) = 5.0×(1-0.0625) = 4.6875
        3.75,   // At y=0.005: 5.0×(1-(0.005/0.01)²) = 5.0×(1-0.25) = 3.75
        2.1875, // At y=0.0075: 5.0×(1-(0.0075/0.01)²) = 5.0×(1-0.5625) = 2.1875
        0.0,    // Wall: 5.0×(1-(0.01/0.01)²) = 5.0×(1-1) = 0.0
    ];

    for (y, u_expected) in y_values.iter().zip(expected_velocities.iter()) {
        let u_computed = flow.evaluate(0.0, *y, 0.0, 0.0)[0];
        assert_relative_eq!(
            u_computed,
            *u_expected,
            epsilon = 1e-10,
            max_relative = 1e-10
        );
    }

    // Validate flow rate calculation
    let PoiseuilleFlowRate::PerWidth(q) = flow.flow_rate() else {
        panic!("plate flow must return a per-width rate");
    };
    let q_expected = 2.0 * 0.01 * (2.0 / 3.0) * 5.0; // 2h * (2/3) * u_max = 2×0.01×(2/3)×5.0
    assert_relative_eq!(q.into_base(), q_expected, epsilon = 1e-10);
}

#[test]
fn validate_taylor_green_decay() {
    // Classic Taylor-Green vortex test
    let vortex = TaylorGreenVortex::create_2d(
        Length::from_base(2.0 * std::f64::consts::PI),
        Velocity::from_base(1.0),
        KinematicViscosity::from_base(0.1),
    );

    // Test velocity field at t=0
    let v0 = vortex.evaluate(0.0, 0.0, 0.0, 0.0);
    assert_relative_eq!(v0[0], 0.0, epsilon = 1e-10); // u=0 at origin
    assert_relative_eq!(v0[1], 0.0, epsilon = 1e-10); // v=0 at origin

    // Test velocity at (0, π), where the canonical 2D solution has u=u0.
    let v1 = vortex.evaluate(0.0, std::f64::consts::PI, 0.0, 0.0);
    assert_relative_eq!(v1[0], 1.0, epsilon = 1e-10); // u=u0
    assert_relative_eq!(v1[1], 0.0, epsilon = 1e-10); // v=0

    // Test energy decay. The implementation pins the k = π/L convention
    // (u = U·cos(kx)·sin(ky), asserted by the velocity checks above), under
    // which velocity decays as e^{−2νk²t} and kinetic energy as e^{−4νk²t}.
    // With L = 2π, ν = 0.1: k = 1/2, so E(t)/E(0) = e^{−0.1·t}.
    let t = 1.0;
    let e0 = vortex.kinetic_energy(Time::from_base(0.0));
    let e1 = vortex.kinetic_energy(Time::from_base(t));
    let k = std::f64::consts::PI / (2.0 * std::f64::consts::PI);
    let expected_ratio = (-4.0 * 0.1 * k * k * t).exp();

    let TaylorGreenKineticEnergy::PerDepth(e0) = e0 else {
        panic!("2D Taylor-Green energy must be reported per depth");
    };
    let TaylorGreenKineticEnergy::PerDepth(e1) = e1 else {
        panic!("2D Taylor-Green energy must be reported per depth");
    };
    assert_relative_eq!(
        e1.into_base() / e0.into_base(),
        expected_ratio,
        epsilon = 1e-10
    );
}

#[test]
fn validate_reynolds_number_calculation() {
    use cfd_core::physics::constants::physics::{dimensionless::reynolds, fluid};

    // Test pipe flow transition
    let d = 0.01; // 1 cm diameter
    let mu = fluid::WATER_VISCOSITY;
    let rho = fluid::WATER_DENSITY;
    let nu = mu / rho;

    // Calculate velocity for Re = 2300 (transition point)
    let re_crit = reynolds::PIPE_LAMINAR_MAX;
    let u_crit = re_crit * nu / d;

    // Verify calculation
    let re_computed = u_crit * d / nu;
    assert_relative_eq!(re_computed, re_crit, epsilon = 1e-10);
}

#[test]
fn validate_prandtl_number() {
    use cfd_core::physics::constants::physics::fluid;

    // Calculate Prandtl number for water
    let mu = fluid::WATER_VISCOSITY;
    let cp = fluid::WATER_SPECIFIC_HEAT;
    let k = fluid::WATER_THERMAL_CONDUCTIVITY;

    let pr_computed = mu * cp / k;

    // Standard Prandtl number for water at 20°C is ~7.0
    let water_prandtl_20c = 7.01;

    assert_relative_eq!(
        pr_computed,
        water_prandtl_20c,
        epsilon = 0.15, // Allow for variations in standard property tables
        max_relative = 0.02
    );
}

#[test]
fn validate_lid_driven_cavity_benchmark() {
    let cavity = LidDrivenCavity::new(1.0_f64, 1.0_f64, 100.0_f64);
    let re100 = cavity.ghia_u_centerline(100.0);
    let re400 = cavity.ghia_u_centerline(400.0);
    let re1000 = cavity.ghia_u_centerline(1000.0);

    assert_eq!(re100.len(), 17);
    assert_eq!(re400.len(), 17);
    assert_eq!(re1000.len(), 17);

    // Boundary values and interior extrema are published Table I values.
    assert_eq!(re100.first(), Some(&(1.0, 1.0)));
    assert_eq!(re100.last(), Some(&(0.0, 0.0)));
    assert_eq!(re100.get(9), Some(&(0.4531, -0.21090)));
    assert_eq!(re400.get(10), Some(&(0.2813, -0.32726)));
    assert_eq!(re1000.get(11), Some(&(0.1719, -0.38289)));
}
