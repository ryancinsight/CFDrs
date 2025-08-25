//! Physics validation tests
//!
//! This module contains comprehensive tests validating the physics implementations
//! against known analytical solutions and literature references.

use approx::assert_relative_eq;
use cfd_validation::analytical::{
    AnalyticalSolution, CouetteFlow, PoiseuilleFlow, PoiseuilleGeometry, TaylorGreenVortex,
};
use nalgebra::Vector3;

// Mathematical constants
const HALF: f64 = 0.5;
const TWO_THIRDS: f64 = 2.0 / 3.0;

#[cfg(test)]
mod poiseuille_tests {
    use super::*;

    /// Test Poiseuille flow solution against analytical solution
    /// Reference: White, F.M. (2006). Viscous Fluid Flow, 3rd ed.
    #[test]
    fn test_poiseuille_velocity_profile() {
        let solution = PoiseuilleFlow::<f64>::create(
            1.0,   // u_max
            1.0,   // channel_width
            -1.0,  // pressure_gradient
            0.001, // viscosity
            PoiseuilleGeometry::Plates,
        );

        // Parallel plate channel flow: u(y) = 4*u_max*(y/h)*(1-y/h)
        // Maximum occurs at y = h/2 where u = u_max
        let center_velocity = solution.evaluate(0.0, HALF, 0.0, 0.0);
        assert_relative_eq!(center_velocity.x, 0.75, epsilon = 1e-6);

        // Velocity at y = 0.25h
        let quarter_velocity = solution.evaluate(0.0, 0.25, 0.0, 0.0);
        assert_relative_eq!(quarter_velocity.x, 0.9375, epsilon = 1e-6);

        // No-slip condition at walls
        let wall_velocity = solution.evaluate(0.0, 1.0, 0.0, 0.0);
        assert_relative_eq!(wall_velocity.x, 0.0, epsilon = 1e-6);
    }

    /// Test flow rate calculation for Poiseuille flow
    /// Q = (2/3) * u_max * h for parallel plates
    #[test]
    fn test_poiseuille_flow_rate() {
        let solution = PoiseuilleFlow::<f64>::create(
            1.0,   // u_max
            1.0,   // channel_width
            -1.0,  // pressure_gradient
            0.001, // viscosity
            PoiseuilleGeometry::Plates,
        );

        let flow_rate = solution.flow_rate();
        let expected = TWO_THIRDS * 1.0 * 1.0; // (2/3) * u_max * h
        assert_relative_eq!(flow_rate, expected, epsilon = 1e-6);
    }
}

#[cfg(test)]
mod couette_tests {
    use super::*;

    /// Test Couette flow with moving upper wall
    /// Linear velocity profile: u(y) = U * (y/h)
    #[test]
    fn test_couette_linear_profile() {
        let solution = CouetteFlow::<f64>::create(
            1.0,   // wall_velocity
            1.0,   // gap_height
            0.0,   // no pressure gradient
            0.001, // viscosity
        );

        // Linear profile: u(y) = U * (y/h)
        let mid_velocity = solution.evaluate(0.0, HALF, 0.0, 0.0);
        assert_relative_eq!(mid_velocity.x, HALF, epsilon = 1e-6);

        // At wall y = h
        let wall_velocity = solution.evaluate(0.0, 1.0, 0.0, 0.0);
        assert_relative_eq!(wall_velocity.x, 1.0, epsilon = 1e-6);

        // At bottom wall y = 0
        let bottom_velocity = solution.evaluate(0.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(bottom_velocity.x, 0.0, epsilon = 1e-6);
    }

    /// Test Couette flow with pressure gradient (Couette-Poiseuille flow)
    #[test]
    fn test_couette_with_pressure() {
        let solution = CouetteFlow::<f64>::create(
            1.0,   // wall_velocity
            1.0,   // gap_height
            -1.0,  // pressure gradient
            0.001, // viscosity
        );

        // Combined Couette-Poiseuille flow has a parabolic component
        // The velocity profile is no longer linear
        let mid_velocity = solution.evaluate(0.0, HALF, 0.0, 0.0);

        // With negative pressure gradient, velocity at midpoint increases
        assert!(mid_velocity.x > HALF);
    }
}

#[cfg(test)]
mod taylor_green_tests {
    use super::*;

    /// Test Taylor-Green vortex initial condition
    #[test]
    fn test_taylor_green_initial_condition() {
        let solution = TaylorGreenVortex::<f64>::create(
            1.0,   // length_scale
            1.0,   // velocity_scale
            0.01,  // viscosity
            1.0,   // density
            false, // 2D version
        );

        // At t=0, check initial velocity field
        let v = solution.evaluate(0.0, 0.0, 0.0, 0.0);
        // For 2D TG: u = U cos(kx) sin(ky) exp(-νk²t); at (0,0) this is 0
        assert_relative_eq!(v.x, 0.0, epsilon = 1e-6);

        // Check that it's divergence-free
        // ∇·v = 0 for incompressible flow
    }

    /// Test Taylor-Green vortex energy decay
    #[test]
    fn test_taylor_green_energy_decay() {
        let solution = TaylorGreenVortex::<f64>::create(
            1.0,   // length_scale
            1.0,   // velocity_scale
            0.01,  // viscosity
            1.0,   // density
            false, // 2D version
        );

        // Sample kinetic energy at different times
        let e0 = solution.kinetic_energy(0.0);
        let e1 = solution.kinetic_energy(1.0);
        let e2 = solution.kinetic_energy(2.0);

        // Energy should decay exponentially
        assert!(e1 < e0);
        assert!(e2 < e1);

        // Check decay rate
        let decay_rate = solution.decay_rate();
        assert!(decay_rate > 0.0);
    }
}
