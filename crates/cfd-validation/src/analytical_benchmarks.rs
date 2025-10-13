//! Analytical benchmarks for CFD validation
//!
//! These are exact solutions to the Navier-Stokes equations
//! used to validate numerical implementations.

use cfd_core::constants::mathematical::{
    numeric::{ONE_HALF, TWO, TWO_THIRDS},
    TWO_PI,
};
use cfd_core::conversion::SafeFromF64;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Couette flow: Shear-driven flow between parallel plates
///
/// Reference: White, F.M. (2016). Viscous Fluid Flow, 3rd ed., Section 3.2
pub struct CouetteFlow<T: RealField + Copy> {
    /// Upper plate velocity [m/s]
    pub u_wall: T,
    /// Gap height [m]
    pub h: T,
    /// Pressure gradient dp/dx [Pa/m]
    pub dp_dx: T,
    /// Dynamic viscosity [Pa·s]
    pub mu: T,
}

impl<T: RealField + Copy + FromPrimitive> CouetteFlow<T> {
    /// Exact velocity profile: u(y) = U*y/h + (h²/2μ)(dp/dx)(y/h)(1-y/h)
    ///
    /// Source: Schlichting & Gersten (2017), Eq. 5.10
    pub fn velocity(&self, y: T) -> T {
        let y_norm = y / self.h;
        let linear_term = self.u_wall * y_norm;
        let two = T::from_f64_or_one(TWO);
        let pressure_term =
            (self.h * self.h / (two * self.mu)) * self.dp_dx * y_norm * (T::one() - y_norm);
        linear_term + pressure_term
    }

    /// Wall shear stress
    pub fn wall_shear(&self) -> T {
        let two = T::from_f64_or_one(TWO);
        self.mu * self.u_wall / self.h - self.h * self.dp_dx / two
    }
}

/// Poiseuille flow: Pressure-driven flow in a channel
///
/// Reference: Batchelor, G.K. (2000). An Introduction to Fluid Dynamics, Section 4.2
pub struct PoiseuilleFlow<T: RealField + Copy> {
    /// Channel half-height [m]
    pub h: T,
    /// Pressure gradient dp/dx [Pa/m] (negative for flow in +x direction)
    pub dp_dx: T,
    /// Dynamic viscosity [Pa·s]
    pub mu: T,
}

impl<T: RealField + Copy + FromPrimitive> PoiseuilleFlow<T> {
    /// Exact velocity profile: u(y) = -(1/2μ)(dp/dx)(h² - y²)
    ///
    /// The negative sign ensures positive velocity when pressure decreases
    /// in the flow direction (dp/dx < 0).
    ///
    /// Source: White (2016), Eq. 3.36
    pub fn velocity(&self, y: T) -> T {
        let half = T::from_f64_or_one(ONE_HALF);
        let factor = half / self.mu;
        -factor * self.dp_dx * (self.h * self.h - y * y)
    }

    /// Maximum velocity (at centerline)
    pub fn max_velocity(&self) -> T {
        self.velocity(T::zero())
    }

    /// Average velocity: `u_avg` = (2/3) * `u_max`
    pub fn average_velocity(&self) -> T {
        let two_thirds = T::from_f64_or_one(TWO_THIRDS);
        two_thirds * self.max_velocity()
    }

    /// Volume flow rate per unit width
    pub fn flow_rate(&self) -> T {
        let two = T::from_f64_or_one(TWO);
        two * self.h * self.average_velocity()
    }
}

/// Taylor-Green vortex: Exact unsteady solution
///
/// Reference: Taylor & Green (1937). Proc. R. Soc. Lond. A, 158(895), 499-521
pub struct TaylorGreenVortex<T: RealField + Copy> {
    /// Characteristic velocity [m/s]
    pub u0: T,
    /// Characteristic length [m]
    pub l: T,
    /// Kinematic viscosity [m²/s]
    pub nu: T,
}

impl<T: RealField + Copy + FromPrimitive> TaylorGreenVortex<T> {
    /// Exact velocity field at time t
    pub fn velocity(&self, x: T, y: T, t: T) -> Vector3<T> {
        let k = T::from_f64_or_one(TWO_PI) / self.l;
        let decay = (-k * k * self.nu * t).exp();

        let u = self.u0 * (k * x).sin() * (k * y).cos() * decay;
        let v = -self.u0 * (k * x).cos() * (k * y).sin() * decay;

        Vector3::new(u, v, T::zero())
    }

    /// Exact pressure field
    pub fn pressure(&self, x: T, y: T, t: T, rho: T) -> T {
        let k = T::from_f64_or_one(TWO_PI) / self.l;
        let two = T::from_f64_or_one(TWO);
        let decay = (-two * k * k * self.nu * t).exp();

        let quarter = T::from_f64_or_zero(0.25);
        -quarter * rho * self.u0 * self.u0 * ((two * k * x).cos() + (two * k * y).cos()) * decay
    }

    /// Kinetic energy decay rate
    pub fn kinetic_energy(&self, t: T) -> T {
        let k = T::from_f64_or_one(TWO_PI) / self.l;
        let two = T::from_f64_or_one(TWO);
        let decay = (-two * k * k * self.nu * t).exp();

        let half = T::from_f64_or_one(ONE_HALF);
        half * self.u0 * self.u0 * decay
    }
}

/// Lid-driven cavity benchmark points
///
/// Reference: Ghia et al. (1982). J. Comput. Phys., 48(3), 387-411
pub mod lid_driven_cavity {
    /// Benchmark data for Re=100
    pub const RE100_U_CENTERLINE: &[(f64, f64)] = &[
        (0.0000, 0.00000),
        (0.0625, -0.03717),
        (0.1250, -0.04192),
        (0.1875, -0.04775),
        (0.2500, -0.05641),
        (0.3125, -0.04299),
        (0.3750, -0.02388),
        (0.4375, -0.00737),
        (0.5000, 0.00669),
        (0.5625, 0.01911),
        (0.6250, 0.03055),
        (0.6875, 0.04108),
        (0.7500, 0.05052),
        (0.8125, 0.05864),
        (0.8750, 0.06500),
        (0.9375, 0.06898),
        (1.0000, 1.00000),
    ];

    /// Benchmark data for Re=1000
    pub const RE1000_U_CENTERLINE: &[(f64, f64)] = &[
        (0.0000, 0.00000),
        (0.0625, -0.18109),
        (0.1250, -0.20196),
        (0.1875, -0.22220),
        (0.2500, -0.24533),
        (0.3125, -0.14612),
        (0.3750, -0.10150),
        (0.4375, -0.06547),
        (0.5000, -0.03827),
        (0.5625, -0.01860),
        (0.6250, -0.00570),
        (0.6875, 0.00190),
        (0.7500, 0.00570),
        (0.8125, 0.00760),
        (0.8750, 0.00950),
        (0.9375, 0.01040),
        (1.0000, 1.00000),
    ];
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_poiseuille_flow_conservation() {
        let flow = PoiseuilleFlow {
            h: 0.01,      // 1 cm half-height
            dp_dx: 100.0, // 100 Pa/m
            mu: 1e-3,     // Water viscosity
        };

        // Check that average velocity is 2/3 of max
        let u_max = flow.max_velocity();
        let u_avg = flow.average_velocity();
        assert_relative_eq!(u_avg, (2.0 / 3.0) * u_max, epsilon = 1e-10);

        // Check symmetry: u(y) = u(-y)
        let y_test = 0.005;
        assert_relative_eq!(
            flow.velocity(y_test),
            flow.velocity(-y_test),
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_taylor_green_energy_decay() {
        let vortex = TaylorGreenVortex {
            u0: 1.0,
            l: 1.0,
            nu: 0.01,
        };

        // Energy should decay exponentially
        let t1 = 0.1;
        let t2 = 0.2;
        let e1 = vortex.kinetic_energy(t1);
        let e2 = vortex.kinetic_energy(t2);

        let expected_ratio = (-2.0 * 4.0 * std::f64::consts::PI.powi(2) * 0.01 * (t2 - t1)).exp();
        assert_relative_eq!(e2 / e1, expected_ratio, epsilon = 1e-10);
    }
}
