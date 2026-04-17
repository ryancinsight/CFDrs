//! Womersley velocity profile calculator.
//!
//! Computes the complete unsteady velocity field for pulsatile pipe flow
//! using the exact analytical Bessel function solution.

use super::WomersleyNumber;
use crate::physics::vascular::bessel::{bessel_j0, bessel_j1};
use nalgebra::{Complex, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Womersley velocity profile calculator
///
/// Computes the complete unsteady velocity field for pulsatile pipe flow.
///
/// # Exact Analytical Form
///
/// The implementation evaluates the closed-form Bessel solution directly.
/// The low-$\alpha$ and high-$\alpha$ expressions are retained only as
/// analytical limits for interpretation and tests, not as runtime branches.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WomersleyProfile<T: RealField + Copy> {
    /// Womersley number parameters
    pub womersley: WomersleyNumber<T>,
    /// Pressure gradient amplitude [Pa/m]
    pub pressure_amplitude: T,
}

impl<T: RealField + FromPrimitive + Copy> WomersleyProfile<T> {
    /// Create velocity profile calculator
    pub fn new(womersley: WomersleyNumber<T>, pressure_amplitude: T) -> Self {
        Self {
            womersley,
            pressure_amplitude,
        }
    }

    /// Calculate velocity at radial position and time using exact Bessel functions
    ///
    /// # Arguments
    /// * `xi` - Dimensionless radial position r/R (0 ≤ xi ≤ 1)
    /// * `t` - Time [s]
    ///
    /// # Returns
    /// Axial velocity u(r,t) [m/s]
    pub fn velocity(&self, xi: T, t: T) -> T {
        let alpha = self.womersley.value();
        let rho = self.womersley.density;
        let omega = self.womersley.omega;
        let p_hat = self.pressure_amplitude;
        let one = T::one();

        // Clamp xi to valid range
        let xi = if xi < T::zero() {
            T::zero()
        } else if xi > one {
            one
        } else {
            xi
        };

        // i^{3/2} = e^{i 3pi/4} = (-1 + i) / sqrt(2)
        let sqrt2 = (T::one() + T::one()).sqrt();
        let i_3_2 = Complex::new(-one / sqrt2, one / sqrt2);

        // z = i^{3/2} * alpha
        let z = i_3_2 * alpha;
        let z_xi = z * xi;

        let j0_z = bessel_j0(z);
        let j0_z_xi = bessel_j0(z_xi);

        let ratio = j0_z_xi / j0_z;
        let term_brackets = Complex::new(one, T::zero()) - ratio;

        // P_hat / (i * rho * omega) = -i * P_hat / (rho * omega)
        let coeff = Complex::new(T::zero(), -p_hat / (rho * omega));

        // e^{i \omega t} = cos(\omega t) + i \sin(\omega t)
        let phase = omega * t;
        let exp_iwt = Complex::new(phase.cos(), phase.sin());

        // Final: Re{ coeff * term_brackets * exp_iwt }
        (coeff * term_brackets * exp_iwt).re
    }

    /// Calculate centerline velocity (maximum velocity)
    pub fn centerline_velocity(&self, t: T) -> T {
        self.velocity(T::zero(), t)
    }

    /// Calculate wall shear stress using exact Bessel functions
    ///
    /// τ_w(t) = -μ · (∂u/∂r)|_{r=R}
    pub fn wall_shear_stress(&self, t: T) -> T {
        let alpha = self.womersley.value();
        let r = self.womersley.radius;
        let rho = self.womersley.density;
        let mu = self.womersley.viscosity;
        let omega = self.womersley.omega;
        let p_hat = self.pressure_amplitude;
        let one = T::one();

        let sqrt2 = (T::one() + T::one()).sqrt();
        let i_3_2 = Complex::new(-one / sqrt2, one / sqrt2);

        let z = i_3_2 * alpha;
        let j0_z = bessel_j0(z);
        let j1_z = bessel_j1(z);

        // z * J_1(z) / J_0(z)
        let term = z * j1_z / j0_z;

        // P_hat / (i * rho * omega)
        let coeff = Complex::new(T::zero(), -p_hat / (rho * omega));

        let phase = omega * t;
        let exp_iwt = Complex::new(phase.cos(), phase.sin());

        // du/dxi at xi=1
        let du_dxi = (coeff * term * exp_iwt).re;

        // tau_w = -mu / R * du/dxi
        -mu / r * du_dxi
    }

    /// Calculate volumetric flow rate Q(t) using exact Bessel functions
    pub fn flow_rate(&self, t: T) -> T {
        let alpha = self.womersley.value();
        let r = self.womersley.radius;
        let rho = self.womersley.density;
        let omega = self.womersley.omega;
        let p_hat = self.pressure_amplitude;
        let one = T::one();
        let two = T::one() + T::one();
        let pi = T::pi();

        let sqrt2 = two.sqrt();
        let i_3_2 = Complex::new(-one / sqrt2, one / sqrt2);

        let z = i_3_2 * alpha;
        let j0_z = bessel_j0(z);
        let j1_z = bessel_j1(z);

        // 2 * J_1(z) / (z * J_0(z))
        let complex_two = Complex::new(two, T::zero());
        let term = complex_two * j1_z / (z * j0_z);
        let bracket = Complex::new(one, T::zero()) - term;

        // P_hat / (i * rho * omega)
        let coeff = Complex::new(T::zero(), -p_hat / (rho * omega));

        let phase = omega * t;
        let exp_iwt = Complex::new(phase.cos(), phase.sin());

        // Q = pi * R^2 * Re{ coeff * bracket * exp_iwt }
        pi * r * r * (coeff * bracket * exp_iwt).re
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use approx::assert_relative_eq;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_womersley_no_slip_invariant(alpha in 0.1..50.0_f64, t in 0.0..100.0_f64) {
            let r = 0.01_f64;
            let rho = 1000.0_f64;
            let mu = 0.003_f64;
            // Back-calculate omega from alpha: alpha = R * sqrt(rho * omega / mu)
            // omega = (alpha^2 * mu) / (R^2 * rho)
            let omega = (alpha * alpha * mu) / (r * r * rho);

            let womersley = WomersleyNumber::new(r, rho, mu, omega);
            let profile = WomersleyProfile::new(womersley, 100.0);

            // Wall location is mathematically defined at radial vector xi = r/R = 1.0
            let u_wall = profile.velocity(1.0, t);

            // Physics Theorem: Velocity uniformly zero at solid boundary (No-Slip Condition)
            assert_relative_eq!(u_wall, 0.0, epsilon = 1e-10);
        }
    }
}
