//! Womersley velocity profile calculator.
//!
//! Computes the complete unsteady velocity field for pulsatile pipe flow
//! using the exact analytical Bessel function solution.

use super::WomersleyNumber;
use crate::physics::vascular::bessel::{bessel_j0, bessel_j0_j1};
use aequitas::systems::si::quantities::{
    Pressure, PressureGradient, Time, Velocity, VolumetricFlowRate,
};
use cfd_core::CfdScalar;
use eunomia::{Complex, FloatElement, NumericElement};
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
pub struct WomersleyProfile<T: CfdScalar + Copy> {
    /// Womersley number parameters
    pub womersley: WomersleyNumber<T>,
    /// Pressure gradient amplitude [Pa/m]
    pub pressure_amplitude: PressureGradient<T>,
}

impl<T: CfdScalar + FloatElement + Copy> WomersleyProfile<T> {
    /// Create velocity profile calculator
    pub fn new(womersley: WomersleyNumber<T>, pressure_amplitude: PressureGradient<T>) -> Self {
        Self {
            womersley,
            pressure_amplitude,
        }
    }

    /// Calculate velocity at radial position and time using exact Bessel functions
    ///
    /// # Arguments
    /// * `xi` - Dimensionless radial position r/R (0 ≤ xi ≤ 1)
    /// * `t` - Time \[s]
    ///
    /// # Returns
    /// Axial velocity u(r,t) \[m/s]
    pub fn velocity(&self, xi: T, t: Time<T>) -> Velocity<T> {
        let alpha = self.womersley.value().into_base();
        let rho = self.womersley.density.into_base();
        let omega = self.womersley.omega.into_base();
        let p_hat = self.pressure_amplitude.into_base();
        let one = T::ONE;

        // Clamp xi to valid range
        let xi = if xi < T::ZERO {
            T::ZERO
        } else if xi > one {
            one
        } else {
            xi
        };

        // i^{3/2} = e^{i 3pi/4} = (-1 + i) / sqrt(2)
        let sqrt2 = <T as NumericElement>::sqrt(T::ONE + T::ONE);
        let i_3_2 = Complex::new(-one / sqrt2, one / sqrt2);

        // z = i^{3/2} * alpha
        let z = i_3_2 * alpha;
        let z_xi = z * xi;

        let j0_z = bessel_j0(z);
        let j0_z_xi = bessel_j0(z_xi);

        let ratio = j0_z_xi / j0_z;
        let term_brackets = Complex::new(one, T::ZERO) - ratio;

        // P_hat / (i * rho * omega) = -i * P_hat / (rho * omega)
        let coeff = Complex::new(T::ZERO, -p_hat / (rho * omega));

        // e^{i \omega t} = cos(\omega t) + i \sin(\omega t)
        let phase = omega * t.into_base();
        let exp_iwt = Complex::new(
            <T as FloatElement>::cos(phase),
            <T as FloatElement>::sin(phase),
        );

        // Final: Re{ coeff * term_brackets * exp_iwt }
        Velocity::from_base((coeff * term_brackets * exp_iwt).re)
    }

    /// Calculate centerline velocity (maximum velocity)
    pub fn centerline_velocity(&self, t: Time<T>) -> Velocity<T> {
        self.velocity(T::ZERO, t)
    }

    /// Calculate wall shear stress using exact Bessel functions
    ///
    /// `τ_w(t)` = -μ · (∂u/∂r)|_{r=R}
    pub fn wall_shear_stress(&self, t: Time<T>) -> Pressure<T> {
        let alpha = self.womersley.value().into_base();
        let r = self.womersley.radius.into_base();
        let rho = self.womersley.density.into_base();
        let mu = self.womersley.viscosity.into_base();
        let omega = self.womersley.omega.into_base();
        let p_hat = self.pressure_amplitude.into_base();
        let one = T::ONE;

        let sqrt2 = <T as NumericElement>::sqrt(T::ONE + T::ONE);
        let i_3_2 = Complex::new(-one / sqrt2, one / sqrt2);

        let z = i_3_2 * alpha;
        let (j0_z, j1_z) = bessel_j0_j1(z);

        // z * J_1(z) / J_0(z)
        let term = z * j1_z / j0_z;

        // P_hat / (i * rho * omega)
        let coeff = Complex::new(T::ZERO, -p_hat / (rho * omega));

        let phase = omega * t.into_base();
        let exp_iwt = Complex::new(
            <T as FloatElement>::cos(phase),
            <T as FloatElement>::sin(phase),
        );

        // du/dxi at xi=1
        let du_dxi = (coeff * term * exp_iwt).re;

        // tau_w = -mu / R * du/dxi
        Pressure::from_base(-mu / r * du_dxi)
    }

    /// Calculate volumetric flow rate Q(t) using exact Bessel functions
    pub fn flow_rate(&self, t: Time<T>) -> VolumetricFlowRate<T> {
        let alpha = self.womersley.value().into_base();
        let r = self.womersley.radius.into_base();
        let rho = self.womersley.density.into_base();
        let omega = self.womersley.omega.into_base();
        let p_hat = self.pressure_amplitude.into_base();
        let one = T::ONE;
        let two = T::ONE + T::ONE;
        let pi = T::pi();

        let sqrt2 = <T as NumericElement>::sqrt(two);
        let i_3_2 = Complex::new(-one / sqrt2, one / sqrt2);

        let z = i_3_2 * alpha;
        let (j0_z, j1_z) = bessel_j0_j1(z);

        // 2 * J_1(z) / (z * J_0(z))
        let complex_two = Complex::new(two, T::ZERO);
        let term = complex_two * j1_z / (z * j0_z);
        let bracket = Complex::new(one, T::ZERO) - term;

        // P_hat / (i * rho * omega)
        let coeff = Complex::new(T::ZERO, -p_hat / (rho * omega));

        let phase = omega * t.into_base();
        let exp_iwt = Complex::new(
            <T as FloatElement>::cos(phase),
            <T as FloatElement>::sin(phase),
        );

        // Q = pi * R^2 * Re{ coeff * bracket * exp_iwt }
        VolumetricFlowRate::from_base(pi * r * r * (coeff * bracket * exp_iwt).re)
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use aequitas::systems::si::quantities::{
        DynamicViscosity, Length, MassDensity, PressureGradient, ReciprocalTime, Time,
    };
    use eunomia::assert_relative_eq;
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

            let womersley = WomersleyNumber::new(
                Length::from_base(r),
                ReciprocalTime::from_base(omega),
                MassDensity::from_base(rho),
                DynamicViscosity::from_base(mu),
            );
            let profile = WomersleyProfile::new(
                womersley,
                PressureGradient::from_base(100.0),
            );

            // Wall location is mathematically defined at radial vector xi = r/R = 1.0
            let u_wall = profile.velocity(1.0, Time::from_base(t));

            // Physics Theorem: Velocity uniformly zero at solid boundary (No-Slip Condition)
            assert_relative_eq!(u_wall.into_base(), 0.0, epsilon = 1e-10);
        }
    }
}
