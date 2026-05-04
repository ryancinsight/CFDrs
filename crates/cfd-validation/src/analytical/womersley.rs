//! Womersley pulsatile flow - oscillatory flow in a circular pipe
//!
//! This is a classic analytical solution for pulsatile blood flow in arteries.
//! The Womersley number (α) characterizes the ratio of unsteady to viscous forces.
//!
//! # References
//! - Womersley, J.R. (1955). "Method for the calculation of velocity, rate of flow
//!   and viscous drag in arteries when the pressure gradient is known"
//!   Journal of Physiology, 127(3):553-563.
//! - Zamir, M. (2000). "The Physics of Pulsatile Flow", Springer.

use super::AnalyticalSolution;
use cfd_1d::physics::vascular::womersley::{
    WomersleyNumber as ExactWomersleyNumber, WomersleyProfile as ExactWomersleyProfile,
};
use cfd_core::conversion::SafeFromF64;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use std::f64::consts::PI;

/// Womersley pulsatile flow parameters
#[derive(Debug, Clone)]
pub struct WomersleyFlow<T: RealField + Copy> {
    /// Pipe radius
    pub radius: T,
    /// Fluid density
    pub density: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Angular frequency of pulsation (rad/s)
    pub omega: T,
    /// Pressure gradient amplitude
    pub pressure_gradient_amplitude: T,
}

impl<T: RealField + Copy + FromPrimitive> WomersleyFlow<T> {
    /// Create new Womersley flow parameters
    pub fn new(
        radius: T,
        density: T,
        viscosity: T,
        omega: T,
        pressure_gradient_amplitude: T,
    ) -> Self {
        Self {
            radius,
            density,
            viscosity,
            omega,
            pressure_gradient_amplitude,
        }
    }

    /// Create physiological blood flow parameters
    /// Typical values for human carotid artery
    pub fn physiological_blood_flow() -> Self {
        Self {
            radius: T::from_f64_or_one(4.0e-3),                      // 4 mm
            density: T::from_f64_or_one(1060.0),                     // Blood density
            viscosity: T::from_f64_or_one(0.0035),                   // Pa·s (apparent)
            omega: T::from_f64_or_one(2.0 * PI * 1.2),               // ~72 bpm
            pressure_gradient_amplitude: T::from_f64_or_one(1000.0), // Pa/m
        }
    }

    /// Calculate Womersley number: α = R * sqrt(ω/ν)
    /// This is the key dimensionless parameter for pulsatile flow
    pub fn womersley_number(&self) -> T {
        let nu = self.viscosity / self.density; // Kinematic viscosity
        self.radius * (self.omega / nu).sqrt()
    }

    /// Get characteristic velocity (steady-state Poiseuille maximum velocity)
    pub fn characteristic_velocity(&self) -> T {
        let dpdx = self.pressure_gradient_amplitude;
        -dpdx * self.radius * self.radius / (T::from_f64_or_one(4.0) * self.viscosity)
    }

    /// Get Stokes layer thickness: δ = sqrt(2ν/ω)
    pub fn stokes_layer_thickness(&self) -> T {
        let nu = self.viscosity / self.density;
        (T::from_f64_or_one(2.0) * nu / self.omega).sqrt()
    }

    /// Check if flow is quasi-steady (α < 1)
    pub fn is_quasi_steady(&self) -> bool {
        self.womersley_number() < T::one()
    }

    /// Check if flow is inertia-dominated (α > 10)
    pub fn is_inertia_dominated(&self) -> bool {
        self.womersley_number() > T::from_f64_or_one(10.0)
    }

    /// Create the canonical exact Womersley profile evaluator.
    fn exact_profile(&self) -> ExactWomersleyProfile<T> {
        ExactWomersleyProfile::new(
            ExactWomersleyNumber::new(self.radius, self.omega, self.density, self.viscosity),
            self.pressure_gradient_amplitude,
        )
    }

    /// Calculate analytical velocity profile at given radius and time.
    ///
    /// # Theorem — Womersley No-Slip Profile
    ///
    /// The axial velocity is evaluated by the canonical complex-Bessel
    /// Womersley solution
    /// `u(r,t) = Re{ P/(i rho omega) [1 - J0(i^(3/2) alpha r/R) /
    /// J0(i^(3/2) alpha)] exp(i omega t) }`.
    ///
    /// **Proof sketch**: Substitution of the harmonic ansatz into the
    /// axisymmetric unsteady Stokes equation yields Bessel's equation in the
    /// radial coordinate. The ratio term enforces `u(R,t)=0`; finite `J0(0)`
    /// gives a bounded centerline value. This wrapper delegates to the
    /// single-source `cfd-1d` implementation that evaluates the same closed
    /// form with shared `J0/J1` recurrence. ∎
    pub fn velocity(&self, r: T, t: T) -> T {
        let xi = (r.abs() / self.radius).min(T::one());
        self.exact_profile().velocity(xi, t)
    }

    /// Calculate wall shear stress at given time
    /// τ_w = μ * (∂u/∂r)|_{r=R}
    pub fn wall_shear_stress(&self, t: T) -> T {
        self.exact_profile().wall_shear_stress(t)
    }

    /// Calculate instantaneous flow rate
    pub fn flow_rate(&self, t: T) -> T {
        self.exact_profile().flow_rate(t)
    }
}

impl<T: RealField + Copy + FromPrimitive> AnalyticalSolution<T> for WomersleyFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, t: T) -> Vector3<T> {
        // Assume pipe axis in x-direction, y is radial coordinate
        // In cylindrical coordinates: r = sqrt(y² + z²)
        // For 2D/axisymmetric: assume z = 0
        let r = y.abs();

        if r > self.radius {
            return Vector3::zeros();
        }

        let u = self.velocity(r, t);
        Vector3::new(u, T::zero(), T::zero())
    }

    fn pressure(&self, x: T, _y: T, _z: T, t: T) -> T {
        // Pressure gradient drives the flow
        // p(x,t) = -dp/dx * x * cos(ωt) + p0
        let dpdx = self.pressure_gradient_amplitude;
        -dpdx * x * (self.omega * t).cos()
    }

    fn name(&self) -> &'static str {
        "Womersley Pulsatile Flow"
    }

    fn domain_bounds(&self) -> [T; 6] {
        let length = T::from_f64_or_one(100.0) * self.radius; // 100R long pipe
        [
            T::zero(),
            length, // x: [0, L]
            -self.radius,
            self.radius, // y: [-R, R]
            -self.radius,
            self.radius, // z: [-R, R]
        ]
    }

    fn length_scale(&self) -> T {
        self.radius
    }

    fn velocity_scale(&self) -> T {
        self.characteristic_velocity().abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_womersley_number() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();
        let alpha = flow.womersley_number();

        // Physiological Womersley number should be ~3-5 for medium arteries
        assert!(
            alpha > 2.0 && alpha < 10.0,
            "Womersley number {} is outside physiological range",
            alpha
        );
    }

    #[test]
    fn test_quasi_steady_limit() {
        // Low frequency flow should be quasi-steady
        let flow = WomersleyFlow::new(
            1.0e-3, // 1 mm radius
            1000.0, // Water density
            0.001,  // Water viscosity
            0.1,    // Very low frequency
            100.0,
        );

        assert!(
            flow.is_quasi_steady(),
            "Should be quasi-steady at low frequency"
        );
    }

    #[test]
    fn test_velocity_profile_symmetry() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();
        let t = 0.0;

        // Velocity should be symmetric about centerline
        let u_center = flow.velocity(0.0, t);
        let u_off_center1 = flow.velocity(0.001, t);
        let u_off_center2 = flow.velocity(-0.001, t);

        assert!(
            (u_off_center1 - u_off_center2).abs() < 1e-10,
            "Velocity profile should be symmetric"
        );
        assert!(
            u_center.abs() >= u_off_center1.abs(),
            "Center velocity should be maximum"
        );
    }

    #[test]
    fn test_exact_womersley_no_slip_wall_condition() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();

        for time_fraction in [0.0, 0.125, 0.25, 0.5, 0.875] {
            let t = time_fraction * 2.0 * std::f64::consts::PI / flow.omega;
            let wall_velocity = flow.velocity(flow.radius, t);
            assert!(
                wall_velocity.abs() < 1e-10,
                "Womersley wall velocity must satisfy no-slip; got {wall_velocity} at t={t}"
            );
        }
    }

    #[test]
    fn test_flow_rate_oscillation() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();

        // Flow rate should oscillate with time
        let q1 = flow.flow_rate(0.0);
        let q2 = flow.flow_rate(std::f64::consts::PI / flow.omega);

        // Half-period separation reverses the harmonic flow-rate phase.
        assert!(
            q1 * q2 < 0.0 || q2.abs() < 1e-10,
            "Flow rate should reverse during cycle"
        );
    }
}
