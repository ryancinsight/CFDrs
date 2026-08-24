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
use crate::scalar;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, Length, MassDensity, Pressure, PressureGradient,
    ReciprocalTime, Time, Velocity, VolumetricFlowRate,
};
use cfd_1d::physics::vascular::womersley::{
    WomersleyNumber as ExactWomersleyNumber, WomersleyProfile as ExactWomersleyProfile,
};
use cfd_core::CfdScalar;
use eunomia::FloatElement;
use leto::geometry::Vector3;
use std::f64::consts::PI;

/// Womersley pulsatile flow parameters
#[derive(Debug, Clone)]
pub struct WomersleyFlow<T: CfdScalar> {
    /// Pipe radius
    pub radius: Length<T>,
    /// Fluid density
    pub density: MassDensity<T>,
    /// Dynamic viscosity
    pub viscosity: DynamicViscosity<T>,
    /// Angular frequency of pulsation (rad/s)
    pub omega: ReciprocalTime<T>,
    /// Pressure gradient amplitude
    pub pressure_gradient_amplitude: PressureGradient<T>,
}

impl<T: CfdScalar> WomersleyFlow<T> {
    /// Create new Womersley flow parameters
    pub fn new(
        radius: Length<T>,
        density: MassDensity<T>,
        viscosity: DynamicViscosity<T>,
        omega: ReciprocalTime<T>,
        pressure_gradient_amplitude: PressureGradient<T>,
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
            radius: Length::from_base(<T as FloatElement>::from_f64(4.0e-3)),
            density: MassDensity::from_base(<T as FloatElement>::from_f64(cfd_core::physics::fluid::blood::constants::BLOOD_DENSITY)),
            viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(cfd_core::physics::fluid::blood::constants::INFINITE_SHEAR_VISCOSITY)),
            omega: ReciprocalTime::from_base(<T as FloatElement>::from_f64(2.0 * PI * 1.2)),
            pressure_gradient_amplitude: PressureGradient::from_base(
                <T as FloatElement>::from_f64(1000.0),
            ),
        }
    }

    /// Calculate Womersley number: α = R * sqrt(ω/ν)
    /// This is the key dimensionless parameter for pulsatile flow
    pub fn womersley_number(&self) -> Dimensionless<T> {
        let nu = self.viscosity.into_base() / self.density.into_base();
        Dimensionless::from_base(
            self.radius.into_base() * scalar::sqrt(self.omega.into_base() / nu),
        )
    }

    /// Get characteristic velocity (steady-state Poiseuille maximum velocity)
    pub fn characteristic_velocity(&self) -> Velocity<T> {
        let dpdx = self.pressure_gradient_amplitude.into_base();
        let radius = self.radius.into_base();
        Velocity::from_base(
            -dpdx * radius * radius
                / (<T as FloatElement>::from_f64(4.0) * self.viscosity.into_base()),
        )
    }

    /// Get Stokes layer thickness: δ = sqrt(2ν/ω)
    pub fn stokes_layer_thickness(&self) -> Length<T> {
        let nu = self.viscosity.into_base() / self.density.into_base();
        Length::from_base(scalar::sqrt(
            <T as FloatElement>::from_f64(2.0) * nu / self.omega.into_base(),
        ))
    }

    /// Check if flow is quasi-steady (α < 1)
    pub fn is_quasi_steady(&self) -> bool {
        self.womersley_number().into_base() < scalar::one::<T>()
    }

    /// Check if flow is inertia-dominated (α > 10)
    pub fn is_inertia_dominated(&self) -> bool {
        self.womersley_number().into_base() > <T as FloatElement>::from_f64(10.0)
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
    pub fn velocity(&self, r: Length<T>, t: Time<T>) -> Velocity<T> {
        let xi = scalar::min(
            scalar::abs(r.into_base()) / self.radius.into_base(),
            scalar::one::<T>(),
        );
        self.exact_profile().velocity(xi, t)
    }

    /// Calculate wall shear stress at given time
    /// τ_w = μ * (∂u/∂r)|_{r=R}
    pub fn wall_shear_stress(&self, t: Time<T>) -> Pressure<T> {
        self.exact_profile().wall_shear_stress(t)
    }

    /// Calculate instantaneous flow rate
    pub fn flow_rate(&self, t: Time<T>) -> VolumetricFlowRate<T> {
        self.exact_profile().flow_rate(t)
    }
}

impl<T: CfdScalar> AnalyticalSolution<T> for WomersleyFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, t: T) -> Vector3<T> {
        // Assume pipe axis in x-direction, y is radial coordinate
        // In cylindrical coordinates: r = sqrt(y² + z²)
        // For 2D/axisymmetric: assume z = 0
        let r = scalar::abs(y);

        if r > self.radius.into_base() {
            return Vector3::zeros();
        }

        let u = self.velocity(Length::from_base(r), Time::from_base(t));
        let u = u.into_base();
        Vector3::new(u, scalar::zero::<T>(), scalar::zero::<T>())
    }

    fn pressure(&self, x: T, _y: T, _z: T, t: T) -> T {
        // Pressure gradient drives the flow
        // p(x,t) = -dp/dx * x * cos(ωt) + p0
        let dpdx = self.pressure_gradient_amplitude.into_base();
        -dpdx * x * scalar::cos(self.omega.into_base() * t)
    }

    fn name(&self) -> &'static str {
        "Womersley Pulsatile Flow"
    }

    fn domain_bounds(&self) -> [T; 6] {
        let length = <T as FloatElement>::from_f64(100.0) * self.radius.into_base();
        [
            scalar::zero::<T>(),
            length,
            -self.radius.into_base(),
            self.radius.into_base(),
            -self.radius.into_base(),
            self.radius.into_base(),
        ]
    }

    fn length_scale(&self) -> T {
        self.radius.into_base()
    }

    fn velocity_scale(&self) -> T {
        scalar::abs(self.characteristic_velocity().into_base())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_womersley_number() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();
        let alpha = flow.womersley_number().into_base();

        // Physiological Womersley number should be ~3-5 for medium arteries
        assert!(
            alpha > 2.0 && alpha < 10.0,
            "Womersley number {alpha} is outside physiological range"
        );
    }

    #[test]
    fn test_quasi_steady_limit() {
        // Low frequency flow should be quasi-steady
        let flow = WomersleyFlow::new(
            Length::from_base(1.0e-3),
            MassDensity::from_base(1000.0),
            DynamicViscosity::from_base(0.001),
            ReciprocalTime::from_base(0.1),
            PressureGradient::from_base(100.0),
        );

        assert!(
            flow.is_quasi_steady(),
            "Should be quasi-steady at low frequency"
        );
    }

    #[test]
    fn test_velocity_profile_symmetry() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();
        let t = Time::from_base(0.0);

        // Velocity should be symmetric about centerline
        let u_center = flow.velocity(Length::from_base(0.0), t).into_base();
        let u_off_center1 = flow.velocity(Length::from_base(0.001), t).into_base();
        let u_off_center2 = flow.velocity(Length::from_base(-0.001), t).into_base();

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
            let t_base = time_fraction * 2.0 * std::f64::consts::PI / flow.omega.into_base();
            let t = Time::from_base(t_base);
            let wall_velocity = flow.velocity(flow.radius, t).into_base();
            assert!(
                wall_velocity.abs() < 1e-10,
                "Womersley wall velocity must satisfy no-slip; got {wall_velocity} at t={t_base}",
            );
        }
    }

    #[test]
    fn test_flow_rate_oscillation() {
        let flow = WomersleyFlow::<f64>::physiological_blood_flow();

        // Flow rate should oscillate with time
        let q1 = flow.flow_rate(Time::from_base(0.0)).into_base();
        let q2 = flow
            .flow_rate(Time::from_base(
                std::f64::consts::PI / flow.omega.into_base(),
            ))
            .into_base();

        // Half-period separation reverses the harmonic flow-rate phase.
        assert!(
            q1 * q2 < 0.0 || q2.abs() < 1e-10,
            "Flow rate should reverse during cycle"
        );
    }
}
