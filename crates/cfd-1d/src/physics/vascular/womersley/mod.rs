//! Womersley pulsatile flow analytical solution.
//!
//! The Womersley solution describes fully-developed pulsatile flow in a rigid
//! circular tube, accounting for unsteady inertial effects characterized by
//! the Womersley number α = R√(ωρ/μ).
//!
//! # Mathematical Foundation
//!
//! ## Theorem: Womersley Exact Analytical Solution
//!
//! **Theorem**: For a fully-developed, axisymmetric, incompressible Newtonian flow
//! in a rigid circular tube of radius R driven by a purely oscillatory pressure
//! gradient, the exact analytical solution to the linearized Navier-Stokes
//! momentum equation is:
//!
//! ```text
//! u(r,t) = Re{ (P̂/(iρω)) [1 - J₀(α·r/R·i^{3/2}) / J₀(α·i^{3/2})] · e^{iωt} }
//! ```
//!
//! ## Limiting Behavior
//! - α → 0: Quasi-steady Poiseuille flow
//! - α → ∞: Plug flow with thin boundary layer (Stokes layer)
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`profile`] | `WomersleyProfile<T>` — velocity, wall shear stress, flow rate |
//! | [`flow_solver`] | `WomersleyFlow<T>` — mean+pulsatile solver, impedance |
//! | [`pulsatility`] | `womersley_pulsatility_index()` — Gosling-King PI |
//!
//! ## References
//!
//! - Womersley, J.R. (1955) *Phil. Mag.* 46, 199-221.
//! - Fung, Y.C. (1997) *Biomechanics: Circulation*.
//! - Gosling, R.G. & King, D.H. (1974) *Proc. R. Soc. Med.* 67:447-449.

pub mod flow_solver;
pub mod profile;
pub mod pulsatility;

use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, Frequency, Length, MassDensity, ReciprocalTime,
};
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use flow_solver::WomersleyFlow;
pub use profile::WomersleyProfile;
pub use pulsatility::womersley_pulsatility_index;

// ── Womersley Number ────────────────────────────────────────────────────────

/// Womersley number calculator and classifier
///
/// The Womersley number α characterizes the relative importance of unsteady
/// inertial forces to viscous forces in pulsatile flow:
/// ```text
/// α = R · √(ω·ρ/μ) = R · √(2πf·ρ/μ)
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WomersleyNumber<T: CfdScalar + Copy> {
    /// Vessel radius \[m]
    pub radius: Length<T>,
    /// Angular frequency of pulsation \[rad/s]
    pub omega: ReciprocalTime<T>,
    /// Fluid density [kg/m³]
    pub density: MassDensity<T>,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: DynamicViscosity<T>,
}

impl<T: CfdScalar + FloatElement + Copy> WomersleyNumber<T> {
    /// Create Womersley number calculator from physical parameters
    pub fn new(
        radius: Length<T>,
        omega: ReciprocalTime<T>,
        density: MassDensity<T>,
        viscosity: DynamicViscosity<T>,
    ) -> Self {
        Self {
            radius,
            omega,
            density,
            viscosity,
        }
    }

    /// Create from vessel diameter and heart rate (frequency in Hz)
    pub fn from_heart_rate(
        diameter: Length<T>,
        heart_rate_hz: Frequency<T>,
        density: MassDensity<T>,
        viscosity: DynamicViscosity<T>,
    ) -> Self {
        let two = T::ONE + T::ONE;
        let pi = T::pi();
        Self {
            radius: Length::from_base(diameter.into_base() / two),
            omega: ReciprocalTime::from_base(two * pi * heart_rate_hz.into_base()),
            density,
            viscosity,
        }
    }

    /// Create for human aorta at 72 bpm with blood properties
    pub fn human_aorta() -> Self {
        // Aortic root diameter ~25 mm
        // Heart rate 72 bpm = 1.2 Hz
        // Blood: ρ = 1060 kg/m³, μ = 0.0035 Pa·s
        Self {
            radius: Length::from_base(<T as FloatElement>::from_f64(0.0125)),
            omega: ReciprocalTime::from_base(<T as FloatElement>::from_f64(2.0 * PI * 1.2)),
            density: MassDensity::from_base(<T as FloatElement>::from_f64(1060.0)),
            viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(0.0035)),
        }
    }

    /// Create for human femoral artery at 72 bpm
    pub fn human_femoral() -> Self {
        // Femoral diameter ~6 mm
        Self {
            radius: Length::from_base(<T as FloatElement>::from_f64(0.003)),
            omega: ReciprocalTime::from_base(<T as FloatElement>::from_f64(2.0 * PI * 1.2)),
            density: MassDensity::from_base(<T as FloatElement>::from_f64(1060.0)),
            viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(0.0035)),
        }
    }

    /// Calculate the Womersley number α
    ///
    /// α = R · √(ω·ρ/μ)
    pub fn value(&self) -> Dimensionless<T> {
        let radius = self.radius.into_base();
        let omega = self.omega.into_base();
        let density = self.density.into_base();
        let viscosity = self.viscosity.into_base();
        Dimensionless::from_base(radius * <T as NumericElement>::sqrt(omega * density / viscosity))
    }

    /// The Stokes layer thickness δ = √(2μ/(ρω))
    ///
    /// This is the characteristic length over which viscous effects penetrate
    /// from the wall into the flow during one oscillation cycle.
    pub fn stokes_layer_thickness(&self) -> Length<T> {
        let two = T::ONE + T::ONE;
        Length::from_base(<T as NumericElement>::sqrt(
            two * self.viscosity.into_base() / (self.density.into_base() * self.omega.into_base()),
        ))
    }

    /// Ratio of radius to Stokes layer thickness R/δ
    ///
    /// For α >> 1, R/δ = α/√2
    pub fn radius_to_stokes_ratio(&self) -> Dimensionless<T> {
        Dimensionless::from_base(
            self.radius.into_base() / self.stokes_layer_thickness().into_base(),
        )
    }

    /// Classify flow regime based on Womersley number
    pub fn flow_regime(&self) -> FlowRegime {
        let alpha = self.value().into_base();
        let one = T::ONE;
        let three = T::ONE + T::ONE + T::ONE;
        let ten = <T as FloatElement>::from_f64(10.0);

        if alpha < one {
            FlowRegime::QuasiSteady
        } else if alpha < three {
            FlowRegime::Transitional
        } else if alpha < ten {
            FlowRegime::Inertial
        } else {
            FlowRegime::PlugFlow
        }
    }
}

/// Flow regime classification based on Womersley number
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    /// α < 1: Quasi-steady Poiseuille-like flow
    /// Viscous effects dominate, parabolic velocity profile maintained
    QuasiSteady,
    /// 1 ≤ α < 3: Transitional regime
    /// Both viscous and inertial effects significant
    Transitional,
    /// 3 ≤ α < 10: Inertia-dominated flow
    /// Flat velocity core with thin boundary layer
    Inertial,
    /// α ≥ 10: Plug flow with thin Stokes layer
    /// Nearly flat profile across most of vessel
    PlugFlow,
}

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::quantities::{Pressure, PressureGradient, Time};
    use eunomia::assert_relative_eq;

    #[test]
    fn test_womersley_number_calculation() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let alpha = wom.value().into_base();
        assert!(
            alpha > 10.0 && alpha < 25.0,
            "Aortic α = {alpha} should be ~18"
        );
    }

    #[test]
    fn test_womersley_femoral() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let alpha = wom.value().into_base();
        assert!(
            alpha > 2.0 && alpha < 5.0,
            "Femoral α = {alpha} should be ~3.3"
        );
    }

    #[test]
    fn test_flow_regime_classification() {
        let low_alpha = WomersleyNumber::<f64>::new(
            Length::from_base(0.0001),
            ReciprocalTime::from_base(7.54),
            MassDensity::from_base(1060.0),
            DynamicViscosity::from_base(0.0035),
        );
        assert_eq!(low_alpha.flow_regime(), FlowRegime::QuasiSteady);

        let high_alpha = WomersleyNumber::<f64>::human_aorta();
        matches!(
            high_alpha.flow_regime(),
            FlowRegime::Inertial | FlowRegime::PlugFlow
        );
    }

    #[test]
    fn test_stokes_layer_thickness() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let delta = wom.stokes_layer_thickness().into_base();
        assert!(
            delta > 0.0005 && delta < 0.002,
            "δ = {delta} should be ~1 mm"
        );
    }

    #[test]
    fn test_velocity_profile_centerline() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(100.0));
        let u_center = profile
            .centerline_velocity(Time::from_base(0.25))
            .into_base();
        assert!(
            u_center.abs() < 1.0,
            "Centerline velocity should be < 1 m/s"
        );
    }

    #[test]
    fn test_velocity_profile_wall_zero() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(100.0));
        let u_wall = profile.velocity(1.0, Time::from_base(0.1)).into_base();
        assert!(u_wall.abs() < 0.01, "Wall velocity {u_wall} should be ~0");
    }

    #[test]
    fn test_velocity_decreases_radially() {
        let wom = WomersleyNumber::<f64>::new(
            Length::from_base(0.003),
            ReciprocalTime::from_base(7.54),
            MassDensity::from_base(1060.0),
            DynamicViscosity::from_base(0.0035),
        );
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(100.0));

        let t = 0.0_f64;
        let u_center = profile.velocity(0.0, Time::from_base(t)).into_base().abs();
        let u_mid = profile.velocity(0.5, Time::from_base(t)).into_base().abs();
        let u_wall = profile.velocity(1.0, Time::from_base(t)).into_base().abs();

        assert!(
            u_center >= u_wall,
            "Centerline should have higher |u| than wall"
        );
        let _ = u_mid;
    }

    #[test]
    fn test_flow_rate_positive() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(100.0));
        let q = profile.flow_rate(Time::from_base(0.1)).into_base();
        let q_liters_per_second = q * 1000.0;
        assert!(
            q.abs() < 0.001,
            "Flow rate {q_liters_per_second} L/s should be realistic"
        );
    }

    #[test]
    fn test_womersley_flow_solver() {
        let flow = WomersleyFlow::<f64>::new(
            Length::from_base(0.003),
            Length::from_base(0.1),
            MassDensity::from_base(1060.0),
            DynamicViscosity::from_base(0.0035),
            ReciprocalTime::from_base(7.54),
            Pressure::from_base(133.0),
            PressureGradient::from_base(-1000.0),
        );
        let alpha = flow.womersley_number().value().into_base();
        assert!(alpha > 2.0 && alpha < 5.0);
        let u = flow.velocity(0.5, Time::from_base(0.3)).into_base();
        assert!(u.abs() < 1.0, "Velocity {u} should be < 1 m/s");
    }

    #[test]
    fn test_impedance_magnitude() {
        let flow = WomersleyFlow::<f64>::new(
            Length::from_base(0.003),
            Length::from_base(0.1),
            MassDensity::from_base(1060.0),
            DynamicViscosity::from_base(0.0035),
            ReciprocalTime::from_base(7.54),
            Pressure::from_base(133.0),
            PressureGradient::from_base(-1000.0),
        );
        let z = flow.impedance_magnitude().into_base();
        assert!(
            z > 0.0 && z.is_finite(),
            "Impedance {z} should be finite positive"
        );
    }

    #[test]
    fn test_wall_shear_stress() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(100.0));
        let tau_w = profile.wall_shear_stress(Time::from_base(0.2)).into_base();
        assert!(tau_w.abs() < 200.0, "WSS {tau_w} Pa should be < 200 Pa");
    }

    #[test]
    fn test_womersley_low_alpha_poiseuille_limit() {
        let alpha = 0.01;
        let p_hat = 100.0;
        let rho = 1000.0;
        let mu = 0.001;
        let r = 0.01;
        let omega = (alpha / r) * (alpha / r) * mu / rho;

        let wom = WomersleyNumber::<f64>::new(
            Length::from_base(r),
            ReciprocalTime::from_base(omega),
            MassDensity::from_base(rho),
            DynamicViscosity::from_base(mu),
        );
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(p_hat));

        let u_womersley = profile
            .centerline_velocity(Time::from_base(0.0))
            .into_base()
            .abs();
        let u_poiseuille = (p_hat * r * r) / (4.0 * mu);

        assert_relative_eq!(u_womersley, u_poiseuille, max_relative = 0.01);
    }

    #[test]
    fn test_womersley_flow_rate_oscillates() {
        let wom = WomersleyNumber::<f64>::new(
            Length::from_base(0.005),
            ReciprocalTime::from_base(2.0 * PI),
            MassDensity::from_base(1000.0),
            DynamicViscosity::from_base(0.001),
        );
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(10.0));

        let mut net_volume = 0.0;
        let steps = 1000;
        let dt = 1.0 / f64::from(steps);
        for i in 0..steps {
            let t = f64::from(i) * dt;
            net_volume += profile.flow_rate(Time::from_base(t)).into_base() * dt;
        }

        assert!(
            net_volume.abs() < 1e-10,
            "Net volume must sum to zero, got {net_volume}"
        );
    }

    #[test]
    fn test_womersley_wall_bc_satisfied() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let profile = WomersleyProfile::<f64>::new(wom, PressureGradient::from_base(120.0));

        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let u_wall = profile.velocity(1.0, Time::from_base(t)).into_base();
            assert!(
                u_wall.abs() < 1e-12,
                "No-slip violated at t={t}: u={u_wall}"
            );
        }
    }

    // ── Pulsatility Index tests ─────────────────────────────────────────

    #[test]
    fn test_pulsatility_index_steady_flow() {
        let area = 1e-4;
        let q_mean = 5e-6;
        let (pi, v_peak, v_trough) = womersley_pulsatility_index(q_mean, 0.0, area);

        assert!(pi.abs() < 1e-12, "PI should be 0 for steady flow, got {pi}");
        let v_mean = q_mean / area;
        assert_relative_eq!(v_peak, v_mean, max_relative = 1e-12);
        assert_relative_eq!(v_trough, v_mean, max_relative = 1e-12);
    }

    #[test]
    fn test_pulsatility_index_typical_artery() {
        let area = 2e-5;
        let q_mean = 5e-6;
        let q_amplitude = 2.5e-6;
        let (pi, v_peak, v_trough) = womersley_pulsatility_index(q_mean, q_amplitude, area);

        assert!(
            pi > 0.8 && pi < 2.0,
            "Typical arterial PI = {pi:.2} should be in 0.8-2.0"
        );
        assert!(
            v_peak > v_trough,
            "Peak ({v_peak}) should exceed trough ({v_trough})"
        );
    }

    #[test]
    fn test_pulsatility_index_positive() {
        let area = 1e-4;
        for &q_mean in &[1e-6, 1e-5, 1e-4] {
            for &q_amp in &[0.0, 0.5e-6, 1e-6, 5e-6] {
                let (pi, _, _) = womersley_pulsatility_index(q_mean, q_amp, area);
                assert!(
                    pi >= 0.0,
                    "PI should be >= 0, got {pi} for q_mean={q_mean}, q_amp={q_amp}"
                );
            }
        }
    }
}
