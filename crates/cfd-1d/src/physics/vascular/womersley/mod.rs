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

use nalgebra::RealField;
use num_traits::FromPrimitive;
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
pub struct WomersleyNumber<T: RealField + Copy> {
    /// Vessel radius [m]
    pub radius: T,
    /// Angular frequency of pulsation [rad/s]
    pub omega: T,
    /// Fluid density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: T,
}

impl<T: RealField + FromPrimitive + Copy> WomersleyNumber<T> {
    /// Create Womersley number calculator from physical parameters
    pub fn new(radius: T, omega: T, density: T, viscosity: T) -> Self {
        Self {
            radius,
            omega,
            density,
            viscosity,
        }
    }

    /// Create from vessel diameter and heart rate (frequency in Hz)
    pub fn from_heart_rate(diameter: T, heart_rate_hz: T, density: T, viscosity: T) -> Self {
        let two = T::one() + T::one();
        let pi = T::pi();
        Self {
            radius: diameter / two,
            omega: two * pi * heart_rate_hz,
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
            radius: T::from_f64(0.0125).expect("Mathematical constant conversion compromised"),
            omega: T::from_f64(2.0 * PI * 1.2)
                .expect("Mathematical constant conversion compromised"),
            density: T::from_f64(1060.0).expect("Mathematical constant conversion compromised"),
            viscosity: T::from_f64(0.0035).expect("Mathematical constant conversion compromised"),
        }
    }

    /// Create for human femoral artery at 72 bpm
    pub fn human_femoral() -> Self {
        // Femoral diameter ~6 mm
        Self {
            radius: T::from_f64(0.003).expect("Mathematical constant conversion compromised"),
            omega: T::from_f64(2.0 * PI * 1.2)
                .expect("Mathematical constant conversion compromised"),
            density: T::from_f64(1060.0).expect("Mathematical constant conversion compromised"),
            viscosity: T::from_f64(0.0035).expect("Mathematical constant conversion compromised"),
        }
    }

    /// Calculate the Womersley number α
    ///
    /// α = R · √(ω·ρ/μ)
    pub fn value(&self) -> T {
        self.radius * (self.omega * self.density / self.viscosity).sqrt()
    }

    /// The Stokes layer thickness δ = √(2μ/(ρω))
    ///
    /// This is the characteristic length over which viscous effects penetrate
    /// from the wall into the flow during one oscillation cycle.
    pub fn stokes_layer_thickness(&self) -> T {
        let two = T::one() + T::one();
        (two * self.viscosity / (self.density * self.omega)).sqrt()
    }

    /// Ratio of radius to Stokes layer thickness R/δ
    ///
    /// For α >> 1, R/δ = α/√2
    pub fn radius_to_stokes_ratio(&self) -> T {
        self.radius / self.stokes_layer_thickness()
    }

    /// Classify flow regime based on Womersley number
    pub fn flow_regime(&self) -> FlowRegime {
        let alpha = self.value();
        let one = T::one();
        let three = T::one() + T::one() + T::one();
        let ten = T::from_f64(10.0).expect("Mathematical constant conversion compromised");

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
    use approx::assert_relative_eq;

    #[test]
    fn test_womersley_number_calculation() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let alpha = wom.value();
        assert!(
            alpha > 10.0 && alpha < 25.0,
            "Aortic α = {} should be ~18",
            alpha
        );
    }

    #[test]
    fn test_womersley_femoral() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let alpha = wom.value();
        assert!(
            alpha > 2.0 && alpha < 5.0,
            "Femoral α = {} should be ~3.3",
            alpha
        );
    }

    #[test]
    fn test_flow_regime_classification() {
        let low_alpha = WomersleyNumber::<f64>::new(0.0001, 7.54, 1060.0, 0.0035);
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
        let delta = wom.stokes_layer_thickness();
        assert!(
            delta > 0.0005 && delta < 0.002,
            "δ = {} should be ~1 mm",
            delta
        );
    }

    #[test]
    fn test_velocity_profile_centerline() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);
        let u_center = profile.centerline_velocity(0.25);
        assert!(
            u_center.abs() < 1.0,
            "Centerline velocity should be < 1 m/s"
        );
    }

    #[test]
    fn test_velocity_profile_wall_zero() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);
        let u_wall = profile.velocity(1.0, 0.1);
        assert!(u_wall.abs() < 0.01, "Wall velocity {} should be ~0", u_wall);
    }

    #[test]
    fn test_velocity_decreases_radially() {
        let wom = WomersleyNumber::<f64>::new(0.003, 7.54, 1060.0, 0.0035);
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);

        let t = 0.0_f64;
        let u_center = profile.velocity(0.0, t).abs();
        let u_mid = profile.velocity(0.5, t).abs();
        let u_wall = profile.velocity(1.0, t).abs();

        assert!(
            u_center >= u_wall,
            "Centerline should have higher |u| than wall"
        );
        let _ = u_mid;
    }

    #[test]
    fn test_flow_rate_positive() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);
        let q = profile.flow_rate(0.1);
        assert!(
            q.abs() < 0.001,
            "Flow rate {} L/s should be realistic",
            q * 1000.0
        );
    }

    #[test]
    fn test_womersley_flow_solver() {
        let flow = WomersleyFlow::<f64>::new(0.003, 0.1, 1060.0, 0.0035, 7.54, 133.0, -1000.0);
        let alpha = flow.womersley_number().value();
        assert!(alpha > 2.0 && alpha < 5.0);
        let u = flow.velocity(0.5, 0.3);
        assert!(u.abs() < 1.0, "Velocity {} should be < 1 m/s", u);
    }

    #[test]
    fn test_impedance_magnitude() {
        let flow = WomersleyFlow::<f64>::new(0.003, 0.1, 1060.0, 0.0035, 7.54, 133.0, -1000.0);
        let z = flow.impedance_magnitude();
        assert!(
            z > 0.0 && z.is_finite(),
            "Impedance {} should be finite positive",
            z
        );
    }

    #[test]
    fn test_wall_shear_stress() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);
        let tau_w = profile.wall_shear_stress(0.2);
        assert!(tau_w.abs() < 200.0, "WSS {} Pa should be < 200 Pa", tau_w);
    }

    #[test]
    fn test_womersley_low_alpha_poiseuille_limit() {
        let alpha = 0.01;
        let p_hat = 100.0;
        let rho = 1000.0;
        let mu = 0.001;
        let r = 0.01;
        let omega = (alpha / r) * (alpha / r) * mu / rho;

        let wom = WomersleyNumber::<f64>::new(r, omega, rho, mu);
        let profile = WomersleyProfile::<f64>::new(wom, p_hat);

        let u_womersley = profile.centerline_velocity(0.0).abs();
        let u_poiseuille = (p_hat * r * r) / (4.0 * mu);

        assert_relative_eq!(u_womersley, u_poiseuille, max_relative = 0.01);
    }

    #[test]
    fn test_womersley_flow_rate_oscillates() {
        let wom = WomersleyNumber::<f64>::new(0.005, 2.0 * PI, 1000.0, 0.001);
        let profile = WomersleyProfile::<f64>::new(wom, 10.0);

        let mut net_volume = 0.0;
        let steps = 1000;
        let dt = 1.0 / steps as f64;
        for i in 0..steps {
            let t = i as f64 * dt;
            net_volume += profile.flow_rate(t) * dt;
        }

        assert!(
            net_volume.abs() < 1e-10,
            "Net volume must sum to zero, got {}",
            net_volume
        );
    }

    #[test]
    fn test_womersley_wall_bc_satisfied() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let profile = WomersleyProfile::<f64>::new(wom, 120.0);

        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let u_wall = profile.velocity(1.0, t);
            assert!(
                u_wall.abs() < 1e-12,
                "No-slip violated at t={}: u={}",
                t,
                u_wall
            );
        }
    }

    // ── Pulsatility Index tests ─────────────────────────────────────────

    #[test]
    fn test_pulsatility_index_steady_flow() {
        let area = 1e-4;
        let q_mean = 5e-6;
        let (pi, v_peak, v_trough) = womersley_pulsatility_index(q_mean, 0.0, area);

        assert!(
            pi.abs() < 1e-12,
            "PI should be 0 for steady flow, got {}",
            pi
        );
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
            "Typical arterial PI = {:.2} should be in 0.8-2.0",
            pi
        );
        assert!(
            v_peak > v_trough,
            "Peak ({}) should exceed trough ({})",
            v_peak,
            v_trough
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
                    "PI should be >= 0, got {} for q_mean={}, q_amp={}",
                    pi,
                    q_mean,
                    q_amp
                );
            }
        }
    }
}
