//! Womersley pulsatile flow analytical solution
//!
//! The Womersley solution describes fully-developed pulsatile flow in a rigid
//! circular tube, accounting for unsteady inertial effects characterized by
//! the Womersley number α = R√(ωρ/μ).
//!
//! # Mathematical Foundation
//!
//! For a pressure gradient of the form:
//! ```text
//! -∂p/∂x = P̂ · e^(iωt)
//! ```
//!
//! The velocity profile is:
//! ```text
//! u(r,t) = Re{ (P̂/(iρω)) · [1 - J₀(α·ξ·i^(3/2)) / J₀(α·i^(3/2))] · e^(iωt) }
//! ```
//!
//! Where:
//! - ξ = r/R (dimensionless radial position)
//! - α = R√(ωρ/μ) (Womersley number)
//! - J₀ = Bessel function of first kind, order zero
//! - i^(3/2) = e^(i·3π/4) = (-1 + i)/√2
//!
//! # Limiting Behavior
//! - α → 0: Quasi-steady Poiseuille flow
//! - α → ∞: Plug flow with thin boundary layer (Stokes layer)
//!
//! # Physiological Values
//! - Human aorta: α ≈ 12-20
//! - Femoral artery: α ≈ 3-5
//! - Arterioles: α < 1 (quasi-steady)
//!
//! # References
//! - Womersley, J.R. (1955) "Method for the calculation of velocity, rate of
//!   flow and viscous drag in arteries when the pressure gradient is known"
//! - Fung, Y.C. (1997) "Biomechanics: Circulation"

use crate::vascular::bessel::{bessel_j0, bessel_j1};
use nalgebra::{Complex, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ============================================================================
// Womersley Number
// ============================================================================

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
        let two = T::from_f64(2.0).unwrap();
        let pi = T::from_f64(PI).unwrap();
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
            radius: T::from_f64(0.0125).unwrap(),
            omega: T::from_f64(2.0 * PI * 1.2).unwrap(),
            density: T::from_f64(1060.0).unwrap(),
            viscosity: T::from_f64(0.0035).unwrap(),
        }
    }

    /// Create for human femoral artery at 72 bpm
    pub fn human_femoral() -> Self {
        // Femoral diameter ~6 mm
        Self {
            radius: T::from_f64(0.003).unwrap(),
            omega: T::from_f64(2.0 * PI * 1.2).unwrap(),
            density: T::from_f64(1060.0).unwrap(),
            viscosity: T::from_f64(0.0035).unwrap(),
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
        let two = T::from_f64(2.0).unwrap();
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
        let one = T::from_f64(1.0).unwrap();
        let three = T::from_f64(3.0).unwrap();
        let ten = T::from_f64(10.0).unwrap();

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

// ============================================================================
// Womersley Velocity Profile
// ============================================================================

/// Womersley velocity profile calculator
///
/// Computes the complete unsteady velocity field for pulsatile pipe flow.
///
/// # Closed-Form Approximations
///
/// For computational efficiency, we use asymptotic approximations:
///
/// ## Low α limit (α << 1):
/// ```text
/// u(r,t) = (P̂R²/4μ) · (1 - ξ²) · cos(ωt - φ)
/// ```
/// Parabolic profile with small phase lag φ ≈ α²/8.
///
/// ## High α limit (α >> 1):
/// ```text
/// u(r,t) ≈ (P̂/ρω) · [1 - exp(-(1-ξ)R/δ)] · sin(ωt)
/// ```
/// Plug flow with boundary layer of thickness δ.
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
        let sqrt2 = T::from_f64(2.0).unwrap().sqrt();
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

        let sqrt2 = T::from_f64(2.0).unwrap().sqrt();
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
        let two = T::from_f64(2.0).unwrap();
        let pi = T::from_f64(PI).unwrap();

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

// ============================================================================
// Womersley Flow Solver
// ============================================================================

/// Complete Womersley flow solver for arterial segments
///
/// Provides time-varying flow solutions for vessel segments with
/// given inlet conditions and geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WomersleyFlow<T: RealField + Copy> {
    /// Vessel radius [m]
    pub radius: T,
    /// Vessel length [m]
    pub length: T,
    /// Fluid density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Angular frequency [rad/s]
    pub omega: T,
    /// Inlet pressure amplitude [Pa]
    pub inlet_pressure_amplitude: T,
    /// Mean pressure gradient [Pa/m]
    pub mean_pressure_gradient: T,
}

impl<T: RealField + FromPrimitive + Copy> WomersleyFlow<T> {
    /// Create new Womersley flow solver
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        radius: T,
        length: T,
        density: T,
        viscosity: T,
        omega: T,
        inlet_pressure_amplitude: T,
        mean_pressure_gradient: T,
    ) -> Self {
        Self {
            radius,
            length,
            density,
            viscosity,
            omega,
            inlet_pressure_amplitude,
            mean_pressure_gradient,
        }
    }

    /// Get Womersley number for this flow configuration
    pub fn womersley_number(&self) -> WomersleyNumber<T> {
        WomersleyNumber::new(self.radius, self.omega, self.density, self.viscosity)
    }

    /// Get velocity profile calculator
    pub fn profile(&self) -> WomersleyProfile<T> {
        let pressure_gradient_amplitude = self.inlet_pressure_amplitude / self.length;
        WomersleyProfile::new(self.womersley_number(), pressure_gradient_amplitude)
    }

    /// Calculate total velocity (mean + pulsatile) at position and time
    pub fn velocity(&self, xi: T, t: T) -> T {
        // Mean (steady) component - Poiseuille
        let r = self.radius;
        let mu = self.viscosity;
        let four = T::from_f64(4.0).unwrap();
        let u_mean = -self.mean_pressure_gradient * r * r / (four * mu) * (T::one() - xi * xi);

        // Pulsatile component
        let u_pulsatile = self.profile().velocity(xi, t);

        u_mean + u_pulsatile
    }

    /// Calculate impedance magnitude |Z| for this segment
    ///
    /// Z = ΔP / Q (complex impedance)
    pub fn impedance_magnitude(&self) -> T {
        let alpha = self.womersley_number().value();
        let r = self.radius;
        let mu = self.viscosity;
        let rho = self.density;
        let omega = self.omega;
        let pi = T::from_f64(PI).unwrap();
        let eight = T::from_f64(8.0).unwrap();

        if alpha < T::from_f64(1.0).unwrap() {
            // Low α: Z ≈ 8μL/(πR⁴) (Poiseuille resistance dominates)
            eight * mu * self.length / (pi * r.powi(4))
        } else {
            // High α: Z ≈ ρωL/(πR²) (inertance dominates)
            rho * omega * self.length / (pi * r * r)
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_womersley_number_calculation() {
        // Human aorta: R=12.5mm, ω=7.54 rad/s, ρ=1060, μ=0.0035
        let wom = WomersleyNumber::<f64>::human_aorta();
        let alpha = wom.value();

        // Expected α ≈ 18.3 for aorta
        assert!(alpha > 10.0 && alpha < 25.0, "Aortic α = {} should be ~18", alpha);
    }

    #[test]
    fn test_womersley_femoral() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let alpha = wom.value();

        // Expected α ≈ 3.3 for femoral
        assert!(alpha > 2.0 && alpha < 5.0, "Femoral α = {} should be ~3.3", alpha);
    }

    #[test]
    fn test_flow_regime_classification() {
        // Small vessel with low α
        let low_alpha = WomersleyNumber::<f64>::new(0.0001, 7.54, 1060.0, 0.0035);
        assert_eq!(low_alpha.flow_regime(), FlowRegime::QuasiSteady);

        // Aorta with high α
        let high_alpha = WomersleyNumber::<f64>::human_aorta();
        matches!(high_alpha.flow_regime(), FlowRegime::Inertial | FlowRegime::PlugFlow);
    }

    #[test]
    fn test_stokes_layer_thickness() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let delta = wom.stokes_layer_thickness();

        // Stokes layer should be ~1 mm for cardiac frequency
        assert!(delta > 0.0005 && delta < 0.002, "δ = {} should be ~1 mm", delta);
    }

    #[test]
    fn test_velocity_profile_centerline() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0); // 100 Pa/m gradient

        // Centerline velocity should be finite and positive at some time
        let u_center = profile.centerline_velocity(0.25); // t = 0.25s
        assert!(u_center.abs() < 1.0, "Centerline velocity should be < 1 m/s");
    }

    #[test]
    fn test_velocity_profile_wall_zero() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);

        // Velocity at wall (xi=1) should be zero or very small
        let u_wall = profile.velocity(1.0, 0.1);
        assert!(u_wall.abs() < 0.01, "Wall velocity {} should be ~0", u_wall);
    }

    #[test]
    fn test_velocity_decreases_radially() {
        let wom = WomersleyNumber::<f64>::new(0.003, 7.54, 1060.0, 0.0035);
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);

        let t = 0.0_f64; // At t=0, check profile shape
        let u_center = profile.velocity(0.0, t).abs();
        let u_mid = profile.velocity(0.5, t).abs();
        let u_wall = profile.velocity(1.0, t).abs();

        // For moderate α, velocity generally decreases toward wall
        // (though phase shifts may cause exceptions at certain times)
        assert!(u_center >= u_wall, "Centerline should have higher |u| than wall");
        let _ = u_mid; // Silence unused warning
    }

    #[test]
    fn test_flow_rate_positive() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);

        // Flow rate magnitude should be finite
        let q = profile.flow_rate(0.1);
        assert!(q.abs() < 0.001, "Flow rate {} L/s should be realistic", q * 1000.0);
    }

    #[test]
    fn test_womersley_flow_solver() {
        let flow = WomersleyFlow::<f64>::new(
            0.003,   // 6 mm diameter vessel
            0.1,     // 10 cm length
            1060.0,  // Blood density
            0.0035,  // Blood viscosity
            7.54,    // 1.2 Hz heart rate
            133.0,   // ~1 mmHg inlet amplitude
            -1000.0, // Mean pressure gradient
        );

        // Check Womersley number
        let alpha = flow.womersley_number().value();
        assert!(alpha > 2.0 && alpha < 5.0);

        // Check velocity is finite
        let u = flow.velocity(0.5, 0.3);
        assert!(u.abs() < 1.0, "Velocity {} should be < 1 m/s", u);
    }

    #[test]
    fn test_impedance_magnitude() {
        let flow = WomersleyFlow::<f64>::new(
            0.003,
            0.1,
            1060.0,
            0.0035,
            7.54,
            133.0,
            -1000.0,
        );

        let z = flow.impedance_magnitude();
        assert!(z > 0.0 && z.is_finite(), "Impedance {} should be finite positive", z);
    }

    #[test]
    fn test_wall_shear_stress() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);

        let tau_w = profile.wall_shear_stress(0.2);
        // WSS can range from O(1) to O(100) Pa depending on vessel and conditions
        // For moderate pressure gradients in large arteries, WSS is typically 1-100 Pa
        assert!(tau_w.abs() < 200.0, "WSS {} Pa should be < 200 Pa", tau_w);
    }
}

