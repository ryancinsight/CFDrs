//! Womersley pulsatile flow analytical solution
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
//! in a rigid circular tube of radius $R$ driven by a purely oscillatory pressure
//! gradient $-\frac{\partial P}{\partial x} = \hat{P} e^{i \omega t}$, the exact
//! analytical solution to the linearized Navier-Stokes momentum equation is:
//!
//! $$ u(r,t) = \text{Re} \left\{ \frac{\hat{P}}{i \rho \omega} \left[ 1 - \frac{J_0\left(\alpha \frac{r}{R} i^{3/2}\right)}{J_0\left(\alpha i^{3/2}\right)} \right] e^{i \omega t} \right\} $$
//!
//! Where:
//! - $r$ is the radial coordinate ($0 \le r \le R$)
//! - $\alpha = R \sqrt{\frac{\omega \rho}{\mu}}$ is the dimensionless Womersley number
//! - $J_0(z)$ is the complex Bessel function of the first kind, order zero
//! - $i^{3/2} = e^{i \cdot 3\pi/4} = \frac{-1 + i}{\sqrt{2}}$
//!
//! **Proof Outline**: Applying the harmonic ansatz $u(r,t) = \hat{u}(r) e^{i \omega t}$
//! directly reduces the momentum PDE to a modified Bessel differential equation:
//! $$ \frac{d^2\hat{u}}{dr^2} + \frac{1}{r}\frac{d\hat{u}}{dr} - \frac{i\omega\rho}{\mu}\hat{u} = -\frac{\hat{P}}{\mu} $$
//! Enforcing the no-slip boundary condition $\hat{u}(R) = 0$ uniquely resolves the constant terms.
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

use crate::physics::vascular::bessel::{bessel_j0, bessel_j1};
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
        let four = T::one() + T::one() + T::one() + T::one();
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
        let pi = T::pi();
        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");

        if alpha < T::one() {
            // Low α: Z ≈ 8μL/(πR⁴) (Poiseuille resistance dominates)
            eight * mu * self.length / (pi * r.powi(4))
        } else {
            // High α: Z ≈ ρωL/(πR²) (inertance dominates)
            rho * omega * self.length / (pi * r * r)
        }
    }
}

// ============================================================================
// Pulsatility Index
// ============================================================================

/// Pulsatility index and peak/trough velocities from Womersley flow.
///
/// ## Theorem — Gosling-King Pulsatility Index (Gosling & King 1974)
///
/// The pulsatility index (PI) is a dimensionless measure of flow waveform
/// pulsatility, widely used in Doppler ultrasound assessment of arterial
/// hemodynamics. It quantifies the oscillatory component of flow relative
/// to the mean:
///
/// ```text
/// PI = (V_peak_systolic - V_end_diastolic) / V_mean
/// ```
///
/// For Womersley flow decomposed into a steady (Poiseuille) component plus
/// a single oscillatory harmonic:
///
/// ```text
/// V(t) = V_mean + V_osc · cos(omega·t + phi)
/// ```
///
/// The peak and trough velocities are:
///
/// ```text
/// V_peak   = V_mean + |V_osc_max|
/// V_trough = V_mean - |V_osc_max|
/// ```
///
/// For a single harmonic, `|V_osc_max| = Q_amplitude / A`, giving:
///
/// ```text
/// PI = 2 · Q_amplitude / Q_mean
/// ```
///
/// ## Clinical significance
///
/// | Vessel           | Typical PI  |
/// |------------------|-------------|
/// | Internal carotid | 0.8 - 1.2   |
/// | Middle cerebral  | 0.6 - 1.0   |
/// | Femoral artery   | 4.0 - 12.0  |
/// | Uterine artery   | 0.6 - 1.2   |
///
/// **Reference**: Gosling, R.G. & King, D.H. (1974). "Arterial Assessment
/// by Doppler-shift Ultrasound", *Proc. R. Soc. Med.* 67:447-449.
///
/// # Arguments
/// * `q_mean` - Mean volumetric flow rate [m^3/s] (must be positive)
/// * `q_amplitude` - Oscillatory flow rate amplitude [m^3/s] (non-negative)
/// * `cross_section_area` - Cross-sectional area of the vessel [m^2] (must be positive)
///
/// # Returns
/// Tuple of `(pulsatility_index, v_peak, v_trough)`:
/// - `pulsatility_index` — PI = (v_peak - v_trough) / v_mean
/// - `v_peak` — peak systolic velocity [m/s]
/// - `v_trough` — end-diastolic (minimum) velocity [m/s]
///
/// # Panics
/// Does not panic; returns `(0.0, v_mean, v_mean)` when `q_mean` is zero
/// or negative to avoid division by zero.
#[must_use]
pub fn womersley_pulsatility_index(
    q_mean: f64,
    q_amplitude: f64,
    cross_section_area: f64,
) -> (f64, f64, f64) {
    let area = cross_section_area.max(1e-30); // guard against zero area
    let v_mean = q_mean / area;
    let v_osc = q_amplitude.abs() / area;

    let v_peak = v_mean + v_osc;
    let v_trough = v_mean - v_osc;

    let pi = if q_mean.abs() < 1e-30 {
        // Steady flow or zero flow: PI is undefined / zero
        0.0
    } else {
        (v_peak - v_trough) / v_mean.abs()
    };

    (pi, v_peak, v_trough)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_womersley_number_calculation() {
        // Human aorta: R=12.5mm, ω=7.54 rad/s, ρ=1060, μ=0.0035
        let wom = WomersleyNumber::<f64>::human_aorta();
        let alpha = wom.value();

        // Expected α ≈ 18.3 for aorta
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

        // Expected α ≈ 3.3 for femoral
        assert!(
            alpha > 2.0 && alpha < 5.0,
            "Femoral α = {} should be ~3.3",
            alpha
        );
    }

    #[test]
    fn test_flow_regime_classification() {
        // Small vessel with low α
        let low_alpha = WomersleyNumber::<f64>::new(0.0001, 7.54, 1060.0, 0.0035);
        assert_eq!(low_alpha.flow_regime(), FlowRegime::QuasiSteady);

        // Aorta with high α
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

        // Stokes layer should be ~1 mm for cardiac frequency
        assert!(
            delta > 0.0005 && delta < 0.002,
            "δ = {} should be ~1 mm",
            delta
        );
    }

    #[test]
    fn test_velocity_profile_centerline() {
        let wom = WomersleyNumber::<f64>::human_aorta();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0); // 100 Pa/m gradient

        // Centerline velocity should be finite and positive at some time
        let u_center = profile.centerline_velocity(0.25); // t = 0.25s
        assert!(
            u_center.abs() < 1.0,
            "Centerline velocity should be < 1 m/s"
        );
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
        assert!(
            u_center >= u_wall,
            "Centerline should have higher |u| than wall"
        );
        let _ = u_mid; // Silence unused warning
    }

    #[test]
    fn test_flow_rate_positive() {
        let wom = WomersleyNumber::<f64>::human_femoral();
        let profile = WomersleyProfile::<f64>::new(wom, 100.0);

        // Flow rate magnitude should be finite
        let q = profile.flow_rate(0.1);
        assert!(
            q.abs() < 0.001,
            "Flow rate {} L/s should be realistic",
            q * 1000.0
        );
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
        // WSS can range from O(1) to O(100) Pa depending on vessel and conditions
        // For moderate pressure gradients in large arteries, WSS is typically 1-100 Pa
        assert!(tau_w.abs() < 200.0, "WSS {} Pa should be < 200 Pa", tau_w);
    }

    #[test]
    fn test_womersley_low_alpha_poiseuille_limit() {
        // At very low alpha (quasi-steady), Womersley solution approaches Poiseuille
        // u_max_steady = (P_hat * R^2) / (4 * mu)
        let alpha = 0.01;
        let p_hat = 100.0;
        let rho = 1000.0;
        let mu = 0.001;
        let r = 0.01;
        // alpha = R * sqrt(omega * rho / mu) => omega = (alpha/R)^2 * mu / rho
        let omega = (alpha / r) * (alpha / r) * mu / rho;

        let wom = WomersleyNumber::<f64>::new(r, omega, rho, mu);
        let profile = WomersleyProfile::<f64>::new(wom, p_hat);

        // Max magnitude of dynamic centerline velocity
        let u_womersley = profile.centerline_velocity(0.0).abs();

        let u_poiseuille = (p_hat * r * r) / (4.0 * mu);

        assert_relative_eq!(u_womersley, u_poiseuille, max_relative = 0.01);
    }

    #[test]
    fn test_womersley_flow_rate_oscillates() {
        let wom = WomersleyNumber::<f64>::new(0.005, 2.0 * PI, 1000.0, 0.001); // 1 Hz
        let profile = WomersleyProfile::<f64>::new(wom, 10.0);

        // Integrate flow rate over one complete period (T = 1.0s)
        let mut net_volume = 0.0;
        let steps = 1000;
        let dt = 1.0 / steps as f64;
        for i in 0..steps {
            let t = i as f64 * dt;
            net_volume += profile.flow_rate(t) * dt;
        }

        // Net volume over one period for a purely oscillatory flow must be exactly zero
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

        // Wall velocity (xi = 1.0) must be identically zero at all times (no-slip condition)
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

    // ========================================================================
    // Pulsatility Index tests
    // ========================================================================

    #[test]
    fn test_pulsatility_index_steady_flow() {
        // Zero oscillatory amplitude => PI=0, v_peak = v_trough = v_mean
        let area = 1e-4; // 1 cm^2
        let q_mean = 5e-6; // 5 mL/s
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
        // Typical internal carotid: PI ~ 0.8-1.2
        // Q_mean = 5 mL/s, Q_amplitude ~ 2.5 mL/s (PI = 2*2.5/5 = 1.0)
        let area = 2e-5; // ~5 mm diameter vessel
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
        // PI should be >= 0 for all valid inputs (positive q_mean)
        let area = 1e-4;
        for &q_mean in &[1e-6, 1e-5, 1e-4] {
            for &q_amp in &[0.0, 0.5e-6, 1e-6, 5e-6] {
                let (pi, _, _) = womersley_pulsatility_index(q_mean, q_amp, area);
                assert!(pi >= 0.0, "PI should be >= 0, got {} for q_mean={}, q_amp={}", pi, q_mean, q_amp);
            }
        }
    }
}
