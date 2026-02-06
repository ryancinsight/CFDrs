//! Non-Newtonian Poiseuille flow in 2D channels
//!
//! Analytical and semi-analytical solutions for non-Newtonian fluid flow
//! between parallel plates, including power-law and Casson models.
//!
//! # Mathematical Foundation
//!
//! ## Power-Law Fluid
//! For a power-law fluid: τ = K·γ̇ⁿ in a 2D channel (width 2H):
//! ```text
//! u(y) = u_c · [1 - |y/H|^((n+1)/n)]^(n/(n+1))
//! ```
//! where:
//! - u_c = centerline velocity
//! - n = power-law index (n<1: shear-thinning, n>1: shear-thickening, n=1: Newtonian)
//! - K = consistency index [Pa·sⁿ]
//!
//! ## Casson Blood Model
//! For Casson fluid: √τ = √τ_y + √(μ_∞·γ̇)
//! The velocity profile must be solved numerically, but satisfies:
//! ```text
//! γ̇(y) = -du/dy > 0 for y > 0
//! τ(y) = (dp/dx)·y (linear stress distribution)
//! ```
//!
//! # References
//! - Bird, R.B., Stewart, W.E., Lightfoot, E.N. (2002) "Transport Phenomena"
//! - Chhabra, R.P., Richardson, J.F. (2008) "Non-Newtonian Flow and Applied Rheology"
//! - Merrill, E.W. et al. (1969) "Rheology of blood"

use super::AnalyticalSolution;
use cfd_core::physics::fluid::blood::CassonBlood;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

// ============================================================================
// Rheological Model Trait
// ============================================================================

/// Trait for rheological models in Poiseuille flow
pub trait RheologicalModel<T: RealField + Copy> {
    /// Compute shear stress from shear rate
    fn shear_stress(&self, shear_rate: T) -> T;
    
    /// Compute viscosity from shear rate
    fn viscosity(&self, shear_rate: T) -> T;
    
    /// Name of the rheological model
    fn model_name(&self) -> &str;
}

// ============================================================================
// Power-Law Model
// ============================================================================

/// Power-law rheological model: τ = K·γ̇ⁿ
#[derive(Debug, Clone, Copy)]
pub struct PowerLawModel<T: RealField + Copy> {
    /// Consistency index K [Pa·sⁿ]
    pub consistency: T,
    /// Power-law index n (dimensionless)
    pub index: T,
}

impl<T: RealField + FromPrimitive + Copy> PowerLawModel<T> {
    /// Create power-law model
    ///
    /// # Arguments
    /// * `consistency` - K [Pa·sⁿ]
    /// * `index` - n (n<1: shear-thinning, n=1: Newtonian, n>1: shear-thickening)
    pub fn new(consistency: T, index: T) -> Self {
        Self { consistency, index }
    }
    
    /// Create Newtonian model (n=1)
    pub fn newtonian(viscosity: T) -> Self {
        Self {
            consistency: viscosity,
            index: T::one(),
        }
    }
    
    /// Create typical blood-like shear-thinning model (n=0.6)
    pub fn blood_like(consistency: T) -> Self {
        Self {
            consistency,
            index: T::from_f64(0.6).unwrap(),
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> RheologicalModel<T> for PowerLawModel<T> {
    fn shear_stress(&self, shear_rate: T) -> T {
        self.consistency * shear_rate.abs().powf(self.index) * shear_rate.signum()
    }
    
    fn viscosity(&self, shear_rate: T) -> T {
        if shear_rate.abs() < T::from_f64(1e-12).unwrap() {
            // Avoid division by zero at centerline
            return T::from_f64(1e12).unwrap();
        }
        self.consistency * shear_rate.abs().powf(self.index - T::one())
    }
    
    fn model_name(&self) -> &str {
        "Power-Law"
    }
}

// ============================================================================
// Non-Newtonian Poiseuille Flow - Power-Law
// ============================================================================

/// Non-Newtonian Poiseuille flow with power-law rheology
///
/// Analytical solution for flow between parallel plates with half-width H.
#[derive(Debug, Clone)]
pub struct PowerLawPoiseuille<T: RealField + Copy> {
    /// Power-law model
    pub model: PowerLawModel<T>,
    /// Channel half-width H [m]
    pub half_width: T,
    /// Pressure gradient magnitude |dp/dx| [Pa/m] (positive value)
    pub pressure_gradient: T,
    /// Channel length for domain [m]
    pub length: T,
}

impl<T: RealField + FromPrimitive + Copy> PowerLawPoiseuille<T> {
    /// Create power-law Poiseuille flow
    ///
    /// # Arguments
    /// * `consistency` - K [Pa·sⁿ]
    /// * `index` - n (power-law index)
    /// * `half_width` - H [m]
    /// * `pressure_gradient` - |dp/dx| [Pa/m]
    /// * `length` - L [m] (domain length)
    pub fn new(
        consistency: T,
        index: T,
        half_width: T,
        pressure_gradient: T,
        length: T,
    ) -> Self {
        Self {
            model: PowerLawModel::new(consistency, index),
            half_width,
            pressure_gradient,
            length,
        }
    }
    
    /// Centerline velocity u_c
    ///
    /// For power-law fluid:
    /// ```text
    /// u_c = [(n/(n+1)) · (H^(n+1)/K) · (dp/dx)^(1/n)]
    /// ```
    pub fn centerline_velocity(&self) -> T {
        let n = self.model.index;
        let one = T::one();
        let h = self.half_width;
        let dp_dx = self.pressure_gradient;
        let k = self.model.consistency;
        
        // u_c = [n/(n+1)] · (H/K)^(1/n) · (dp/dx·H)^(1/n)
        let factor = n / (n + one);
        let term = (dp_dx * h) / k;
        
        factor * h * term.powf(one / n)
    }
    
    /// Velocity profile u(y) for power-law fluid
    ///
    /// ```text
    /// u(y) = u_c · [1 - |y/H|^((n+1)/n)]
    /// ```
    pub fn velocity_at(&self, y: T) -> T {
        let n = self.model.index;
        let h = self.half_width;
        let u_c = self.centerline_velocity();
        
        let y_normalized = (y / h).abs();
        
        if y_normalized >= T::one() {
            // Outside channel
            return T::zero();
        }
        
        // Exponent: (n+1)/n
        let exponent = (n + T::one()) / n;
        u_c * (T::one() - y_normalized.powf(exponent))
    }
    
    /// Flow rate per unit depth Q' [m²/s]
    ///
    /// ```text
    /// Q' = 2·∫[0 to H] u(y) dy
    /// ```
    pub fn flow_rate_per_depth(&self) -> T {
        let n = self.model.index;
        let h = self.half_width;
        let u_c = self.centerline_velocity();
        
        // Q' = 2·u_c·H·[n/(2n+1)]
        let two = T::from_f64(2.0).unwrap();
        let factor = n / (two * n + T::one());
        
        two * u_c * h * factor
    }
    
    /// Wall shear stress τ_w [Pa]
    ///
    /// ```text
    /// τ_w = H · dp/dx
    /// ```
    pub fn wall_shear_stress(&self) -> T {
        self.half_width * self.pressure_gradient
    }
    
    /// Wall shear rate γ̇_w [1/s]
    ///
    /// From τ_w = K·γ̇_w^n:
    /// ```text
    /// γ̇_w = (τ_w / K)^(1/n)
    /// ```
    pub fn wall_shear_rate(&self) -> T {
        let tau_w = self.wall_shear_stress();
        let k = self.model.consistency;
        let n = self.model.index;
        
        (tau_w / k).powf(T::one() / n)
    }
    
    /// Reynolds number for power-law fluid
    ///
    /// Defined using generalized Reynolds number:
    /// ```text
    /// Re_PL = ρ·u_c^(2-n)·H^n / K
    /// ```
    pub fn reynolds_number(&self, density: T) -> T {
        let u_c = self.centerline_velocity();
        let n = self.model.index;
        let h = self.half_width;
        let k = self.model.consistency;
        
        density * u_c.powf(T::from_f64(2.0).unwrap() - n) * h.powf(n) / k
    }
}

impl<T: RealField + FromPrimitive + Copy> AnalyticalSolution<T> for PowerLawPoiseuille<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        let u = self.velocity_at(y);
        Vector3::new(u, T::zero(), T::zero())
    }
    
    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop: p(x) = p0 - (dp/dx)·x
        -self.pressure_gradient * x
    }
    
    fn name(&self) -> &str {
        "Power-Law Poiseuille (2D Channel)"
    }
    
    fn domain_bounds(&self) -> [T; 6] {
        [
            T::zero(),
            self.length,                 // x: [0, L]
            -self.half_width,
            self.half_width,             // y: [-H, H]
            T::zero(),
            T::zero(),                   // z: 0 (2D)
        ]
    }
    
    fn length_scale(&self) -> T {
        self.half_width
    }
    
    fn velocity_scale(&self) -> T {
        self.centerline_velocity()
    }
}

// ============================================================================
// Casson Blood Model Wrapper
// ============================================================================

/// Wrapper for Casson blood rheology in Poiseuille flow
impl<T: RealField + FromPrimitive + Copy> RheologicalModel<T> for CassonBlood<T> {
    fn shear_stress(&self, shear_rate: T) -> T {
        let gamma_dot = shear_rate.abs();
        let sqrt_tau_y = self.yield_stress.sqrt();
        let sqrt_mu_inf = self.infinite_shear_viscosity.sqrt();
        let sqrt_gamma = gamma_dot.sqrt();
        
        let sqrt_tau = sqrt_tau_y + sqrt_mu_inf * sqrt_gamma;
        sqrt_tau * sqrt_tau * shear_rate.signum()
    }
    
    fn viscosity(&self, shear_rate: T) -> T {
        self.apparent_viscosity(shear_rate.abs())
    }
    
    fn model_name(&self) -> &str {
        "Casson"
    }
}

// ============================================================================
// Non-Newtonian Poiseuille - Casson Blood (Numerical Solution)
// ============================================================================

/// Casson blood Poiseuille flow (numerical solution required)
///
/// The velocity profile for Casson fluid must be solved numerically.
/// We provide numerical integration and validation methods.
#[derive(Debug, Clone)]
pub struct CassonPoiseuille<T: RealField + Copy> {
    /// Casson blood model
    pub model: CassonBlood<T>,
    /// Channel half-width H [m]
    pub half_width: T,
    /// Pressure gradient magnitude |dp/dx| [Pa/m]
    pub pressure_gradient: T,
    /// Channel length [m]
    pub length: T,
    /// Plug flow radius (where τ < τ_y) [m]
    pub plug_radius: T,
}

impl<T: RealField + FromPrimitive + Copy> CassonPoiseuille<T> {
    /// Create Casson Poiseuille flow
    ///
    /// # Arguments
    /// * `model` - Casson blood rheology
    /// * `half_width` - H [m]
    /// * `pressure_gradient` - |dp/dx| [Pa/m]
    /// * `length` - L [m]
    pub fn new(
        model: CassonBlood<T>,
        half_width: T,
        pressure_gradient: T,
        length: T,
    ) -> Self {
        // Calculate plug radius: y_p where τ(y_p) = τ_y
        // τ(y) = (dp/dx)·y, so y_p = τ_y / (dp/dx)
        let tau_y = model.yield_stress;
        let y_p = tau_y / pressure_gradient;
        
        // Clamp to channel width
        let plug_radius = if y_p > half_width {
            half_width
        } else {
            y_p
        };
        
        Self {
            model,
            half_width,
            pressure_gradient,
            length,
            plug_radius,
        }
    }
    
    /// Velocity in plug region (|y| ≤ y_p)
    ///
    /// In the plug region where τ < τ_y, the fluid moves as a rigid body
    /// with velocity equal to the boundary velocity at y = y_p.
    pub fn plug_velocity(&self) -> T {
        // Integrate from plug boundary to wall
        // This requires numerical integration of du/dy = γ̇(y)
        self.velocity_at_numerical(self.plug_radius)
    }
    
    /// Velocity profile u(y) - numerical integration
    ///
    /// Integrate du/dy from y to H:
    /// ```text
    /// u(y) = ∫[y to H] γ̇(η) dη
    /// ```
    /// where γ̇(η) is obtained from Casson equation inverted.
    pub fn velocity_at_numerical(&self, y: T) -> T {
        let y_abs = y.abs();
        
        if y_abs >= self.half_width {
            return T::zero();
        }
        
        if y_abs <= self.plug_radius {
            // Plug region: constant velocity
            return self.plug_velocity();
        }
        
        // Numerical integration using Simpson's rule
        let n_points = 100;
        let dy = (self.half_width - y_abs) / T::from_usize(n_points).unwrap();
        let mut integral = T::zero();
        
        for i in 0..=n_points {
            let eta = y_abs + T::from_usize(i).unwrap() * dy;
            let tau = self.pressure_gradient * eta;
            
            // Casson equation: √τ = √τ_y + √(μ_∞·γ̇)
            // Solve for γ̇: γ̇ = [(√τ - √τ_y)² / μ_∞]
            let sqrt_tau_y = self.model.yield_stress.sqrt();
            let sqrt_tau = tau.sqrt();
            
            let gamma_dot = if tau > self.model.yield_stress {
                let sqrt_mu_inf = self.model.infinite_shear_viscosity.sqrt();
                let diff = sqrt_tau - sqrt_tau_y;
                (diff / sqrt_mu_inf) * (diff / sqrt_mu_inf)
            } else {
                T::zero()
            };
            
            // Simpson's rule weights
            let weight = if i == 0 || i == n_points {
                T::one()
            } else if i % 2 == 0 {
                T::from_f64(2.0).unwrap()
            } else {
                T::from_f64(4.0).unwrap()
            };
            
            integral = integral + weight * gamma_dot;
        }
        
        integral * dy / T::from_f64(3.0).unwrap()
    }
    
    /// Wall shear stress [Pa]
    pub fn wall_shear_stress(&self) -> T {
        self.half_width * self.pressure_gradient
    }
    
    /// Flow rate per unit depth (numerical integration)
    pub fn flow_rate_per_depth(&self) -> T {
        // Q' = 2·∫[0 to H] u(y) dy
        let n_points = 100;
        let dy = self.half_width / T::from_usize(n_points).unwrap();
        let mut integral = T::zero();
        
        for i in 0..=n_points {
            let y = T::from_usize(i).unwrap() * dy;
            let u = self.velocity_at_numerical(y);
            
            let weight = if i == 0 || i == n_points {
                T::one()
            } else if i % 2 == 0 {
                T::from_f64(2.0).unwrap()
            } else {
                T::from_f64(4.0).unwrap()
            };
            
            integral = integral + weight * u;
        }
        
        T::from_f64(2.0).unwrap() * integral * dy / T::from_f64(3.0).unwrap()
    }
}

impl<T: RealField + FromPrimitive + Copy> AnalyticalSolution<T> for CassonPoiseuille<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        let u = self.velocity_at_numerical(y);
        Vector3::new(u, T::zero(), T::zero())
    }
    
    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        -self.pressure_gradient * x
    }
    
    fn name(&self) -> &str {
        "Casson Blood Poiseuille (2D Channel)"
    }
    
    fn domain_bounds(&self) -> [T; 6] {
        [
            T::zero(),
            self.length,
            -self.half_width,
            self.half_width,
            T::zero(),
            T::zero(),
        ]
    }
    
    fn length_scale(&self) -> T {
        self.half_width
    }
    
    fn velocity_scale(&self) -> T {
        self.plug_velocity()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_power_law_newtonian_limit() {
        // For n=1, power-law should match Newtonian Poiseuille
        let mu = 0.001_f64; // Pa·s
        let h = 0.001_f64;  // 1 mm
        let dp_dx = 100.0_f64; // 100 Pa/m
        let l = 0.1_f64;    // 10 cm
        
        let flow = PowerLawPoiseuille::new(mu, 1.0, h, dp_dx, l);
        
        // Newtonian: u_c = H²·dp/dx / (2μ)
        let u_c_expected = h * h * dp_dx / (2.0 * mu);
        let u_c_actual = flow.centerline_velocity();
        
        assert_relative_eq!(u_c_actual, u_c_expected, epsilon = 1e-10);
        
        // Check parabolic profile: u(y) = u_c·(1 - y²/H²)
        let y = 0.5 * h;
        let u_expected = u_c_expected * (1.0 - (y / h).powi(2));
        let u_actual = flow.velocity_at(y);
        
        assert_relative_eq!(u_actual, u_expected, epsilon = 1e-10);
    }
    
    #[test]
    fn test_power_law_shear_thinning() {
        // Shear-thinning (n=0.6) should have higher centerline velocity than Newtonian
        let k = 0.01_f64;   // Pa·s^n
        let n = 0.6_f64;
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;
        
        let flow = PowerLawPoiseuille::<f64>::new(k, n, h, dp_dx, l);
        let u_c = flow.centerline_velocity();
        
        // Should be finite and positive
        assert!(u_c > 0.0 && u_c.is_finite());
        
        // Velocity should decrease from center to wall
        let u_center = flow.velocity_at(0.0);
        let u_mid = flow.velocity_at(0.5 * h);
        let u_wall = flow.velocity_at(h);
        
        assert!(u_center > u_mid);
        assert!(u_mid > u_wall);
        assert_relative_eq!(u_wall, 0.0, epsilon = 1e-10);
    }
    
    #[test]
    fn test_power_law_wall_shear() {
        let k = 0.01_f64;
        let n = 0.6_f64;
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;
        
        let flow = PowerLawPoiseuille::new(k, n, h, dp_dx, l);
        let tau_w = flow.wall_shear_stress();
        
        // τ_w = H·dp/dx
        assert_relative_eq!(tau_w, h * dp_dx, epsilon = 1e-10);
    }
    
    #[test]
    fn test_casson_plug_flow() {
        // High yield stress should create large plug region
        let tau_y = 0.01_f64; // 10 mPa (yield stress)
        let mu_inf = 0.0035_f64; // 3.5 cP
        
        let model = CassonBlood::new(
            1060.0,
            tau_y,
            mu_inf,
            0.45,
        );
        
        let h = 0.001_f64;      // 1 mm
        let dp_dx = 100.0_f64;  // 100 Pa/m
        let l = 0.1_f64;
        
        let flow = CassonPoiseuille::new(model, h, dp_dx, l);
        
        // Plug radius: y_p = τ_y / (dp/dx)
        let y_p_expected = tau_y / dp_dx;
        assert_relative_eq!(flow.plug_radius, y_p_expected, epsilon = 1e-10);
        
        // In plug region, velocity should be constant
        let u_plug = flow.plug_velocity();
        let u_at_half_plug = flow.velocity_at_numerical(0.5 * flow.plug_radius);
        
        // Should be approximately equal (within numerical error)
        assert_relative_eq!(u_at_half_plug, u_plug, epsilon = 0.01);
    }
    
    #[test]
    fn test_casson_wall_velocity_zero() {
        let model = CassonBlood::new(
            1060.0,
            0.0056,
            0.00345,
            0.45,
        );
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;
        
        let flow = CassonPoiseuille::<f64>::new(model, h, dp_dx, l);
        
        // Wall velocity should be zero
        let u_wall = flow.velocity_at_numerical(h);
        assert!(u_wall.abs() < 1e-6, "Wall velocity {} should be ~0", u_wall);
    }
    
    #[test]
    fn test_casson_flow_rate() {
        let model = CassonBlood::new(
            1060.0,
            0.0056,
            0.00345,
            0.45,
        );
        let h = 0.001_f64;
        let dp_dx = 100.0_f64;
        let l = 0.1_f64;
        
        let flow = CassonPoiseuille::<f64>::new(model, h, dp_dx, l);
        let q = flow.flow_rate_per_depth();
        
        // Flow rate should be positive and finite
        assert!(q > 0.0 && q.is_finite(), "Flow rate {} invalid", q);
    }
}
