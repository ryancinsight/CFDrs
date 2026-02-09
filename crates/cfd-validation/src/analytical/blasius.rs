//! Blasius boundary layer - self-similar solution for laminar flow over flat plate
//!
//! The Blasius solution describes the velocity profile in a laminar boundary layer
//! over a semi-infinite flat plate. It's a self-similar solution to the
//! Prandtl boundary layer equations.
//!
//! # Mathematical Formulation
//! The solution uses a similarity variable η = y * sqrt(U/(νx))
//! and a stream function ψ = sqrt(νxU) * f(η)
//! where f satisfies: f''' + (1/2)ff'' = 0
//! with boundary conditions: f(0) = f'(0) = 0, f'(∞) = 1
//!
//! # References
//! - Blasius, H. (1908). "Grenzschichten in Flüssigkeiten mit kleiner Reibung"
//!   Z. Math. Phys., 56:1-37.
//! - Schlichting, H. (1979). "Boundary Layer Theory", 7th ed., McGraw-Hill.

use super::AnalyticalSolution;
use cfd_core::conversion::SafeFromF64;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Blasius boundary layer solution
#[derive(Debug, Clone)]
pub struct BlasiusBoundaryLayer<T: RealField + Copy> {
    /// Free-stream velocity
    pub u_inf: T,
    /// Kinematic viscosity
    pub nu: T,
    /// Position along plate (distance from leading edge)
    pub x: T,
}

/// Tabulated Blasius solution values (η, f, f', f'')
/// f' is the normalized velocity u/U
const BLASIUS_TABLE: [(f64, f64, f64, f64); 41] = [
    (0.0, 0.0, 0.0, 0.33206),
    (0.2, 0.00664, 0.06641, 0.33199),
    (0.4, 0.02656, 0.13277, 0.33147),
    (0.6, 0.05974, 0.19894, 0.33008),
    (0.8, 0.10611, 0.26471, 0.32739),
    (1.0, 0.16557, 0.32979, 0.32301),
    (1.2, 0.23795, 0.39378, 0.31659),
    (1.4, 0.32298, 0.45627, 0.30787),
    (1.6, 0.42032, 0.51676, 0.29667),
    (1.8, 0.52952, 0.57477, 0.28293),
    (2.0, 0.65003, 0.62977, 0.26675),
    (2.2, 0.78120, 0.68132, 0.24835),
    (2.4, 0.92230, 0.72899, 0.22809),
    (2.6, 1.07252, 0.77246, 0.20646),
    (2.8, 1.23099, 0.81152, 0.18401),
    (3.0, 1.39682, 0.84605, 0.16136),
    (3.2, 1.56911, 0.87609, 0.13913),
    (3.4, 1.74696, 0.90177, 0.11788),
    (3.6, 1.92954, 0.92333, 0.09809),
    (3.8, 2.11605, 0.94112, 0.08013),
    (4.0, 2.30576, 0.95552, 0.06424),
    (4.2, 2.49806, 0.96696, 0.05052),
    (4.4, 2.69238, 0.97587, 0.03897),
    (4.6, 2.88826, 0.98269, 0.02948),
    (4.8, 3.08534, 0.98779, 0.02187),
    (5.0, 3.28329, 0.99155, 0.01591),
    (5.2, 3.48189, 0.99425, 0.01134),
    (5.4, 3.68094, 0.99616, 0.00793),
    (5.6, 3.88031, 0.99748, 0.00543),
    (5.8, 4.07990, 0.99838, 0.00365),
    (6.0, 4.27964, 0.99898, 0.00240),
    (6.2, 4.47948, 0.99937, 0.00155),
    (6.4, 4.67938, 0.99961, 0.00098),
    (6.6, 4.87931, 0.99977, 0.00061),
    (6.8, 5.07928, 0.99987, 0.00037),
    (7.0, 5.27926, 0.99992, 0.00022),
    (7.2, 5.47925, 0.99996, 0.00013),
    (7.4, 5.67924, 0.99998, 0.00007),
    (7.6, 5.87924, 0.99999, 0.00004),
    (7.8, 6.07923, 1.00000, 0.00002),
    (8.0, 6.27923, 1.00000, 0.00001),
];

impl<T: RealField + Copy + FromPrimitive + num_traits::ToPrimitive> BlasiusBoundaryLayer<T> {
    /// Create new Blasius boundary layer
    pub fn new(u_inf: T, nu: T, x: T) -> Self {
        Self { u_inf, nu, x }
    }

    /// Create air flow over flat plate at standard conditions
    pub fn air_flow(u_inf: f64) -> Self {
        Self {
            u_inf: T::from_f64_or_one(u_inf),
            nu: T::from_f64_or_one(1.5e-5), // Air kinematic viscosity
            x: T::from_f64_or_one(1.0),
        }
    }

    /// Create water flow over flat plate
    pub fn water_flow(u_inf: f64) -> Self {
        Self {
            u_inf: T::from_f64_or_one(u_inf),
            nu: T::from_f64_or_one(1.0e-6), // Water kinematic viscosity
            x: T::from_f64_or_one(1.0),
        }
    }

    /// Calculate local Reynolds number: Re_x = U*x/ν
    pub fn local_reynolds(&self) -> T {
        self.u_inf * self.x / self.nu
    }

    /// Calculate boundary layer thickness (99% of free-stream velocity)
    /// δ ≈ 5.0 * sqrt(νx/U) = 5.0x / sqrt(Re_x)
    pub fn boundary_layer_thickness(&self) -> T {
        let five = T::from_f64_or_one(5.0);
        five * (self.nu * self.x / self.u_inf).sqrt()
    }

    /// Calculate displacement thickness
    /// δ* ≈ 1.72 * sqrt(νx/U) = 1.72x / sqrt(Re_x)
    pub fn displacement_thickness(&self) -> T {
        let factor = T::from_f64_or_one(1.7208);
        factor * (self.nu * self.x / self.u_inf).sqrt()
    }

    /// Calculate momentum thickness
    /// θ ≈ 0.664 * sqrt(νx/U) = 0.664x / sqrt(Re_x)
    pub fn momentum_thickness(&self) -> T {
        let factor = T::from_f64_or_one(0.664);
        factor * (self.nu * self.x / self.u_inf).sqrt()
    }

    /// Calculate shape factor H = δ*/θ
    /// For Blasius: H ≈ 2.59
    pub fn shape_factor(&self) -> T {
        let h_blasius = T::from_f64_or_one(2.591);
        h_blasius
    }

    /// Calculate wall shear stress
    /// τ_w = 0.332 * μ * U * sqrt(U/(νx))
    pub fn wall_shear_stress(&self) -> T {
        let factor = T::from_f64_or_one(0.332);
        let mu = self.nu; // Assuming unit density for simplicity
        factor * mu * self.u_inf * (self.u_inf / (self.nu * self.x)).sqrt()
    }

    /// Calculate skin friction coefficient
    /// Cf = 0.664 / sqrt(Re_x)
    pub fn skin_friction_coefficient(&self) -> T {
        let factor = T::from_f64_or_one(0.664);
        factor / self.local_reynolds().sqrt()
    }

    /// Calculate similarity variable η = y * sqrt(U/(νx))
    pub fn similarity_variable(&self, y: T) -> T {
        y * (self.u_inf / (self.nu * self.x)).sqrt()
    }

    /// Get normalized velocity u/U at similarity variable η
    /// Uses linear interpolation from tabulated Blasius solution
    fn velocity_ratio_at_eta(&self, eta: T) -> T {
        let eta_f64 = eta.to_f64().unwrap_or(0.0);
        let zero = T::zero();
        let one = T::one();
        let eight = T::from_f64_or_one(8.0);

        if eta <= zero {
            return zero;
        }
        if eta >= eight {
            return one;
        }

        // Find bracketing points
        let mut i = 0;
        while i < BLASIUS_TABLE.len() - 1 && BLASIUS_TABLE[i + 1].0 < eta_f64 {
            i += 1;
        }

        let (eta1, _, f1, _) = BLASIUS_TABLE[i];
        let (eta2, _, f2, _) = BLASIUS_TABLE[i + 1];

        // Linear interpolation
        let frac = T::from_f64_or_one((eta_f64 - eta1) / (eta2 - eta1));
        T::from_f64_or_one(f1) + frac * (T::from_f64_or_one(f2) - T::from_f64_or_one(f1))
    }

    /// Get velocity at physical coordinates (x, y)
    /// Note: x should be >= self.x (downstream of leading edge)
    pub fn velocity_at(&self, x: T, y: T) -> T {
        // Adjust for different x positions
        let local_eta = y * (self.u_inf / (self.nu * x)).sqrt();
        self.velocity_ratio_at_eta(local_eta) * self.u_inf
    }

    /// Get wall-normal velocity component v/U at similarity variable
    /// v/U = (1/2) * sqrt(ν/(Ux)) * (ηf' - f)
    fn normal_velocity_ratio_at_eta(&self, eta: T) -> T {
        let eta_f64 = eta.to_f64().unwrap_or(0.0);
        let zero = T::zero();
        let eight = T::from_f64_or_one(8.0);

        if eta <= zero {
            return zero;
        }
        if eta >= eight {
            return zero; // Goes to zero at edge
        }

        // Find bracketing points
        let mut i = 0;
        while i < BLASIUS_TABLE.len() - 1 && BLASIUS_TABLE[i + 1].0 < eta_f64 {
            i += 1;
        }

        let (eta1, f1, fp1, _) = BLASIUS_TABLE[i];
        let (eta2, f2, fp2, _) = BLASIUS_TABLE[i + 1];

        // Linear interpolation
        let frac = (eta_f64 - eta1) / (eta2 - eta1);
        let f = f1 + frac * (f2 - f1);
        let fp = fp1 + frac * (fp2 - fp1);

        // Calculate (ηf' - f) / 2
        T::from_f64_or_one((eta_f64 * fp - f) * 0.5)
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::ToPrimitive> AnalyticalSolution<T> for BlasiusBoundaryLayer<T> {
    fn evaluate(&self, x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        if x <= T::zero() {
            return Vector3::new(self.u_inf, T::zero(), T::zero());
        }

        let eta = self.similarity_variable(y);
        let u = self.velocity_ratio_at_eta(eta) * self.u_inf;

        // v-velocity (wall-normal)
        let sqrt_term = (self.nu / (self.u_inf * x)).sqrt();
        let v_ratio = self.normal_velocity_ratio_at_eta(eta);
        let v = sqrt_term * v_ratio * self.u_inf;

        Vector3::new(u, v, T::zero())
    }

    fn pressure(&self, _x: T, _y: T, _z: T, _t: T) -> T {
        // Constant pressure across boundary layer (first-order boundary layer theory)
        T::zero()
    }

    fn name(&self) -> &str {
        "Blasius Boundary Layer"
    }

    fn domain_bounds(&self) -> [T; 6] {
        let delta = self.boundary_layer_thickness();
        let length = T::from_f64_or_one(10.0) * self.x;
        [
            T::zero(),
            length,           // x: [0, 10x]
            T::zero(),
            delta * T::from_f64_or_one(5.0), // y: [0, 5δ]
            T::zero(),
            T::zero(),        // z: 0 (2D flow)
        ]
    }

    fn length_scale(&self) -> T {
        self.boundary_layer_thickness()
    }

    fn velocity_scale(&self) -> T {
        self.u_inf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blasius_boundary_conditions() {
        let blasius = BlasiusBoundaryLayer::<f64>::air_flow(1.0);

        // At wall (η = 0): u = 0, v = 0
        let u_wall = blasius.velocity_ratio_at_eta(0.0_f64);
        assert!(u_wall.abs() < 1e-6, "Velocity at wall should be zero");

        // At infinity (η → ∞): u → U
        let u_inf = blasius.velocity_ratio_at_eta(8.0_f64);
        assert!((u_inf - 1.0).abs() < 0.001, "Velocity at infinity should be 1.0");
    }

    #[test]
    fn test_boundary_layer_growth() {
        let blasius1 = BlasiusBoundaryLayer::new(1.0, 1.5e-5, 1.0);
        let blasius2 = BlasiusBoundaryLayer::new(1.0, 1.5e-5, 4.0);

        // Boundary layer grows with sqrt(x)
        let delta1 = blasius1.boundary_layer_thickness();
        let delta2 = blasius2.boundary_layer_thickness();

        let ratio = delta2 / delta1;
        assert!((ratio - 2.0_f64).abs() < 0.01, "Boundary layer should grow as sqrt(x)");
    }

    #[test]
    fn test_shape_factor() {
        let blasius = BlasiusBoundaryLayer::<f64>::air_flow(10.0);

        // Calculate shape factor from definitions
        let delta_star = blasius.displacement_thickness();
        let theta = blasius.momentum_thickness();
        let h_calc = delta_star / theta;

        // Should match theoretical value of 2.59
        assert!((h_calc - 2.591).abs() < 0.01,
            "Shape factor should be ~2.59 for Blasius");
    }

    #[test]
    fn test_skin_friction_scaling() {
        // Cf should decrease with Reynolds number
        let blasius_low_re = BlasiusBoundaryLayer::new(1.0, 1.5e-5, 0.1);
        let blasius_high_re = BlasiusBoundaryLayer::new(1.0, 1.5e-5, 10.0);

        let cf_low = blasius_low_re.skin_friction_coefficient();
        let cf_high = blasius_high_re.skin_friction_coefficient();

        assert!(cf_low > cf_high, "Skin friction should decrease with Re");
    }
}