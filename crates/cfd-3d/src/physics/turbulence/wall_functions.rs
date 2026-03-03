//! Near-wall treatment using the law of the wall for RANS models.
//!
//! # Theorem — Law of the Wall (Millikan 1938, Prandtl 1925)
//!
//! For turbulent flow over a smooth wall, the mean velocity profile in
//! wall units u⁺ = u/u_τ (where u_τ = √(τ_w/ρ) is the friction velocity)
//! follows:
//!
//! ```text
//! u⁺ = y⁺                          (viscous sublayer,  y⁺ < 11.53)
//! u⁺ = (1/κ) ln(y⁺) + B           (log-law region,    y⁺ > 30)
//! ```
//!
//! where:
//! - y⁺ = y · u_τ / ν  (dimensionless wall distance)
//! - κ = 0.41          (von Kármán constant, empirical)
//! - B = 5.2           (smooth-wall additive constant, Millikan 1938)
//!
//! **Proof sketch.** In the viscous sublayer, inertia is negligible and
//! molecular viscosity dominates: u/u_τ = y u_τ/ν = y⁺.  In the log-law
//! region, dimensional analysis (Buckingham π) requires the velocity gradient
//! ∂u/∂y ∝ u_τ/(κy), integrating to u⁺ = (1/κ)ln(y⁺) + B (Prandtl 1925).
//!
//! ## Algorithm — Wall Shear Stress from Velocity
//!
//! ```text
//! Given: velocity u at node P (distance y_P from wall), ν, ρ
//! 1. Estimate y⁺ by Newton iteration:
//!    a. Initial guess: u_τ⁰ = sqrt(ν · u / y_P)  (viscous sublayer)
//!    b. For iter in 0..MAX_ITER:
//!       y_plus = y_P · u_τ / ν
//!       u_plus_law = law_of_wall(y_plus)
//!       residual = u_τ · u_plus_law − u
//!       Update u_τ via Newton step
//! 2. τ_w = ρ · u_τ²
//! 3. y_plus = y_P · u_τ / ν
//! ```
//!
//! ## References
//!
//! - Millikan, C.B. (1938). "A critical discussion of turbulent flows in
//!   channels and circular tubes." *Proc. Fifth Int. Congr. Appl. Mech.*
//!   386–392.
//! - Prandtl, L. (1925). "Über die ausgebildete Turbulenz." *Z. Angew.
//!   Math. Mech.* 5:136–139.
//! - Pope, S.B. (2000). *Turbulent Flows*. Cambridge University Press, §7.3.

use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{WALL_B, WALL_KAPPA, WALL_Y_PLUS_TRANSITION};

/// Wall functions for near-wall turbulence treatment.
///
/// Provides y⁺ calculation, law-of-the-wall velocity profile, and wall shear
/// stress estimation for RANS solvers using wall functions instead of
/// resolving the viscous sublayer.
pub struct WallFunctions;

impl WallFunctions {
    /// Compute the dimensionless velocity u⁺ for a given y⁺ (law of the wall).
    ///
    /// Returns the composite profile:
    /// - Viscous sublayer (y⁺ < 11.53): u⁺ = y⁺
    /// - Log-law region (y⁺ ≥ 11.53): u⁺ = (1/κ) ln(y⁺) + B
    ///
    /// # Arguments
    /// * `y_plus` — dimensionless wall distance y⁺ = y u_τ / ν (≥ 0)
    ///
    /// # Returns
    /// Dimensionless velocity u⁺ (≥ 0)
    pub fn u_plus<T: RealField + Copy + FromPrimitive + num_traits::Float>(y_plus: T) -> T {
        let y_plus_tr = <T as FromPrimitive>::from_f64(WALL_Y_PLUS_TRANSITION).unwrap_or_else(T::one);
        if y_plus < y_plus_tr {
            y_plus // viscous sublayer
        } else {
            let kappa = <T as FromPrimitive>::from_f64(WALL_KAPPA).unwrap_or_else(T::one);
            let b     = <T as FromPrimitive>::from_f64(WALL_B    ).unwrap_or_else(T::zero);
            num_traits::Float::ln(y_plus) / kappa + b
        }
    }

    /// Compute dimensionless wall distance y⁺ = y · u_τ / ν.
    ///
    /// # Arguments
    /// * `y` — physical wall distance [m]
    /// * `u_tau` — friction velocity u_τ = √(τ_w/ρ) [m/s]
    /// * `nu` — kinematic viscosity [m²/s]
    pub fn y_plus<T: RealField + Copy>(y: T, u_tau: T, nu: T) -> T {
        y * u_tau / nu
    }

    /// Estimate friction velocity u_τ from near-wall velocity using the law of
    /// the wall (Newton-Raphson iteration).
    ///
    /// # Arguments
    /// * `u_wall` — magnitude of velocity at wall-adjacent node P [m/s]
    /// * `y_p` — distance from wall to node P [m]
    /// * `nu` — kinematic viscosity [m²/s]
    ///
    /// # Returns
    /// Friction velocity u_τ [m/s]; returns √(ν·u/y) on failure.
    pub fn friction_velocity<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        u_wall: T,
        y_p: T,
        nu: T,
    ) -> T {
        let eps = <T as FromPrimitive>::from_f64(1e-15).unwrap_or_else(T::zero);
        let half = <T as FromPrimitive>::from_f64(0.5).unwrap_or_else(T::one);

        // Initial guess: viscous sublayer (u_τ = √(ν u / y))
        let mut u_tau = num_traits::Float::sqrt(nu * u_wall / (y_p + eps));

        // Newton-Raphson: residual f = u_τ · u⁺(y⁺) − u_wall = 0
        for _ in 0..20 {
            let y_p_val = Self::y_plus(y_p, u_tau, nu + eps);
            let u_p_val = Self::u_plus(y_p_val);
            let f = u_tau * u_p_val - u_wall;

            // Approximate Jacobian: df/du_τ ≈ u_p (linear approximation)
            let df = u_p_val + eps;
            let du_tau = f / df;
            u_tau = u_tau - half * du_tau;
            u_tau = num_traits::Float::max(u_tau, eps);

            if num_traits::Float::abs(du_tau) < eps * u_tau { break; }
        }
        u_tau
    }

    /// Compute wall shear stress τ_w = ρ u_τ².
    ///
    /// # Arguments
    /// * `u_wall` — near-wall velocity magnitude [m/s]
    /// * `y_p` — wall-adjacent node distance [m]
    /// * `nu` — kinematic viscosity [m²/s]
    /// * `rho` — density [kg/m³]
    pub fn wall_shear_stress<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        u_wall: T,
        y_p: T,
        nu: T,
        rho: T,
    ) -> T {
        let u_tau = Self::friction_velocity(u_wall, y_p, nu);
        rho * u_tau * u_tau
    }

    /// Determine wall-function BC region: returns true if y⁺ is in the
    /// log-law region (11.53 ≤ y⁺ ≤ 300) where wall functions are valid.
    ///
    /// Wall functions are **not** appropriate in the viscous sublayer (y⁺ < 11.53)
    /// or the outer layer (y⁺ > 300) where the log-law breaks down.
    pub fn in_log_law_region<T: RealField + Copy + FromPrimitive + num_traits::Float>(y_plus: T) -> bool {
        let y_tr  = <T as FromPrimitive>::from_f64(WALL_Y_PLUS_TRANSITION).unwrap_or_else(T::one);
        let y_max = <T as FromPrimitive>::from_f64(300.0).unwrap_or_else(T::one);
        y_plus >= y_tr && y_plus <= y_max
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_viscous_sublayer() {
        // At y⁺ = 5 (viscous sublayer): u⁺ = y⁺ = 5.0
        let u_plus: f64 = WallFunctions::u_plus(5.0_f64);
        assert!((u_plus - 5.0).abs() < 1e-10, "Viscous sublayer: u⁺={u_plus:.4}, expected 5.0");
    }

    #[test]
    fn test_log_law_region() {
        // At y⁺ = 100: u⁺ = (1/0.41) ln(100) + 5.2 ≈ 18.44
        let u_plus: f64 = WallFunctions::u_plus(100.0_f64);
        let expected = f64::ln(100.0) / WALL_KAPPA + WALL_B;
        assert!((u_plus - expected).abs() < 1e-10, "Log law: u⁺={u_plus:.4}, expected {expected:.4}");
    }

    #[test]
    fn test_friction_velocity_round_trip() {
        // u_τ = 0.05 m/s, y = 0.001 m, ν = 1e-6 m²/s
        let nu = 1e-6_f64;
        let y  = 1e-3_f64;
        let u_tau_true = 0.05_f64;
        let y_plus_true = WallFunctions::y_plus(y, u_tau_true, nu);
        let u_plus_true = WallFunctions::u_plus(y_plus_true);
        let u_wall = u_tau_true * u_plus_true;

        let u_tau_est = WallFunctions::friction_velocity(u_wall, y, nu);
        assert!((u_tau_est - u_tau_true).abs() / u_tau_true < 1e-4,
            "u_τ round-trip: estimated={u_tau_est:.6}, true={u_tau_true}");
    }
}
