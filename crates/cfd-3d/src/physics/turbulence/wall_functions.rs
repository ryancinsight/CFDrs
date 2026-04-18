//! Near-wall treatment using a smooth universal wall law.
//!
//! # Theorem — Spalding Wall Law
//!
//! Spalding's single-formula wall law blends the viscous sublayer, buffer
//! layer, and logarithmic region into one differentiable relation:
//!
//! ```text
//! y⁺ = u⁺ + e^(−κB) [e^(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6]
//! ```
//!
//! where `κ` is the von Kármán constant and `B` is the additive constant.
//! The friction velocity `u_τ` is recovered by solving the scalar nonlinear
//! equation `u_wall / u_τ = u⁺(y⁺)` with Newton's method and the exact
//! derivative of the Spalding profile.
//!
//! ## References
//!
//! - Spalding, D.B. (1961). "A single formula for the law of the wall."
//! - Recent wall-function formulations using the Spalding law appear in
//!   wall-modelled LES and RANS literature.

use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{WALL_B, WALL_KAPPA, WALL_Y_PLUS_TRANSITION};

/// Wall functions for near-wall turbulence treatment.
///
/// Provides `y⁺` calculation, a smooth wall-law inversion, and wall shear
/// stress estimation for RANS solvers.
pub struct WallFunctions;

impl WallFunctions {
    fn spalding_e<T: RealField + Copy + FromPrimitive + num_traits::Float>() -> T {
        let kappa = <T as FromPrimitive>::from_f64(WALL_KAPPA)
            .expect("WALL_KAPPA is an IEEE 754 representable f64 constant");
        let b = <T as FromPrimitive>::from_f64(WALL_B)
            .expect("WALL_B is an IEEE 754 representable f64 constant");
        num_traits::Float::exp(-(kappa * b))
    }

    /// Evaluate the Spalding wall law `y⁺(u⁺)`.
    fn spalding_y_plus<T: RealField + Copy + FromPrimitive + num_traits::Float>(u_plus: T) -> T {
        if u_plus <= T::zero() {
            return T::zero();
        }

        let kappa = <T as FromPrimitive>::from_f64(WALL_KAPPA)
            .expect("WALL_KAPPA is an IEEE 754 representable f64 constant");
        let e = Self::spalding_e::<T>();
        let ku = kappa * u_plus;
        let ku_sq = ku * ku;
        let ku_cu = ku_sq * ku;
        u_plus
            + e * (num_traits::Float::exp(ku)
                - T::one()
                - ku
                - ku_sq / (T::one() + T::one())
                - ku_cu / (T::one() + T::one() + T::one() + T::one() + T::one() + T::one()))
    }

    /// Evaluate `dy⁺/du⁺` for the Spalding wall law.
    fn spalding_dy_plus_du_plus<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        u_plus: T,
    ) -> T {
        if u_plus <= T::zero() {
            return T::one();
        }

        let kappa = <T as FromPrimitive>::from_f64(WALL_KAPPA)
            .expect("WALL_KAPPA is an IEEE 754 representable f64 constant");
        let e = Self::spalding_e::<T>();
        let ku = kappa * u_plus;
        let exp_ku = num_traits::Float::exp(ku);
        let kappa_sq = kappa * kappa;
        let kappa_cu = kappa_sq * kappa;
        T::one()
            + e * (kappa * exp_ku
                - kappa
                - kappa_sq * u_plus
                - (kappa_cu * u_plus * u_plus) / (T::one() + T::one()))
    }

    fn solve_u_plus_from_y_plus<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        y_plus: T,
    ) -> T {
        let eps = <T as FromPrimitive>::from_f64(1e-14)
            .expect("1e-14 is an IEEE 754 representable f64 constant");
        if y_plus <= eps {
            return T::zero();
        }

        let kappa = <T as FromPrimitive>::from_f64(WALL_KAPPA)
            .expect("WALL_KAPPA is an IEEE 754 representable f64 constant");
        let b = <T as FromPrimitive>::from_f64(WALL_B)
            .expect("WALL_B is an IEEE 754 representable f64 constant");
        let transition = <T as FromPrimitive>::from_f64(WALL_Y_PLUS_TRANSITION)
            .expect("WALL_Y_PLUS_TRANSITION is an IEEE 754 representable f64 constant");

        let mut u_plus = if y_plus < transition {
            y_plus
        } else {
            num_traits::Float::ln(y_plus) / kappa + b
        };
        u_plus = num_traits::Float::max(u_plus, eps);

        for _ in 0..32 {
            let f = Self::spalding_y_plus(u_plus) - y_plus;
            let df = Self::spalding_dy_plus_du_plus(u_plus);
            let du = f / df;
            u_plus -= du;
            u_plus = num_traits::Float::max(u_plus, eps);

            if num_traits::Float::abs(du) <= eps * num_traits::Float::max(T::one(), u_plus) {
                break;
            }
        }

        u_plus
    }

    /// Compute the dimensionless velocity `u⁺` from `y⁺`.
    pub fn u_plus<T: RealField + Copy + FromPrimitive + num_traits::Float>(y_plus: T) -> T {
        Self::solve_u_plus_from_y_plus(y_plus)
    }

    /// Compute dimensionless wall distance `y⁺ = y · u_τ / ν`.
    pub fn y_plus<T: RealField + Copy>(y: T, u_tau: T, nu: T) -> T {
        y * u_tau / nu
    }

    /// Estimate friction velocity `u_τ` from near-wall velocity using the
    /// Spalding wall law.
    pub fn friction_velocity<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        u_wall: T,
        y_p: T,
        nu: T,
    ) -> T {
        let eps = <T as FromPrimitive>::from_f64(1e-14)
            .expect("1e-14 is an IEEE 754 representable f64 constant");
        if u_wall <= eps || y_p <= eps || nu <= eps {
            return T::zero();
        }

        let u_tau_guess = num_traits::Float::sqrt(nu * u_wall / y_p);
        let mut u_plus = num_traits::Float::max(u_wall / u_tau_guess, eps);
        let target = y_p * u_wall / nu;

        for _ in 0..32 {
            let f = Self::spalding_y_plus(u_plus) - target / u_plus;
            let df = Self::spalding_dy_plus_du_plus(u_plus) + target / (u_plus * u_plus);
            let du = f / df;
            u_plus -= du;
            u_plus = num_traits::Float::max(u_plus, eps);

            if num_traits::Float::abs(du) <= eps * num_traits::Float::max(T::one(), u_plus) {
                break;
            }
        }

        let u_tau = u_wall / u_plus;
        num_traits::Float::max(u_tau, T::zero())
    }

    /// Compute wall shear stress `τ_w = ρ u_τ²`.
    pub fn wall_shear_stress<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        u_wall: T,
        y_p: T,
        nu: T,
        rho: T,
    ) -> T {
        let u_tau = Self::friction_velocity(u_wall, y_p, nu);
        rho * u_tau * u_tau
    }

    /// Determine whether `y⁺` lies in the traditional log-law region.
    pub fn in_log_law_region<T: RealField + Copy + FromPrimitive + num_traits::Float>(
        y_plus: T,
    ) -> bool {
        let y_tr = <T as FromPrimitive>::from_f64(WALL_Y_PLUS_TRANSITION)
            .expect("WALL_Y_PLUS_TRANSITION is an IEEE 754 representable f64 constant");
        let y_max = <T as FromPrimitive>::from_f64(300.0)
            .expect("300.0 is representable in all IEEE 754 types");
        y_plus >= y_tr && y_plus <= y_max
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_u_plus_round_trip() {
        let y_plus = 100.0_f64;
        let u_plus = WallFunctions::u_plus(y_plus);
        let y_plus_round_trip = super::WallFunctions::spalding_y_plus(u_plus);
        assert!(
            (y_plus_round_trip - y_plus).abs() < 1e-10,
            "round-trip through Spalding law must be exact to floating-point tolerance"
        );
    }

    #[test]
    fn test_friction_velocity_round_trip() {
        let nu = 1e-6_f64;
        let y = 1e-3_f64;
        let u_tau_true = 0.05_f64;
        let y_plus_true = WallFunctions::y_plus(y, u_tau_true, nu);
        let u_plus_true = WallFunctions::u_plus(y_plus_true);
        let u_wall = u_tau_true * u_plus_true;

        let u_tau_est = WallFunctions::friction_velocity(u_wall, y, nu);
        assert!(
            (u_tau_est - u_tau_true).abs() / u_tau_true < 1e-8,
            "u_τ round-trip: estimated={u_tau_est:.6}, true={u_tau_true}"
        );
    }

    #[test]
    fn test_wall_law_monotonicity() {
        let low = WallFunctions::u_plus(5.0_f64);
        let high = WallFunctions::u_plus(100.0_f64);
        assert!(high > low);
    }
}
