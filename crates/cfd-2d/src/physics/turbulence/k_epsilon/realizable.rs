//! Realizable k-ε variant: strain-rate-dependent C_μ.
//!
//! ## Theorem — Realizable C_μ (Shih, Zhu & Lumley 1995)
//!
//! The standard k-ε model uses a fixed C_μ = 0.09. The Realizable variant
//! makes C_μ a function of the mean strain rate:
//!
//! ```text
//! C_μ = 1 / (A₀ + A_s · S̃ · k / ε)
//! ```
//!
//! where:
//! - A₀ = 4.04 (calibrated constant ensuring C_μ ≤ 1/4.04 ≈ 0.247)
//! - A_s = √6 · cos(φ/3), φ = ⅓ arccos(√6 · W)
//! - W = S_ij S_jk S_ki / S̃³ (third invariant of the strain tensor)
//! - S̃ = √(2 S_ij S_ij) (strain rate magnitude)
//!
//! **Proof sketch**: The denominator A₀ + A_s·S̃·k/ε ≥ A₀ > 0 by
//! construction, so C_μ ∈ (0, 1/A₀]. For S̃ → 0, C_μ → 1/A₀ ≈ 0.247.
//! For S̃ → ∞, C_μ → 0, preventing stagnation-region over-prediction.
//!
//! ## Reference
//!
//! Shih, T.-H., Zhu, J., & Lumley, J. L. (1995). A New Reynolds Stress
//! Algebraic Equation Model. *Computers & Fluids*, 24(3), 227–238.

use super::model::KEpsilonModel;
use crate::physics::turbulence::constants::{EPSILON_MIN, REALIZABLE_A0};
use cfd_core::physics::constants::mathematical::numeric::{ONE_HALF, TWO};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Compute the Realizable C_μ from local strain rate, k, and ε.
///
/// Returns the local realizable C_μ value, bounded in (0, 1/A₀].
pub fn realizable_c_mu<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive>(
    model: &KEpsilonModel<T>,
    velocity_gradient: &[[T; 2]; 2],
    k: T,
    epsilon: T,
) -> T {
    let _ = model; // model fields not needed beyond this point; API future-proofing.
    let eps_safe = epsilon.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero));
    let half = T::from_f64(ONE_HALF).unwrap_or_else(T::one);

    // Build symmetric strain rate tensor S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
    let s00 = velocity_gradient[0][0];
    let s11 = velocity_gradient[1][1];
    let s01 = (velocity_gradient[0][1] + velocity_gradient[1][0]) * half;

    // S̃ = sqrt(2 * S_ij * S_ij)
    let two = T::from_f64(TWO).unwrap_or_else(T::one);
    let s_sq = s00 * s00 + s11 * s11 + two * s01 * s01;
    let s_tilde = (two * s_sq).sqrt();

    // Compute W = S_ij * S_jk * S_ki / S̃³  (third invariant)
    let ss00 = s00 * s00 + s01 * s01;
    let ss01 = s00 * s01 + s01 * s11;
    let ss11 = s01 * s01 + s11 * s11;

    let trace_s3 = ss00 * s00 + two * ss01 * s01 + ss11 * s11;

    let s_tilde_cubed = s_tilde * s_tilde * s_tilde;
    let s_tilde_min = T::from_f64(1e-30).unwrap_or_else(T::zero);

    let w = if s_tilde_cubed > s_tilde_min {
        trace_s3 / s_tilde_cubed
    } else {
        T::zero()
    };

    // φ = (1/3) * arccos(√6 * W)
    let sqrt6 = T::from_f64(6.0_f64.sqrt()).unwrap_or_else(T::one);
    let one = T::one();
    let arg = (sqrt6 * w).max(-one).min(one);
    let phi = arg.acos();
    let third = T::from_f64(1.0 / 3.0).unwrap_or_else(T::one);

    // A_s = √6 * cos(φ / 3)
    let a_s = sqrt6 * (phi * third).cos();

    // C_μ = 1 / (A₀ + A_s * S̃ * k / ε)
    let a0 = T::from_f64(REALIZABLE_A0).unwrap_or_else(T::one);
    let denom = a0 + a_s * s_tilde * k / eps_safe;

    // Ensure denominator is at least A₀ (C_μ ≤ 1/A₀) and positive
    let denom_safe = denom.max(a0);
    one / denom_safe
}
