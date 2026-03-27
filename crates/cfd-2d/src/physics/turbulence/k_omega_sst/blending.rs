//! SST blending functions and cross-diffusion computation.
//!
//! # SST Blending Functions (Menter 1994)
//!
//! ## Theorem: F1 Blending Function
//!
//! The primary blending function determines the transition between
//! k-ω (near-wall) and k-ε (freestream) behavior:
//!
//! ```text
//! F₁ = tanh(min[max(√k / (β* ω y), 500ν / (y² ω)), 4ρσ_ω₂ k / (CD_kω y²)]⁴)
//! ```
//!
//! where the cross-diffusion term is:
//! ```text
//! CD_kω = max(2ρ σ_ω₂ (1/ω) (∂k/∂xⱼ)(∂ω/∂xⱼ), 10⁻²⁰)
//! ```
//!
//! **Proof**: F₁ → 1 near the wall (activating k-ω constants) because
//! `√k / (β* ω y)` ~ `O(1)` in the viscous sublayer, while F₁ → 0 in the
//! freestream where `y → ∞`. The `tanh` ensures smooth transition.
//!
//! ## Theorem: F2 Blending Function
//!
//! ```text
//! F₂ = tanh[max(2√k / (β* ω y), 500ν / (y² ω))²]
//! ```
//!
//! F₂ controls the eddy viscosity limiter in the Bradshaw assumption region.
//!
//! # References
//! - Menter, F.R. (1994). AIAA Journal, 32(8), 1598-1605.

use super::super::constants::{OMEGA_MIN, SST_BETA_STAR, SST_SIGMA_OMEGA2};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Cross-diffusion term `CD_kω` for the SST blending function F1.
///
/// Evaluates `2 σ_ω₂ (∇k · ∇ω) / ω`, clamped to `≥ ω_min` for numerical stability.
pub fn cross_diffusion<T: RealField + FromPrimitive + Copy>(
    k: &[T],
    omega: &[T],
    idx: usize,
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> T {
    let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);
    let i = idx % nx;
    let j = idx / nx;

    if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
        return omega_min;
    }

    let two = T::from_f64(2.0).unwrap_or_else(T::one);

    // ∇k
    let dk_dx = (k[idx + 1] - k[idx - 1]) / (two * dx);
    let dk_dy = (k[idx + nx] - k[idx - nx]) / (two * dy);

    // ∇ω
    let domega_dx = (omega[idx + 1] - omega[idx - 1]) / (two * dx);
    let domega_dy = (omega[idx + nx] - omega[idx - nx]) / (two * dy);

    // ∇k · ∇ω
    let grad_dot = dk_dx * domega_dx + dk_dy * domega_dy;

    let sigma_omega2 = T::from_f64(SST_SIGMA_OMEGA2).unwrap_or_else(T::one);
    (two * sigma_omega2 * grad_dot / omega[idx].max(omega_min)).max(omega_min)
}

/// Compute SST blending functions F1 and F2 for all grid points.
///
/// F1 and F2 are stored in the provided mutable slices.
pub fn compute_blending_functions<T: RealField + FromPrimitive + Copy>(
    f1: &mut [T],
    f2: &mut [T],
    k: &[T],
    omega: &[T],
    wall_distance: &[T],
    molecular_viscosity: T,
    density: T,
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) {
    let beta_star = T::from_f64(SST_BETA_STAR).unwrap_or_else(T::one);
    let sigma_omega2 = T::from_f64(SST_SIGMA_OMEGA2).unwrap_or_else(T::one);
    let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);

    for idx in 0..k.len() {
        let y = wall_distance[idx];
        let k_val = k[idx].max(omega_min);
        let omega_val = omega[idx].max(omega_min);
        let nu = molecular_viscosity / density;

        let cd_kw = cross_diffusion(k, omega, idx, nx, ny, dx, dy);

        // F1 arguments
        let sqrt_k = k_val.sqrt();
        let arg1_1 = sqrt_k / (beta_star * omega_val * y);
        let arg1_2 = T::from_f64(500.0).unwrap_or_else(T::one) * nu / (y * y * omega_val);
        let arg1_3 = T::from_f64(4.0).unwrap_or_else(T::one) * density * sigma_omega2 * k_val
            / (cd_kw * y * y);
        let arg1 = arg1_1.min(arg1_2).max(arg1_3);
        f1[idx] = (T::from_f64(4.0).unwrap_or_else(T::one) * arg1).tanh();

        // F2 arguments
        let arg2_1 =
            T::from_f64(2.0).unwrap_or_else(T::one) * sqrt_k / (beta_star * omega_val * y);
        let arg2_2 = T::from_f64(500.0).unwrap_or_else(T::one) * nu / (y * y * omega_val);
        let arg2 = arg2_1.max(arg2_2);
        f2[idx] = (arg2 * arg2).tanh();
    }
}

/// Blend coefficient between k-ω (set 1) and k-ε (set 2) using F1.
///
/// `φ = F₁ · φ₁ + (1 − F₁) · φ₂`
#[inline]
pub fn blend_coefficient<T: RealField + FromPrimitive + Copy>(coef1: f64, coef2: f64, f1: T) -> T {
    let c1 = T::from_f64(coef1).unwrap_or_else(T::one);
    let c2 = T::from_f64(coef2).unwrap_or_else(T::one);
    f1 * c1 + (T::one() - f1) * c2
}
