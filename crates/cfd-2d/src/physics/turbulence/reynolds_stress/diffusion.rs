//! Dissipation and turbulent transport terms for the RSTM.
//!
//! ### Dissipation Tensor ε_ij (Launder, 1975)
//! Either uses the anisotropic dissipation components stored in
//! [`ReynoldsStressTensor`] (when enabled), or falls back to the isotropic
//! approximation ε_ij = (2/3) ε δ_ij.
//!
//! ### Turbulent Transport T_ij (Launder et al., 1975)
//! Triple-correlation modelling via Daly & Harlow (1970):
//! `⟨u_i'u_j'u_k'⟩ ≈ −C_s (k³/ε²) ∂⟨u_i'u_j'⟩/∂x_k`

use super::tensor::ReynoldsStressTensor;
use nalgebra::RealField;
use num_traits::FromPrimitive;

fn c<T: RealField + Copy + FromPrimitive>(v: f64) -> T {
    T::from_f64(v).expect("diffusion constant must be representable")
}

/// Compute the dissipation tensor component `ε_ij` at grid point `(x, y)`.
///
/// Falls back to isotropic approximation when the optional anisotropic
/// components are not initialised.
#[inline]
pub fn dissipation_tensor<T: RealField + Copy + FromPrimitive>(
    rs: &ReynoldsStressTensor<T>,
    i: usize,
    j: usize,
    x: usize,
    y: usize,
) -> T {
    if let (Some(eps_xx), Some(eps_xy), Some(eps_yy)) =
        (&rs.epsilon_xx, &rs.epsilon_xy, &rs.epsilon_yy)
    {
        match (i, j) {
            (0, 0) => eps_xx[(x, y)],
            (0, 1) | (1, 0) => eps_xy[(x, y)],
            (1, 1) => eps_yy[(x, y)],
            _ => T::zero(),
        }
    } else {
        // Isotropic: ε_ij = (2/3) ε δ_ij
        let epsilon = rs.epsilon[(x, y)];
        match (i, j) {
            (0, 0) | (1, 1) => c::<T>(2.0 / 3.0) * epsilon,
            _ => T::zero(),
        }
    }
}

/// Optimised dissipation tensor (avoids Option unwrap for hot-path).
#[inline]
pub fn dissipation_tensor_optimized<T: RealField + Copy + FromPrimitive>(
    rs: &ReynoldsStressTensor<T>,
    i: usize,
    j: usize,
    x: usize,
    y: usize,
    two_thirds: T,
    epsilon: T,
) -> T {
    if let (Some(eps_xx), Some(eps_xy), Some(eps_yy)) =
        (&rs.epsilon_xx, &rs.epsilon_xy, &rs.epsilon_yy)
    {
        match (i, j) {
            (0, 0) => eps_xx[(x, y)],
            (0, 1) | (1, 0) => eps_xy[(x, y)],
            (1, 1) => eps_yy[(x, y)],
            _ => T::zero(),
        }
    } else {
        match (i, j) {
            (0, 0) | (1, 1) => two_thirds * epsilon,
            _ => T::zero(),
        }
    }
}

/// Compute turbulent transport `T_ij` via triple-correlation modelling.
///
/// `T_ij = −∂⟨u_i'u_j'u_k'⟩/∂x_k ≈ C_s (k³/ε²) ∇²⟨u_i'u_j'⟩`
#[inline]
pub fn turbulent_transport<T: RealField + Copy + FromPrimitive>(
    k: T,
    epsilon: T,
    stress_gradient: &[[T; 2]; 2],
    i: usize,
    j: usize,
) -> T {
    let c_s = c::<T>(0.11);
    let diffusion_coeff = c_s * k * k * k / (epsilon * epsilon);

    match (i, j) {
        (0, 0) => -diffusion_coeff * stress_gradient[0][0],
        (0, 1) | (1, 0) => {
            -diffusion_coeff * c::<T>(0.5) * (stress_gradient[0][1] + stress_gradient[1][0])
        }
        (1, 1) => -diffusion_coeff * stress_gradient[1][1],
        _ => T::zero(),
    }
}
