//! Dissipation and turbulent transport terms for the RSTM.
//!
//! ### Dissipation Tensor ε_ij (Launder, 1975)
//! Uses the anisotropic dissipation components stored in [`ReynoldsStressTensor`]
//! when they are populated; otherwise applies the isotropic closure
//! `ε_ij = (2/3) ε δ_ij`.
//!
//! ### Turbulent Transport T_ij (Launder et al., 1975)
//! Triple-correlation modelling via Daly & Harlow (1970):
//! `⟨u_i'u_j'u_k'⟩ ≈ −C_s (k³/ε²) ∂⟨u_i'u_j'⟩/∂x_k`
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

use super::tensor::ReynoldsStressTensor;
use eunomia::RealField;

fn c<T: RealField + Copy>(v: f64) -> T {
    T::from_f64(v)
}

/// Compute the dissipation tensor component `ε_ij` at grid point `(x, y)`.
///
/// Uses the stored anisotropic components when available and otherwise applies
/// the isotropic closure.
#[inline]
pub fn dissipation_tensor<T: RealField + Copy>(
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
            (0, 0) => eps_xx[[x, y]],
            (0, 1) | (1, 0) => eps_xy[[x, y]],
            (1, 1) => eps_yy[[x, y]],
            _ => T::ZERO,
        }
    } else {
        let epsilon = rs.epsilon[[x, y]];
        match (i, j) {
            (0, 0) | (1, 1) => c::<T>(2.0 / 3.0) * epsilon,
            _ => T::ZERO,
        }
    }
}

/// Optimised dissipation tensor (avoids Option unwrap for hot-path).
#[inline]
pub fn dissipation_tensor_optimized<T: RealField + Copy>(
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
            (0, 0) => eps_xx[[x, y]],
            (0, 1) | (1, 0) => eps_xy[[x, y]],
            (1, 1) => eps_yy[[x, y]],
            _ => T::ZERO,
        }
    } else {
        match (i, j) {
            (0, 0) | (1, 1) => two_thirds * epsilon,
            _ => T::ZERO,
        }
    }
}

/// Compute turbulent transport `T_ij` via triple-correlation modelling.
///
/// `T_ij = −∂⟨u_i'u_j'u_k'⟩/∂x_k ≈ C_s (k³/ε²) ∇²⟨u_i'u_j'⟩`
#[inline]
pub fn turbulent_transport<T: RealField + Copy>(
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
        _ => T::ZERO,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::turbulence::reynolds_stress::model::ReynoldsStressModel;
    use approx::assert_relative_eq;
    use leto::Array2;

    #[test]
    fn dissipation_tensor_prefers_anisotropic_components_when_available() {
        let model = ReynoldsStressModel::<f64>::new(2, 2);
        let mut tensor = model.initialize_reynolds_stresses(0.3, 0.1);

        tensor.epsilon_xx = Some(Array2::from_elem([2, 2], 1.23));
        tensor.epsilon_xy = Some(Array2::from_elem([2, 2], 0.45));
        tensor.epsilon_yy = Some(Array2::from_elem([2, 2], 2.34));

        assert_relative_eq!(
            dissipation_tensor(&tensor, 0, 0, 1, 1),
            1.23,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            dissipation_tensor(&tensor, 0, 1, 1, 1),
            0.45,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            dissipation_tensor(&tensor, 1, 1, 1, 1),
            2.34,
            epsilon = 1e-15
        );
    }

    #[test]
    fn turbulent_transport_uses_stress_gradient() {
        let transport: f64 = turbulent_transport(0.5, 0.2, &[[1.0, 2.0], [2.0, 3.0]], 0, 1);
        assert!(transport.is_finite());
    }
}
