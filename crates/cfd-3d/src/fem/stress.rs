//! Stress and strain tensor operations for Finite Element Methods.
//!
//! Provides helper methods for computing stress and strain tensors from
//! velocity gradients and fluid properties, specifically for FEM assembly.
//!
//! # Theorem — Generalised-Newtonian Constitutive Relation
//!
//! For an incompressible generalised-Newtonian fluid, the Cauchy stress tensor is
//!
//! ```text
//! σ = −p I + 2μ(γ̇) ε̇
//! ```
//!
//! where $p$ is the mechanical pressure and the apparent viscosity $\mu(\dot{\gamma})$
//! depends on the scalar shear rate $\dot{\gamma} = \sqrt{2\,\dot{\varepsilon}:\dot{\varepsilon}}$.
//! The strain-rate tensor is the symmetric part of the velocity gradient:
//!
//! ```text
//! ε̇ = ½(∇u + (∇u)ᵀ)
//! ```
//!
//! The Newtonian case is recovered when $\mu(\dot{\gamma}) = \mu_0 = \text{const}$.
//!
//! For the Carreau-Yasuda model (Yasuda 1979):
//!
//! ```text
//! μ(γ̇) = μ_∞ + (μ_0 − μ_∞)(1 + (λγ̇)^a)^{(n−1)/a}
//! ```
//!
//! **Theorem (Symmetry).** $\dot{\varepsilon}$ is symmetric by construction:
//! $\dot{\varepsilon}_{ij} = \dot{\varepsilon}_{ji}$.
//!
//! **Theorem (Trace-free for incompressible flow).** $\mathrm{tr}(\dot{\varepsilon})
//! = \nabla \cdot \mathbf{u} = 0$, so the deviatoric and total strain rates coincide.
//!
//! **References:**
//! - Cauchy, A.-L. (1827). Mémoire sur les dilatations, les condensations et les
//!   rotations produites par un changement de forme dans un système de points matériels.
//! - Stokes, G.G. (1845). On the theories of the internal friction of fluids in motion.
//! - Yasuda, K. (1979). Investigation of the analogies between viscometric and
//!   linear viscoelastic properties of polystyrene fluids. PhD thesis, MIT.

use cfd_core::physics::fluid::ConstantPropertyFluid;
use nalgebra::{Matrix3, RealField};
use num_traits::FromPrimitive;

/// Calculate the Cauchy stress tensor from a scalar viscosity, pressure,
/// and strain-rate tensor.
///
/// # Generalised-Newtonian Constitutive Relation
///
/// ```text
/// σ = −p I + 2μ ε̇
/// ```
///
/// The viscosity `mu` may be spatially varying (e.g. from a Picard iteration
/// over a Carreau-Yasuda model) or constant.
///
/// # Arguments
/// * `mu` — Dynamic viscosity [Pa·s] (scalar, possibly element-local).
/// * `pressure` — Mechanical pressure [Pa].
/// * `strain_rate` — Symmetric strain-rate tensor ε̇ = ½(∇u + (∇u)ᵀ).
pub fn stress_tensor<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    mu: T,
    pressure: T,
    strain_rate: &Matrix3<T>,
) -> Matrix3<T> {
    let two =
        <T as FromPrimitive>::from_f64(2.0).expect("2.0 is representable in all IEEE 754 types");
    let mut stress = strain_rate * (two * mu);

    // Add pressure contribution to diagonal: σ_ij += −p δ_ij
    stress[(0, 0)] -= pressure;
    stress[(1, 1)] -= pressure;
    stress[(2, 2)] -= pressure;

    stress
}

/// Calculate the Cauchy stress tensor using a `ConstantPropertyFluid`.
///
/// Convenience wrapper over [`stress_tensor`] that extracts
/// `fluid.viscosity` as the scalar viscosity parameter.
pub fn stress_tensor_with_fluid<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    fluid: &ConstantPropertyFluid<T>,
    pressure: T,
    strain_rate: &Matrix3<T>,
) -> Matrix3<T> {
    stress_tensor(fluid.viscosity, pressure, strain_rate)
}

/// Calculate strain rate tensor from velocity gradient.
///
/// ```text
/// ε̇ = ½(∇u + (∇u)ᵀ)
/// ```
///
/// # Theorem (Symmetry)
///
/// The output is symmetric by construction: ε̇_{ij} = ε̇_{ji}.
pub fn strain_rate_tensor<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    velocity_gradient: &Matrix3<T>,
) -> Matrix3<T> {
    let half =
        <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");
    (velocity_gradient + velocity_gradient.transpose()) * half
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Matrix3;

    /// Simple shear flow: u = (γ̇·y, 0, 0).
    /// ∇u has only (0,1) entry = γ̇.
    /// ε̇ = ½[[0, γ̇, 0],[γ̇, 0, 0],[0, 0, 0]]
    /// σ = −pI + 2μ ε̇
    ///
    /// Analytical result for p=100, μ=3.5e-3, γ̇=200 s⁻¹:
    ///   σ_01 = σ_10 = 2 × 3.5e-3 × 0.5 × 200 = 0.7 Pa
    ///   σ_00 = σ_11 = σ_22 = −100 Pa
    #[test]
    fn stress_tensor_simple_shear_analytically_derived() {
        let gamma_dot: f64 = 200.0;
        let mu: f64 = 3.5e-3; // blood-like viscosity
        let pressure: f64 = 100.0;

        // Velocity gradient for simple shear: du_x/dy = γ̇
        let grad_u = Matrix3::new(
            0.0, gamma_dot, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        );
        let eps = strain_rate_tensor(&grad_u);
        let sigma = stress_tensor(mu, pressure, &eps);

        // Diagonal: −p
        assert!((sigma[(0, 0)] - (-pressure)).abs() < 1e-12);
        assert!((sigma[(1, 1)] - (-pressure)).abs() < 1e-12);
        assert!((sigma[(2, 2)] - (-pressure)).abs() < 1e-12);

        // Off-diagonal shear stress: 2μ × ε̇_01 = 2 × 3.5e-3 × (200/2) = 0.7
        let expected_shear = 2.0 * mu * (gamma_dot / 2.0);
        assert!(
            (sigma[(0, 1)] - expected_shear).abs() < 1e-12,
            "σ_01 = {}, expected {expected_shear}",
            sigma[(0, 1)]
        );
        assert!(
            (sigma[(1, 0)] - expected_shear).abs() < 1e-12,
            "σ_10 = {}, expected {expected_shear}",
            sigma[(1, 0)]
        );

        // Zero off-diagonal entries
        assert!((sigma[(0, 2)]).abs() < 1e-12);
        assert!((sigma[(2, 0)]).abs() < 1e-12);
        assert!((sigma[(1, 2)]).abs() < 1e-12);
        assert!((sigma[(2, 1)]).abs() < 1e-12);
    }

    /// Verify that `stress_tensor_with_fluid` produces identical output
    /// to `stress_tensor(fluid.viscosity, ...)`.
    #[test]
    fn stress_tensor_with_fluid_matches_scalar_api() {
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let pressure = 50.0;
        let grad_u = Matrix3::new(
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        );
        let eps = strain_rate_tensor(&grad_u);
        let sigma_scalar = stress_tensor(fluid.viscosity, pressure, &eps);
        let sigma_fluid = stress_tensor_with_fluid(&fluid, pressure, &eps);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (sigma_scalar[(i, j)] - sigma_fluid[(i, j)]).abs() < 1e-15,
                    "Mismatch at ({i},{j}): {} vs {}",
                    sigma_scalar[(i, j)],
                    sigma_fluid[(i, j)]
                );
            }
        }
    }

    /// Strain rate tensor must be symmetric for an arbitrary velocity gradient.
    #[test]
    fn strain_rate_tensor_symmetry() {
        let grad_u: Matrix3<f64> = Matrix3::new(
            1.0, 3.0, 5.0,
            2.0, 4.0, 7.0,
            9.0, 6.0, 8.0,
        );
        let eps = strain_rate_tensor(&grad_u);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (eps[(i, j)] - eps[(j, i)]).abs() < 1e-15,
                    "ε̇ not symmetric at ({i},{j}): {} ≠ {}",
                    eps[(i, j)],
                    eps[(j, i)]
                );
            }
        }
    }

    /// For pure shear flow, tr(ε̇) = 0 (incompressibility).
    /// u = (γ̇ y, 0, 0) ⟹ ∇·u = 0, tr(ε̇) = 0.
    #[test]
    fn strain_rate_tensor_trace_free_for_incompressible_shear() {
        let gamma_dot: f64 = 500.0;
        let grad_u: Matrix3<f64> = Matrix3::new(
            0.0, gamma_dot, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        );
        let eps = strain_rate_tensor(&grad_u);
        let trace = eps[(0, 0)] + eps[(1, 1)] + eps[(2, 2)];
        assert!(
            trace.abs() < 1e-15,
            "tr(ε̇) = {trace}, expected 0 for incompressible shear"
        );
    }
}

