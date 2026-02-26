//! Stress and strain tensor operations for Finite Element Methods.
//!
//! Provides helper methods for computing stress and strain tensors from
//! velocity gradients and fluid properties, specifically for FEM assembly.
//!
//! # Theorem — Newtonian Constitutive Relation (Cauchy 1827, Stokes 1845)
//!
//! For an incompressible Newtonian fluid, the Cauchy stress tensor is
//!
//! ```text
//! σ = −p I + 2μ ε̇
//! ```
//!
//! where $p$ is the mechanical pressure, $\mu$ is the dynamic viscosity, and
//! the strain-rate tensor is the symmetric part of the velocity gradient:
//!
//! ```text
//! ε̇ = ½(∇u + (∇u)ᵀ)
//! ```
//!
//! **Theorem (Symmetry).** $\dot{\varepsilon}$ is symmetric by construction:
//! $\dot{\varepsilon}_{ij} = \dot{\varepsilon}_{ji}$.
//!
//! **Theorem (Trace-free for incompressible flow).** $\mathrm{tr}(\dot{\varepsilon})
//! = \nabla \cdot \mathbf{u} = 0$, so the deviatoric and total strain rates coincide.

use cfd_core::physics::fluid::ConstantPropertyFluid;
use nalgebra::{Matrix3, RealField};
use num_traits::FromPrimitive;

/// Calculate stress tensor from strain rate and pressure.
///
/// $\sigma = -pI + 2\mu \dot{\epsilon}$
pub fn stress_tensor<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    fluid: &ConstantPropertyFluid<T>,
    pressure: T,
    strain_rate: &Matrix3<T>,
) -> Matrix3<T> {
    let two = <T as FromPrimitive>::from_f64(2.0).unwrap_or_else(T::zero);
    let mut stress = strain_rate * (two * fluid.viscosity);

    // Add pressure contribution to diagonal
    stress[(0, 0)] -= pressure;
    stress[(1, 1)] -= pressure;
    stress[(2, 2)] -= pressure;

    stress
}

/// Calculate strain rate tensor from velocity gradient.
///
/// $\dot{\epsilon} = 0.5 (\nabla u + (\nabla u)^T)$
pub fn strain_rate_tensor<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(velocity_gradient: &Matrix3<T>) -> Matrix3<T> {
    let half = <T as FromPrimitive>::from_f64(0.5).unwrap_or_else(T::zero);
    (velocity_gradient + velocity_gradient.transpose()) * half
}
