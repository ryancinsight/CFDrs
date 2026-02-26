//! Pressure-strain correlation model implementations.
//!
//! Three variants are provided, selected at runtime via [`super::wall_reflection::PressureStrainModel`]:
//! - **Linear** (Rotta, 1951): `Φ_ij = −C₁(ε/k) b_ij`
//! - **Quadratic** (Speziale et al., 1991): slow + rapid quadratic terms
//! - **SSG** (Speziale-Sarkar-Gatski, 1991): full non-linear formulation
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

pub mod linear;
pub mod quadratic;
pub mod ssg;

pub use linear::pressure_strain_linear;
pub use quadratic::pressure_strain_quadratic;
pub use ssg::pressure_strain_ssg;
