//! Pressure-strain correlation model implementations.
//!
//! Three variants are provided, selected at runtime via [`super::wall_reflection::PressureStrainModel`]:
//! - **Linear** (Rotta, 1951): `Φ_ij = −C₁(ε/k) b_ij`
//! - **Quadratic** (Speziale et al., 1991): slow + rapid quadratic terms
//! - **SSG** (Speziale-Sarkar-Gatski, 1991): full non-linear formulation

pub mod linear;
pub mod quadratic;
pub mod ssg;

pub use linear::pressure_strain_linear;
pub use quadratic::pressure_strain_quadratic;
pub use ssg::pressure_strain_ssg;
