//! Reynolds Stress Transport Model (RSTM) — module hierarchy.
//!
//! ## Structure
//! ```text
//! reynolds_stress/
//!   mod.rs            — PressureStrainModel enum + canonical pub use re-exports
//!   tensor.rs         — ReynoldsStressTensor storage type
//!   model.rs          — ReynoldsStressModel struct + constructor + initialisation
//!   production.rs     — P_ij exact production term (not Boussinesq)
//!   diffusion.rs      — ε_ij dissipation tensor + T_ij turbulent transport
//!   curvature.rs      — Suga-Craft (2003) streamline curvature correction
//!   wall_reflection.rs — Gibson-Launder (1978) wall-reflection correction
//!   pressure_strain/
//!     linear.rs       — Rotta (1951) linear return-to-isotropy
//!     quadratic.rs    — Speziale et al. (1991) quadratic model
//!     ssg.rs          — Full SSG model
//!   transport.rs      — Full time-advancement + TurbulenceModel trait impl
//! ```
//!
//! ## References
//! - Pope, S. B. (2000). *Turbulent Flows*. Cambridge University Press.
//! - Launder, B. E., Reece, G. J., & Rodi, W. (1975). J. Fluid Mech., 68(3), 537–566.
//! - Speziale, C. G., Sarkar, S., & Gatski, T. B. (1991). J. Fluid Mech., 227, 245–272.
//! - Gibson, M. M., & Launder, B. E. (1978). J. Fluid Mech., 86(3), 491–511.
//! - Suga, K., & Craft, T. J. (2003). Flow Turbulence Combust., 70, 143–162.
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

pub mod curvature;
pub mod diffusion;
pub mod model;
pub mod pressure_strain;
pub mod production;
pub mod tensor;
pub mod transport;
pub mod wall_reflection;

pub use model::ReynoldsStressModel;
pub use tensor::ReynoldsStressTensor;

/// Pressure-strain correlation model selector.
///
/// Passed to [`ReynoldsStressModel`] at construction time and used at every
/// pressure-strain evaluation to dispatch to the correct kernel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PressureStrainModel {
    /// Linear return-to-isotropy (Rotta, 1951). Φ_ij = −C₁(ε/k) b_ij.
    LinearReturnToIsotropy,
    /// Quadratic slow + rapid (Speziale, Sarkar & Gatski, 1991).
    Quadratic,
    /// Full SSG model: non-linear in b_ij with strain and vorticity coupling.
    SSG,
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
