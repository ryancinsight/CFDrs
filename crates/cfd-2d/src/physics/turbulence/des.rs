//! Detached Eddy Simulation (DES) model
//!
//! DES is a hybrid RANS-LES approach where RANS models are used near walls
//! and LES is used in the detached regions away from walls. This provides
//! accurate boundary layer prediction with reduced computational cost.
//!
//! ## DES97 Formulation
//!
//! The DES length scale is defined as:
//!
//! L_DES = min(L_RANS, C_DES * Δ)
//!
//! where L_RANS is the RANS length scale, C_DES ≈ 0.65, and Δ is the grid spacing.
//!
//! ## DDES (Delayed DES)
//!
//! Delayed DES prevents premature switching to LES in boundary layers by
//! modifying the length scale computation to account for the wall distance.
//!
//! ## IDDES (Improved DDES)
//!
//! IDDES further improves the shielding function and adds a wall-modeled LES
//! capability for very fine grids near walls.
//!
//! ## References
//!
//! - Spalart, P. R., et al. (1997). Comments on the feasibility of LES for wings.
//! - Spalart, P. R., et al. (2006). A new version of detached-eddy simulation.
//! - Shur, M. L., et al. (2008). A hybrid RANS-LES approach with delayed-DES.
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

// Dependencies: RANS model integration, wall distance computation, length scale consistency
// Mathematical Foundation: Spalart et al. (1997) DES97, Shur et al. (2008) IDDES

mod config;
mod length_scale;
mod les;
mod model;

pub use config::{DESConfig, DESVariant};
pub use model::DetachedEddySimulation;

#[cfg(test)]
mod tests;
