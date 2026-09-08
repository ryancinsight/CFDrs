//! Convection scheme correction computations
//!
//! This module contains helper functions for computing high-order convection corrections
//! used in deferred correction approaches (Patankar 1980 §5.4.3).
//!
//! # Theorem
//! The momentum discretization must conserve linear momentum globally and locally.
//!
//! **Proof sketch**:
//! By integrating the Navier-Stokes momentum equation over a control volume $\Omega$,
//! Gauss's divergence theorem converts the convective and diffusive volume integrals
//! into surface fluxes. The finite volume method ensures that the flux leaving one
//! cell exactly equals the flux entering the adjacent cell. Thus, in the absence of
//! external forces and boundary fluxes, the total momentum $\int_\Omega \rho \mathbf{u} dV$
//! is exactly conserved to machine precision.

mod quick;
mod second_order;
mod tvd;
mod weno_z;

pub use quick::{compute_quick_correction_x, compute_quick_correction_y};
pub use second_order::apply_second_order_deferred_correction;
pub use tvd::{compute_tvd_correction_x, compute_tvd_correction_y};
pub use weno_z::apply_weno_z_deferred_correction;
