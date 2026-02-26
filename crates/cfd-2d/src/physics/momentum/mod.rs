//! Momentum equation solver with clean domain separation
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

mod boundary;
mod coefficient_corrections;
mod coefficients;
mod discretization;
mod interpolation;
pub mod muscl;
mod setup;
mod solver;
pub mod tvd_limiters;

#[cfg(test)]
mod tvd_limiter_edge_cases;

#[cfg(test)]
mod algorithm_validation_tests;

pub use boundary::apply_momentum_boundaries;
pub use coefficients::{ConvectionScheme, MomentumCoefficients};
pub use discretization::{CentralDifference, DiscretizationScheme, MusclDiscretization, Upwind};
pub use interpolation::rhie_chow_interpolation;
pub use muscl::{schemes, MusclOrder, MusclReconstruction, MusclScheme};
pub use setup::BoundarySetup;
pub use solver::{MomentumComponent, MomentumSolver};
pub use tvd_limiters::{Minmod, MonotonizedCentral, Superbee, TvdLimiter, VanLeer};
