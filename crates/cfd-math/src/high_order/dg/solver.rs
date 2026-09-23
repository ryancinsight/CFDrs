//! DG time-dependent PDE solver.
//!
//! Combines a [`DGOperator`] with a [`TimeIntegrator`] to evolve DG
//! solutions forward in time using configurable explicit or implicit schemes.

//! This module root is a passthrough facade: the `DGSolver` carrier and its
//! inherent implementation live in `model`, with the test module beside it.

mod model;
mod time_integration;

pub use model::DGSolver;
pub use time_integration::*;

#[cfg(test)]
mod tests;
