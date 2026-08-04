//! Core operators for Discontinuous Galerkin methods.
//!
//! This module implements the spatial discretization operators for DG methods,
//! including volume and surface integrals, weak form operators, and boundary conditions.

mod dg_operator;
mod params;

pub use dg_operator::DGOperator;
pub use params::{BoundaryCondition, BoundaryFlux, DGOperatorParams};

#[cfg(test)]
mod tests;
