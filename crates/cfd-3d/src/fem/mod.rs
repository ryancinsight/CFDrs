//! Finite Element Method (FEM) solver for 3D incompressible flows.
//!
//! This module implements a mixed finite element formulation for the Stokes
//! and Navier-Stokes equations with stabilization.
//!
//! ## Literature References
//! - Hughes, T.J.R., Franca, L.P., Balestra, M. (1986). "A new finite element formulation
//!   for computational fluid dynamics: V. Circumventing the Babuška-Brezzi condition:
//!   A stable Petrov-Galerkin formulation of the Stokes problem accommodating equal-order interpolations"
//! - Donea, J., Huerta, A. (2003). "Finite Element Methods for Flow Problems"

pub mod config;
pub mod constants;
pub mod element;
pub mod fluid;
pub mod problem;
pub mod solution;
pub mod solver;
pub mod stabilization;

// Re-export main types for convenience
pub use config::FemConfig;
pub use element::{ElementMatrices, FluidElement};
pub use fluid::FluidProperties;
pub use problem::StokesFlowProblem;
pub use solution::StokesFlowSolution;
pub use solver::FemSolver;
