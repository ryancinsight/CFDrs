//! Extended preconditioners for iterative linear solvers
//!
//! This module provides state-of-the-art preconditioning techniques
//! essential for practical CFD applications.

pub mod cholesky;
pub mod ilu;
pub mod multigrid;
pub mod ssor;

pub use cholesky::IncompleteCholesky;
pub use ilu::IncompleteLU;
pub use multigrid::AlgebraicMultigrid;
pub use ssor::SSOR;

// Edge case tests module for comprehensive coverage
#[cfg(test)]
mod edge_case_tests;
