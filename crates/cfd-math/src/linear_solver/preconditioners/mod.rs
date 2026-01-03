//! Preconditioners for iterative linear solvers
//!
//! This module provides various preconditioning techniques essential for 
//! accelerating the convergence of iterative solvers in CFD applications.

pub mod basic;
pub mod cholesky;
pub mod deflation;
pub mod ilu;
pub mod multigrid;
pub mod schwarz;
pub mod ssor;

pub use basic::{IdentityPreconditioner, JacobiPreconditioner, SORPreconditioner};
pub use cholesky::IncompleteCholesky;
pub use deflation::DeflationPreconditioner;
pub use ilu::IncompleteLU;
pub use multigrid::AlgebraicMultigrid;
pub use schwarz::SerialSchwarzPreconditioner;
pub use ssor::SSOR;

// Edge case tests module for comprehensive coverage
#[cfg(test)]
mod edge_case_tests;
