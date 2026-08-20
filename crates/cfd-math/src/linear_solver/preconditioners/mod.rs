//! CFD-specific preconditioners.
//!
//! Athena owns the preconditioner seam (Atlas ADR 0033): a preconditioner is
//! anything implementing [`athena_core::Preconditioner`] for the backend the
//! solve runs on. The general-purpose implementations - identity, Jacobi,
//! incomplete LU, successive over-relaxation - come from `athena-leto` and are
//! used directly at their canonical paths rather than re-exported here.
//!
//! What lives in this module is CFD-domain-specific: algebraic multigrid and
//! the level-k incomplete LU built on the CFD sparsity patterns.

pub mod ilu;
pub mod multigrid;

#[cfg(test)]
mod ssor_tests;

pub use ilu::IncompleteLU;
pub use multigrid::AlgebraicMultigrid;
