//! Linear system solvers for CFD applications.
//!
//! This module provides iterative linear solvers for sparse matrices
//! arising from CFD discretizations. It emphasizes the use of validated
//! implementations and efficient memory management.
//!
//! ## Design Philosophy
//!
//! Rather than re-implementing complex numerical algorithms from scratch, this module
//! focuses on providing efficient wrappers and using established implementations where
//! available. This approach reduces maintenance overhead, eliminates numerical bugs,
//! and provides performance through validated algorithms.
//!
//! ## Performance Considerations
//!
//! All solvers prioritize memory efficiency by:
//! - Using in-place operations to avoid unnecessary allocations
//! - Pre-allocating workspace vectors outside iteration loops
//! - Leveraging nalgebra's BLAS-like operations
//! - Providing efficient preconditioner APIs

mod bicgstab;
mod conjugate_gradient;
mod preconditioners;
mod traits;

pub use bicgstab::BiCGSTAB;
pub use conjugate_gradient::ConjugateGradient;
pub use preconditioners::{IdentityPreconditioner, JacobiPreconditioner, SORPreconditioner};
pub use traits::{LinearSolver, Preconditioner};

// Re-export the unified configuration from cfd-core
pub use cfd_core::solver::SolverConfiguration;
pub use cfd_core::solvers::LinearSolverConfig;

#[cfg(test)]
mod tests;
