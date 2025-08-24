//! Linear system solvers for CFD applications.
//!
//! This module provides iterative linear solvers optimized for sparse matrices
//! arising from CFD discretizations. It emphasizes the use of well-tested
//! implementations and efficient memory management.
//!
//! ## Design Philosophy
//! 
//! Rather than re-implementing complex numerical algorithms from scratch, this module
//! focuses on providing efficient wrappers and using robust implementations where
//! available. This approach reduces maintenance overhead, eliminates numerical bugs,
//! and provides better performance through optimized algorithms.
//!
//! ## Performance Considerations
//!
//! All solvers prioritize memory efficiency by:
//! - Using in-place operations to avoid unnecessary allocations
//! - Pre-allocating workspace vectors outside iteration loops
//! - Leveraging nalgebra's optimized BLAS-like operations
//! - Providing efficient preconditioner APIs

mod traits;
mod preconditioners;
mod conjugate_gradient;
mod bicgstab;

pub use traits::{LinearSolver, Preconditioner};
pub use preconditioners::{IdentityPreconditioner, JacobiPreconditioner, SORPreconditioner};
pub use conjugate_gradient::ConjugateGradient;
pub use bicgstab::BiCGSTAB;

// Re-export the unified configuration from cfd-core
pub use cfd_core::solver::{LinearSolverConfig, SolverConfiguration};

#[cfg(test)]
mod tests;