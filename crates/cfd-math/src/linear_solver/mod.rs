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
mod config;
mod conjugate_gradient;
mod traits;
pub mod operators;
pub mod gmres;
pub mod matrix_free;
pub mod preconditioners;

pub use bicgstab::BiCGSTAB;
pub use config::IterativeSolverConfig;
pub use conjugate_gradient::ConjugateGradient;
pub use gmres::GMRES;
#[cfg(feature = "mpi")]
pub use matrix_free::{
    CommunicationOptimization, CommunicationOptimizer, CommunicationOverlap,
    LoadBalancingRecommendations, LoadBalancingStrategy, ParallelLoadBalancer,
    ParallelMatrixFreeBiCGSTAB,
};
pub use matrix_free::{
    EnergyOperator2D, LaplacianOperator2D, LinearOperator,
    MomentumOperator1D, MomentumOperator2D, PoissonOperator3D,
};

pub use preconditioners::{
    AlgebraicMultigrid, DeflationPreconditioner, IdentityPreconditioner, IncompleteLU,
    JacobiPreconditioner, SORPreconditioner, SerialSchwarzPreconditioner, SSOR,
};
pub use preconditioners::multigrid::AMGConfig;

pub use traits::{Configurable, IterativeLinearSolver, LinearSolver, Preconditioner};

// REMOVED: Dependencies on cfd-core break the dependency hierarchy.
// cfd-math should be a foundational library that doesn't depend on application code.
// Linear solvers should define their own configuration types.

#[cfg(test)]
mod tests;

// Edge case tests module for comprehensive coverage
#[cfg(test)]
mod edge_case_tests;

// Extended edge case tests with property-based testing
#[cfg(test)]
mod extended_edge_case_tests;
