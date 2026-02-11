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
pub mod block_preconditioner;
mod direct_solver;
mod config;
mod conjugate_gradient;
pub mod gmres;
pub mod matrix_free;
pub mod operators;
pub mod preconditioners;
mod traits;

pub use bicgstab::BiCGSTAB;
pub use config::IterativeSolverConfig;
pub use conjugate_gradient::ConjugateGradient;
pub use direct_solver::DirectSparseSolver;
pub use gmres::GMRES;
#[cfg(feature = "mpi")]
pub use matrix_free::{
    CommunicationOptimization, CommunicationOptimizer, CommunicationOverlap,
    LoadBalancingRecommendations, LoadBalancingStrategy, ParallelLoadBalancer,
    ParallelMatrixFreeBiCGSTAB,
};
pub use matrix_free::{
    EnergyOperator2D, LaplacianOperator2D, LinearOperator, MomentumOperator1D, MomentumOperator2D,
    PoissonOperator3D,
};

pub use preconditioners::multigrid::AMGConfig;
pub use block_preconditioner::{
    BlockDiagonalPreconditioner, DiagonalPreconditioner, SimplePreconditioner,
};
pub use preconditioners::{
    AlgebraicMultigrid, DeflationPreconditioner, IdentityPreconditioner, IncompleteLU,
    JacobiPreconditioner, SORPreconditioner, SerialSchwarzPreconditioner, SSOR,
};

pub use traits::{Configurable, IterativeLinearSolver, LinearSolver, Preconditioner};

#[cfg(test)]
mod tests;
