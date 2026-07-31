//! CFD-specific linear solver extensions.
//!
//! The iterative solver algorithms (CG, BiCGSTAB, GMRES) are canonical in
//! `leto-ops` and re-exported via `cfd_math::iterative`. This module retains
//! CFD-domain-specific code: multigrid, ILU, block preconditioners, direct
//! solver bridge, and the tiered solver chain for saddle-point systems.

pub mod block_preconditioner;
pub mod chain;
pub mod dense_bridge;
pub mod direct_solver;
pub mod krylov;
pub mod preconditioners;

// Re-export the SSOT `Preconditioner` trait so internal sub-modules can use
// `crate::linear_solver::Preconditioner` without taking a direct `leto_ops` dependency.
pub use crate::iterative::Preconditioner;

// Re-export `IterativeSolverConfig` for legacy paths inside this module tree.
pub use crate::iterative::IterativeSolverConfig;

// Re-export solver types for convenience
pub use crate::iterative::{BiCGSTAB, GMRES};
pub use block_preconditioner::{
    BlockDiagonalPreconditioner, ComponentBlockPreconditioner, DiagonalPreconditioner,
    SimplePreconditioner,
};
pub use chain::{LinearSolverChain, LinearSolverState};
pub use direct_solver::DirectSparseSolver;
pub use preconditioners::multigrid::AMGConfig;
pub use preconditioners::{AlgebraicMultigrid, IncompleteLU};
