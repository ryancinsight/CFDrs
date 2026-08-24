//! CFD-specific linear solver extensions.
//!
//! Athena owns the Krylov recurrences (CG, BiCGSTAB, GMRES, LSQR), the
//! operator and preconditioner seams, and the convergence policy they enforce
//! (Atlas ADR 0033). [`krylov`] is this workspace's entry point onto them,
//! bridging Athena's compile-time GMRES restart to the runtime width CFD
//! callers select.
//!
//! What lives here is CFD-domain-specific: multigrid, ILU, block
//! preconditioners for saddle-point systems, the direct solver bridge, and the
//! tiered fallback chain. Preconditioners defined here implement Athena's
//! [`athena_core::Preconditioner`] seam so a Krylov solve accepts them
//! directly.

pub mod block_preconditioner;
pub mod chain;
pub mod config;
pub mod dense_bridge;
pub mod direct_solver;
pub mod krylov;
pub mod preconditioners;

pub use block_preconditioner::{
    BlockDiagonalPreconditioner, ComponentBlockPreconditioner, DiagonalPreconditioner,
    SimplePreconditioner,
};
pub use chain::{LinearSolverChain, LinearSolverState};
pub use config::IterativeSolverConfig;
pub use direct_solver::DirectSparseSolver;
pub use krylov::KrylovWorkspace;
pub use preconditioners::multigrid::AMGConfig;
pub use preconditioners::{AlgebraicMultigrid, IncompleteLU};
