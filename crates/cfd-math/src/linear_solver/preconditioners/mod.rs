//! CFD-specific preconditioners.
//!
//! Simple preconditioners (Identity, Jacobi, SOR, SSOR) are SSOT in
//! `leto-ops` and re-exported via `cfd_math::iterative::preconditioners`.
//!
//! CFD-domain-specific implementations (AMG, ILU) live here.
//! IncompleteCholesky, Deflation, and Schwarz are available in the git
//! history pending full migration to the leto-ops Preconditioner trait.

pub mod ilu;
pub mod multigrid;

pub use ilu::IncompleteLU;
pub use multigrid::AlgebraicMultigrid;
