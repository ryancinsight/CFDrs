//! CFD-specific linear solver extensions.
//!
//! Athena owns the Krylov recurrences (CG, BiCGSTAB, GMRES, LSQR), the
//! operator and preconditioner seams, and the convergence policy they enforce
//! (Atlas ADR 0033). [`krylov`] is this workspace's entry point onto them: it
//! translates the CFD configuration vocabulary into a validated Athena policy
//! and wraps a caller-owned CSR matrix in Athena's operator seam, while the
//! one runtime restart bridge, [`athena_leto::KrylovWorkspace`], lives in
//! Athena (Atlas ADR 0062).
//!
//! What lives here is CFD-domain-specific: the CFD configuration
//! vocabulary, the translation into Athena's convergence policy, and the
//! interpretation of a solve outcome. The tiered fallback chain itself is
//! [`athena_leto::solver_chain`] (Atlas ADR 0062). Preconditioners defined
//! here implement Athena's
//! [`athena_core::Preconditioner`] seam so a Krylov solve accepts them
//! directly.

pub mod config;
pub mod krylov;
pub mod preconditioners;

pub use config::IterativeSolverConfig;
