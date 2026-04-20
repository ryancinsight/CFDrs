//! Level Set Method for interface tracking in 3D multiphase flows
//!
//! The Level Set method represents interfaces as the zero level set of a
//! signed distance function, providing accurate interface tracking.
//!
//! # Theorem — Signed-Distance Invariant Under Reinitialization
//!
//! The reinitialization PDE converges to a steady state satisfying
//! $|\nabla \phi| = 1$ while preserving the zero-level set location, so interface
//! geometry is maintained while restoring numerical conditioning.

pub mod config;
mod advection;
pub mod solver;
mod weno;

// Re-export main types
pub use config::LevelSetConfig;
pub use solver::LevelSetSolver;
