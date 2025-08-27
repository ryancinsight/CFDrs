//! Level Set Method for interface tracking in 3D multiphase flows
//!
//! The Level Set method represents interfaces as the zero level set of a
//! signed distance function, providing accurate interface tracking.

pub mod config;
pub mod solver;

// Re-export main types
pub use config::LevelSetConfig;
pub use solver::LevelSetSolver;
