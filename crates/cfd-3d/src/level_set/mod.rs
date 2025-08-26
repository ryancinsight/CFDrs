//! Level Set Method for interface tracking in 3D multiphase flows
//!
//! The Level Set method represents interfaces as the zero level set of a
//! signed distance function, providing accurate interface tracking.

pub mod config;
pub mod solver;
pub mod advection;
pub mod reinitialization;
pub mod narrow_band;
pub mod derivatives;

// Re-export main types
pub use config::LevelSetConfig;
pub use solver::LevelSetSolver;
pub use advection::AdvectionScheme;
pub use reinitialization::ReinitializationMethod;
pub use narrow_band::NarrowBandTracker;