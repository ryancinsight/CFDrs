//! Finite Volume Method (FVM) solver for 2D CFD simulations

pub mod config;
pub mod flux;
pub mod geometry;
pub mod solver;

// Re-export main types
pub use config::FvmConfig;
pub use flux::{FluxScheme, FluxSchemeFactory};
pub use geometry::Face;
pub use solver::FvmSolver;
