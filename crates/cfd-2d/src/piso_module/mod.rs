//! Pressure-Implicit with Splitting of Operators (PISO) algorithm
//!
//! Reference: Issa, R.I. (1986) "Solution of the implicitly discretised fluid flow
//! equations by operator-splitting"

pub mod config;
pub mod constants;
pub mod predictor;
pub mod corrector;
pub mod solver;

pub use config::PisoConfig;
pub use constants::PisoConstants;
pub use solver::PisoSolver;

// Re-export the main solver for backward compatibility
pub use solver::PisoSolver as Solver;