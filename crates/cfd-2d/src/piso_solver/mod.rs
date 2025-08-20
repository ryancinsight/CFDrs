//! PISO (Pressure-Implicit with Splitting of Operators) algorithm
//!
//! Reference: Issa, R.I. (1986) "Solution of the implicitly discretised fluid flow
//! equations by operator-splitting", Journal of Computational Physics, 62(1), 40-65

pub mod config;
pub mod constants;
pub mod predictor;
pub mod corrector;
pub mod convergence;
pub mod solver;

pub use config::PisoConfig;
pub use constants::PisoConstants;
pub use solver::PisoSolver;

// Re-export for backward compatibility (will be removed)
pub use solver::PisoSolver as Solver;