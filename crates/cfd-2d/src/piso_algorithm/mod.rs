//! Pressure-Implicit with Splitting of Operators (PISO) algorithm
//!
//! Reference: Issa, R.I. (1986) "Solution of the implicitly discretised fluid flow
//! equations by operator-splitting", Journal of Computational Physics, 62(1), 40-65

pub mod config;
pub mod convergence;
pub mod corrector;
pub mod predictor;
pub mod solver;
pub use config::PisoConfig;
pub use convergence::ConvergenceCriteria;
pub use solver::PisoSolver;
