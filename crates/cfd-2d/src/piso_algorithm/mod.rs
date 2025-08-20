//! Pressure-Implicit with Splitting of Operators (PISO) algorithm
//!
//! Reference: Issa, R.I. (1986) "Solution of the implicitly discretised fluid flow
//! equations by operator-splitting", Journal of Computational Physics, 62(1), 40-65

pub mod config;
pub mod predictor;
pub mod corrector;
pub mod convergence;
pub mod solver;

pub use config::PisoConfig;
pub use solver::PisoSolver;
pub use convergence::ConvergenceCriteria;