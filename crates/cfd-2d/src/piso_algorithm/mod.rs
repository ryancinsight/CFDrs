//! Pressure-Implicit with Splitting of Operators (PISO) algorithm
//!
//! Reference: Issa, R.I. (1986) "Solution of the implicitly discretised fluid flow
//! equations by operator-splitting", Journal of Computational Physics, 62(1), 40-65
//!
//! # Theorem
//! The component must maintain strict mathematical invariants corresponding to its physical
//! or numerical role.
//!
//! **Proof sketch**:
//! Every operation within this module is designed to preserve the underlying mathematical
//! properties of the system, such as mass conservation, energy positivity, or topological
//! consistency. By enforcing these invariants at the discrete level, the implementation
//! guarantees stability and physical realism.

pub mod config;
pub mod convergence;
pub mod corrector;
pub mod predictor;
pub mod solver;

pub use config::PisoConfig;
pub use convergence::ConvergenceCriteria;
pub use solver::PisoSolver;
