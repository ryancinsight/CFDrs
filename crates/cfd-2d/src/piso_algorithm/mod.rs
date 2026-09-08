//! Pressure-Implicit with Splitting of Operators (PISO) algorithm
//!
//! Reference: Issa, R.I. (1986) "Solution of the implicitly discretised fluid flow
//! equations by operator-splitting", Journal of Computational Physics, 62(1), 40-65
//!
//! # Theorem (PISO Splitting Error — Issa 1986)
//!
//! For a linearized incompressible problem with an accurate pressure-correction
//! solve, $m$ corrector steps reduce the operator-splitting error to
//! $O(\Delta t^{m+1})$. Two correctors therefore preserve second-order temporal
//! accuracy in the small-CFL regime.
//!
//! **Proof sketch**:
//! See the detailed proof in [`solver`]. The estimate assumes the predictor and
//! correctors are applied consistently to the same discretized momentum stencil.

pub mod config;
pub mod convergence;
pub mod corrector;
pub mod predictor;
pub mod solver;

pub use config::PisoConfig;
pub use convergence::ConvergenceCriteria;
pub use solver::PisoSolver;
