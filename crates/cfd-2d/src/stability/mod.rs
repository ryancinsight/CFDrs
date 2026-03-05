//! Numerical stability analysis and monitoring
//!
//! # Theorem (Von Neumann Stability — Necessary Condition)
//!
//! A finite difference scheme is stable if and only if all eigenvalues of its
//! amplification matrix $G$ satisfy $|\lambda_j| \le 1 + O(\Delta t)$.
//! The CFL, diffusion, and combined conditions in this module are derived
//! from the von Neumann analysis applied to the respective operators.

pub mod cfl;

pub use cfl::CFLCalculator;
