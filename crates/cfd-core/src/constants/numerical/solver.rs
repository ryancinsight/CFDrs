//! Solver convergence and iteration constants
//!
//! Reference: Ferziger, J.H. & Perić, M. (2002). Computational Methods for Fluid Dynamics

/// Convergence tolerance for iterative solvers
pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;

/// Machine epsilon tolerance for numerical comparisons
pub const EPSILON_TOLERANCE: f64 = 1e-10;

/// Maximum outer iterations for pressure-velocity coupling
pub const MAX_ITERATIONS_OUTER: usize = 100;

/// Maximum inner iterations for linear solvers
pub const MAX_ITERATIONS_INNER: usize = 1000;

/// Logging interval for convergence monitoring
pub const LOG_INTERVAL: usize = 10;
