//! Numerical method constants

/// Default convergence tolerance for iterative solvers (deprecated)
#[deprecated(note = "Use CONVERGENCE_TOLERANCE instead")]
pub const DEFAULT_CONVERGENCE_TOLERANCE: f64 = 1e-6;

/// Maximum iterations for iterative solvers (deprecated)
#[deprecated(note = "Use MAX_ITERATIONS_DEFAULT instead")]
pub const DEFAULT_MAX_ITERATIONS: usize = 1000;

/// Relaxation factor for pressure correction
pub const DEFAULT_PRESSURE_RELAXATION: f64 = 0.7;

/// Relaxation factor for velocity correction
pub const DEFAULT_VELOCITY_RELAXATION: f64 = 0.3;

/// CFL number for stability
pub const DEFAULT_CFL_NUMBER: f64 = 0.5;

/// Small number to prevent division by zero
pub const EPSILON: f64 = 1e-14;