//! Solver-related numerical constants
//!
//! Convergence tolerances, iteration limits, and solver parameters.

/// Default convergence tolerance for iterative solvers
pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;

/// Stricter convergence tolerance for high-accuracy simulations
pub const CONVERGENCE_TOLERANCE_STRICT: f64 = 1e-10;

/// Machine epsilon safety factor
pub const EPSILON_TOLERANCE: f64 = 1e-14;

/// Maximum iterations for outer loops
pub const MAX_ITERATIONS_OUTER: usize = 1000;

/// Maximum iterations for inner loops
pub const MAX_ITERATIONS_INNER: usize = 100;

/// Minimum iterations before convergence check
pub const MIN_ITERATIONS: usize = 3;

/// Default relaxation factor for iterative solvers
pub const DEFAULT_RELAXATION: f64 = 1.0;

/// Under-relaxation factor for pressure
pub const PRESSURE_RELAXATION: f64 = 0.3;

/// Under-relaxation factor for velocity
pub const VELOCITY_RELAXATION: f64 = 0.7;

/// Stagnation detection window size
pub const STAGNATION_WINDOW: usize = 10;

/// Stagnation tolerance ratio
pub const STAGNATION_TOLERANCE: f64 = 1e-3;