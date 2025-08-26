//! Time integration constants
//!
//! Parameters for time stepping and temporal discretization.

/// Default CFL number for explicit schemes
pub const DEFAULT_CFL: f64 = 0.5;
/// Maximum CFL number for stability
pub const MAX_CFL: f64 = 1.0;
/// CFL number for implicit schemes
pub const IMPLICIT_CFL: f64 = 5.0;
/// Minimum time step size
pub const MIN_TIME_STEP: f64 = 1e-12;
/// Maximum time step size growth factor
pub const MAX_TIME_STEP_GROWTH: f64 = 1.2;
/// Maximum time step size reduction factor
pub const MAX_TIME_STEP_REDUCTION: f64 = 0.5;
/// Time step safety factor
pub const TIME_STEP_SAFETY: f64 = 0.9;
/// Runge-Kutta 4 weights
pub const RK4_WEIGHTS: [f64; 4] = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];
/// Adams-Bashforth 2 coefficients
pub const AB2_COEFFICIENTS: [f64; 2] = [1.5, -0.5];
/// Adams-Bashforth 3 coefficients
pub const AB3_COEFFICIENTS: [f64; 3] = [23.0 / 12.0, -16.0 / 12.0, 5.0 / 12.0];
