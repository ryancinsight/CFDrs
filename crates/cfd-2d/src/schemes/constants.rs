//! Constants for numerical schemes

/// Minimum threshold for flux limiters to avoid division by zero
pub const MIN_LIMITER_THRESHOLD: f64 = 1e-10;

/// WENO scheme epsilon for smoothness indicators
pub const WENO_EPSILON: f64 = 1e-6;

/// WENO scheme power parameter
pub const WENO_POWER: i32 = 2;

/// WENO5 optimal weights for uniform grid
pub const WENO5_WEIGHTS: [f64; 3] = [0.1, 0.6, 0.3];

/// MUSCL scheme kappa parameter (1/3 for third-order)
pub const MUSCL_KAPPA: f64 = 1.0 / 3.0;

/// CFL number for stability
pub const DEFAULT_CFL: f64 = 0.5;

/// Courant number limit for explicit schemes
pub const COURANT_LIMIT_EXPLICIT: f64 = 1.0;

/// Courant number limit for implicit schemes
pub const COURANT_LIMIT_IMPLICIT: f64 = 10.0;

// Numerical coefficients for schemes
/// Central difference divisor
pub const CENTRAL_DIFF_DIVISOR: f64 = 2.0;

/// Fourth-order central coefficient
pub const FOURTH_ORDER_COEFF_8: f64 = 8.0;

/// Fourth-order central divisor
pub const FOURTH_ORDER_DIVISOR: f64 = 12.0;

/// TVD coefficient two
pub const TVD_TWO: f64 = 2.0;
