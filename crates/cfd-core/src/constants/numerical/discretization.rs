//! Discretization scheme constants
//!
//! Parameters for various numerical discretization schemes.

/// Central difference scheme weight
pub const CENTRAL_DIFFERENCE_WEIGHT: f64 = 0.5;

/// Upwind scheme blending factor
pub const UPWIND_BLEND_FACTOR: f64 = 1.0;

/// QUICK scheme weight
pub const QUICK_WEIGHT: f64 = 0.875;

/// TVD limiter threshold
pub const TVD_LIMITER_THRESHOLD: f64 = 2.0;

/// Minimum limiter value to prevent division by zero
pub const LIMITER_EPSILON: f64 = 1e-10;

/// Van Leer limiter smoothness parameter
pub const VAN_LEER_KAPPA: f64 = 1.0 / 3.0;

/// MUSCL scheme parameter
pub const MUSCL_KAPPA: f64 = 1.0 / 3.0;

/// Flux limiter small number
pub const FLUX_LIMITER_EPSILON: f64 = 1e-20;

/// Grid quality thresholds
pub const MIN_ASPECT_RATIO: f64 = 0.1;
pub const MAX_ASPECT_RATIO: f64 = 10.0;
pub const MAX_SKEWNESS: f64 = 0.95;
pub const ORTHOGONALITY_THRESHOLD: f64 = 0.1;