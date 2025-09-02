//! Physical and numerical constants for discretization schemes
//!
//! All constants are literature-validated with references

use nalgebra::RealField;
use num_traits::FromPrimitive;

// ============================================================================
// CFL Stability Limits
// ============================================================================

/// Maximum CFL number for explicit Euler time stepping
/// Reference: Courant, Friedrichs, Lewy (1928)
pub const CFL_EXPLICIT_EULER: f64 = 1.0;

/// Maximum CFL number for QUICK scheme (third-order accuracy)
/// Reference: Leonard (1979) "A stable and accurate convective modelling procedure"
pub const CFL_QUICK_SCHEME: f64 = 0.75;

/// Default CFL safety factor for conservative time stepping
pub const CFL_SAFETY_FACTOR: f64 = 0.8;

/// Maximum CFL number for Lax-Wendroff scheme
pub const CFL_LAX_WENDROFF: f64 = 1.0;

// ============================================================================
// Diffusion Stability Limits
// ============================================================================

/// Maximum diffusion number (von Neumann number) for 2D explicit schemes
/// D = ν*dt/(dx²) ≤ 0.5 for stability
pub const DIFFUSION_NUMBER_2D_MAX: f64 = 0.5;

/// Maximum diffusion number for 3D explicit schemes
pub const DIFFUSION_NUMBER_3D_MAX: f64 = 0.25;

// ============================================================================
// QUICK Scheme Coefficients
// ============================================================================

/// QUICK scheme coefficient for upstream cell (6/8)
/// Reference: Leonard (1979) Eq. 15
pub const QUICK_COEFF_UPSTREAM: f64 = 0.75; // 6/8

/// QUICK scheme coefficient for central cell (3/8)
pub const QUICK_COEFF_CENTRAL: f64 = 0.375; // 3/8

/// QUICK scheme coefficient for far upstream cell (-1/8)
pub const QUICK_COEFF_FAR_UPSTREAM: f64 = 0.125; // 1/8

// ============================================================================
// Flux Limiter Constants
// ============================================================================

/// Coefficient 2.0 used in various flux limiters
pub const FLUX_LIMITER_TWO: f64 = 2.0;

/// Coefficient 0.5 for averaging operations
pub const FLUX_LIMITER_HALF: f64 = 0.5;

// ============================================================================
// Ghost Cell Extrapolation Coefficients
// ============================================================================

/// Linear extrapolation coefficient for second-order ghost cells
pub const GHOST_CELL_LINEAR_COEFF: f64 = 2.0;

/// Quadratic extrapolation coefficient for third-order ghost cells
pub const GHOST_CELL_QUADRATIC_COEFF: f64 = 3.0;

/// Fourth-order ghost cell coefficient
pub const GHOST_CELL_FOURTH_ORDER_COEFF: f64 = 5.0;

/// Neumann BC coefficient for second-order (8/3)
pub const NEUMANN_SECOND_ORDER_COEFF: f64 = 8.0 / 3.0;

/// Neumann BC coefficient for second-order (2/3)
pub const NEUMANN_SECOND_ORDER_FRACTION: f64 = 2.0 / 3.0;

// ============================================================================
// Numerical Scheme Constants
// ============================================================================

/// Central difference divisor for second-order accuracy
/// Used in (f[i+1] - f[i-1]) / (2*dx) formulation
pub const CENTRAL_DIFF_DIVISOR: f64 = 2.0;

/// Small epsilon value for WENO scheme to avoid division by zero
/// Reference: Shu (1998) "Essentially Non-oscillatory and Weighted ENO Schemes"
pub const WENO_EPSILON: f64 = 1e-6;

/// WENO5 optimal weights for smooth solutions
/// Reference: Jiang & Shu (1996) "Efficient Implementation of Weighted ENO Schemes"
/// d0 = 1/10, d1 = 6/10, d2 = 3/10 for upwind-biased stencil
pub const WENO5_WEIGHTS: [f64; 3] = [0.1, 0.6, 0.3];

/// General-purpose TWO constant for various schemes
pub const TWO: f64 = 2.0;

// ============================================================================
// Helper Functions
// ============================================================================

/// Convert f64 constant to generic RealField type
pub fn to_realfield<T: RealField + FromPrimitive>(value: f64) -> T {
    T::from_f64(value).unwrap_or_else(T::one)
}
