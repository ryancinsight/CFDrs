//! Named constants for 1D CFD and microfluidics
//!
//! This module provides a single source of truth for all numerical constants
//! used in 1D network simulations, following the SSOT principle.

use nalgebra::RealField;

// Numerical tolerances
/// Default convergence tolerance
pub const DEFAULT_TOLERANCE: f64 = 1e-6;
/// Machine epsilon for comparisons
pub const EPSILON: f64 = 1e-10;
/// Small number to avoid division by zero
pub const SMALL_NUMBER: f64 = 1e-14;

// Resistance calculation constants
/// Poiseuille constant for circular channels (64)
pub const POISEUILLE_CIRCULAR: f64 = 64.0;
/// Poiseuille constant for rectangular channels (wide limit)
pub const POISEUILLE_RECTANGULAR_BASE: f64 = 96.0;
/// Aspect ratio correction factor
pub const ASPECT_RATIO_CORRECTION: f64 = 0.63;
/// Minimum aspect ratio for numerical stability
pub const MIN_ASPECT_RATIO: f64 = 0.1;
/// Maximum aspect ratio for accurate approximation
pub const MAX_ASPECT_RATIO: f64 = 10.0;

// Entrance length correlations
/// Shah correlation coefficient for laminar flow
pub const SHAH_COEFFICIENT: f64 = 0.06;
/// Entrance effect exponent
pub const ENTRANCE_EXPONENT: f64 = 1.0;

// Network analysis constants
/// Default maximum iterations for network solver
pub const DEFAULT_MAX_ITERATIONS: usize = 1000;
/// Default relaxation factor for iterative solvers
pub const DEFAULT_RELAXATION: f64 = 0.8;

// Physical properties (SI units)
/// Standard atmospheric pressure [Pa]
pub const STANDARD_PRESSURE: f64 = 101325.0;
/// Standard temperature [K]
pub const STANDARD_TEMPERATURE: f64 = 293.15; // 20°C
/// Water density at standard conditions [kg/m³]
pub const WATER_DENSITY: f64 = 998.2;
/// Water dynamic viscosity at 20°C [Pa·s]
pub const WATER_VISCOSITY: f64 = 1.002e-3;
/// Air density at standard conditions [kg/m³]
pub const AIR_DENSITY: f64 = 1.204;
/// Air dynamic viscosity at 20°C [Pa·s]
pub const AIR_VISCOSITY: f64 = 1.81e-5;

// Microfluidics specific
/// Mean free path of air at STP [m]
pub const AIR_MEAN_FREE_PATH: f64 = 68e-9;
/// Slip flow regime threshold (Kn < 0.1)
pub const KNUDSEN_SLIP_THRESHOLD: f64 = 0.1;
/// Continuum flow threshold (Kn < 0.01)
pub const KNUDSEN_CONTINUUM_THRESHOLD: f64 = 0.01;

// Non-Newtonian fluid parameters
/// Typical yield stress for Bingham plastics [Pa]
pub const TYPICAL_YIELD_STRESS: f64 = 1.0;
/// Power-law index for shear-thinning fluids
pub const SHEAR_THINNING_INDEX: f64 = 0.8;
/// Power-law index for shear-thickening fluids
pub const SHEAR_THICKENING_INDEX: f64 = 1.2;

/// Get default tolerance for type T
pub fn default_tolerance<T: RealField>() -> T {
    T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::default_epsilon() * T::from_f64(100.0).unwrap())
}

/// Get epsilon for type T
pub fn epsilon<T: RealField>() -> T {
    T::from_f64(EPSILON).unwrap_or_else(T::default_epsilon)
}

/// Get Poiseuille constant for circular channels
pub fn poiseuille_circular<T: RealField>() -> T {
    T::from_f64(POISEUILLE_CIRCULAR).unwrap()
}

/// Get Poiseuille constant for rectangular channels (base)
pub fn poiseuille_rectangular_base<T: RealField>() -> T {
    T::from_f64(POISEUILLE_RECTANGULAR_BASE).unwrap()
}

/// Get aspect ratio correction factor
pub fn aspect_ratio_correction<T: RealField>() -> T {
    T::from_f64(ASPECT_RATIO_CORRECTION).unwrap()
}

/// Get minimum safe aspect ratio
pub fn min_aspect_ratio<T: RealField>() -> T {
    T::from_f64(MIN_ASPECT_RATIO).unwrap()
}

/// Get maximum accurate aspect ratio
pub fn max_aspect_ratio<T: RealField>() -> T {
    T::from_f64(MAX_ASPECT_RATIO).unwrap()
}