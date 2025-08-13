//! Centralized constants module for the entire CFD simulation suite
//!
//! This module provides a single source of truth (SSOT) for all numerical and physical
//! constants used throughout the CFD codebase, eliminating duplication and ensuring
//! consistency across all modules.

use nalgebra::RealField;

// ============================================================================
// NUMERICAL CONSTANTS
// ============================================================================

/// Machine epsilon for floating point comparisons
pub const EPSILON: f64 = 1e-10;

/// Small number to prevent division by zero
pub const SMALL_NUMBER: f64 = 1e-14;

/// Default convergence tolerance for iterative solvers
pub const DEFAULT_TOLERANCE: f64 = 1e-6;

/// Tight tolerance for high-precision calculations
pub const TIGHT_TOLERANCE: f64 = 1e-9;

/// Loose tolerance for preliminary calculations
pub const LOOSE_TOLERANCE: f64 = 1e-3;

// ============================================================================
// MATHEMATICAL CONSTANTS
// ============================================================================

/// Common numerical factors
/// Zero value
pub const ZERO: f64 = 0.0;
/// Unity value
pub const ONE: f64 = 1.0;
/// Two value
pub const TWO: f64 = 2.0;
/// Three value  
pub const THREE: f64 = 3.0;
/// Four value
pub const FOUR: f64 = 4.0;
/// Five value
pub const FIVE: f64 = 5.0;
/// Six value
pub const SIX: f64 = 6.0;
/// Half value
pub const HALF: f64 = 0.5;
/// One third value
pub const ONE_THIRD: f64 = 1.0 / 3.0;
/// Two thirds value
pub const TWO_THIRDS: f64 = 2.0 / 3.0;
/// One fourth value
pub const ONE_FOURTH: f64 = 0.25;
/// One sixth value
pub const ONE_SIXTH: f64 = 1.0 / 6.0;

// ============================================================================
// SOLVER PARAMETERS
// ============================================================================

/// Default maximum iterations for iterative solvers
pub const DEFAULT_MAX_ITERATIONS: usize = 1000;

/// Default CFL number for stability
pub const DEFAULT_CFL_NUMBER: f64 = 0.5;

/// Maximum safe CFL number
pub const MAX_CFL_NUMBER: f64 = 1.0;

/// Default time step factor for SIMPLE solver
/// This is multiplied with CFL number to get initial time step
pub const DEFAULT_TIME_STEP_FACTOR: f64 = 0.02;

/// Default time step for transient simulations (seconds)
pub const DEFAULT_TIME_STEP: f64 = 0.01;

/// Default relaxation factor for iterative methods
pub const DEFAULT_RELAXATION_FACTOR: f64 = 0.7;

/// Optimal SOR relaxation factor for Poisson equation on uniform grid
pub const SOR_OPTIMAL_FACTOR: f64 = 1.85;

/// Default under-relaxation for pressure in SIMPLE algorithm
pub const PRESSURE_UNDER_RELAXATION: f64 = 0.3;

/// Default under-relaxation for velocity in SIMPLE algorithm  
pub const VELOCITY_UNDER_RELAXATION: f64 = 0.7;

// ============================================================================
// PHYSICAL PROPERTIES (SI Units)
// ============================================================================

/// Standard atmospheric pressure [Pa]
pub const STANDARD_PRESSURE: f64 = 101325.0;

/// Standard temperature [K] (20°C)
pub const STANDARD_TEMPERATURE: f64 = 293.15;

/// Water density at 20°C [kg/m³]
pub const WATER_DENSITY: f64 = 998.2;

/// Water dynamic viscosity at 20°C [Pa·s]
pub const WATER_VISCOSITY: f64 = 1.002e-3;

/// Water kinematic viscosity at 20°C [m²/s] (derived from dynamic viscosity and density)
pub const WATER_KINEMATIC_VISCOSITY: f64 = WATER_VISCOSITY / WATER_DENSITY;

/// Air density at 20°C, 1 atm [kg/m³]
pub const AIR_DENSITY: f64 = 1.204;

/// Air dynamic viscosity at 20°C [Pa·s]
pub const AIR_VISCOSITY: f64 = 1.81e-5;

/// Air kinematic viscosity at 20°C [m²/s]
pub const AIR_KINEMATIC_VISCOSITY: f64 = 1.50e-5;

/// Gravitational acceleration [m/s²]
pub const GRAVITY: f64 = 9.81;

// ============================================================================
// FLUID DYNAMICS CONSTANTS
// ============================================================================

/// Reynolds number thresholds
pub const LAMINAR_THRESHOLD: f64 = 2300.0;
pub const TURBULENT_THRESHOLD: f64 = 4000.0;

/// Prandtl number for air at 20°C
pub const AIR_PRANDTL: f64 = 0.71;

/// Prandtl number for water at 20°C
pub const WATER_PRANDTL: f64 = 7.0;

/// Von Karman constant for wall functions
pub const VON_KARMAN: f64 = 0.41;

/// Wall function constants
pub const WALL_FUNCTION_E: f64 = 9.8;
pub const Y_PLUS_LAMINAR: f64 = 11.63;
pub const Y_PLUS_BUFFER_START: f64 = 5.0;
pub const Y_PLUS_BUFFER_END: f64 = 30.0;

// ============================================================================
// CHANNEL FLOW CONSTANTS
// ============================================================================

/// Poiseuille number for circular channels
pub const POISEUILLE_CIRCULAR: f64 = 64.0;

/// Poiseuille number for rectangular channels (base value for square)
pub const POISEUILLE_RECTANGULAR_BASE: f64 = 96.0;

/// Aspect ratio correction factor for rectangular channels
pub const ASPECT_RATIO_CORRECTION: f64 = 0.63;

/// Minimum aspect ratio for numerical stability
pub const MIN_ASPECT_RATIO: f64 = 0.1;

/// Maximum aspect ratio for mesh element validity
pub const MAX_ASPECT_RATIO: f64 = 10.0;

/// Shah correlation coefficient for entrance length
pub const SHAH_COEFFICIENT: f64 = 0.06;

/// Entrance effect exponent
pub const ENTRANCE_EXPONENT: f64 = 1.0;

// ============================================================================
// LATTICE BOLTZMANN METHOD CONSTANTS
// ============================================================================

/// Speed of sound squared in lattice units (cs² = 1/3)
pub const LBM_CS_SQUARED: f64 = ONE_THIRD;

/// Lattice speed in lattice units
pub const LBM_LATTICE_SPEED: f64 = 1.0;

/// D2Q9 lattice weights
pub const D2Q9_WEIGHT_CENTER: f64 = 4.0 / 9.0;
pub const D2Q9_WEIGHT_FACE: f64 = 1.0 / 9.0;
pub const D2Q9_WEIGHT_CORNER: f64 = 1.0 / 36.0;

// ============================================================================
// TURBULENCE MODEL CONSTANTS
// ============================================================================

/// Standard Smagorinsky constant
pub const SMAGORINSKY_CONSTANT: f64 = 0.17;

/// k-epsilon model constants
pub const K_EPSILON_CMU: f64 = 0.09;
pub const K_EPSILON_C1: f64 = 1.44;
pub const K_EPSILON_C2: f64 = 1.92;
pub const K_EPSILON_SIGMA_K: f64 = 1.0;
pub const K_EPSILON_SIGMA_E: f64 = 1.3;

/// Test filter ratio for Germano-Lilly model
pub const TEST_FILTER_RATIO: f64 = 2.0;

// ============================================================================
// NUMERICAL SCHEME CONSTANTS
// ============================================================================

/// QUICK scheme coefficients
pub const QUICK_COEFF_UPSTREAM: f64 = 3.0 / 8.0;
pub const QUICK_COEFF_DOWNSTREAM: f64 = 6.0 / 8.0;
pub const QUICK_COEFF_FAR_UPSTREAM: f64 = -1.0 / 8.0;

/// WENO5 scheme order
pub const WENO_ORDER: usize = 5;
pub const WENO_EPSILON: f64 = 1e-6;

/// Gradient calculation factor (for central differences)
pub const GRADIENT_FACTOR: f64 = 2.0;

/// Laplacian stencil center coefficient (2D)
pub const LAPLACIAN_CENTER_2D: f64 = -4.0;

/// Laplacian stencil center coefficient (3D)
pub const LAPLACIAN_CENTER_3D: f64 = -6.0;

// ============================================================================
// MICROFLUIDICS CONSTANTS
// ============================================================================

/// Mean free path of air at STP [m]
pub const AIR_MEAN_FREE_PATH: f64 = 68e-9;

/// Knudsen number thresholds
pub const KNUDSEN_CONTINUUM_THRESHOLD: f64 = 0.01;
pub const KNUDSEN_SLIP_THRESHOLD: f64 = 0.1;
pub const KNUDSEN_TRANSITION_THRESHOLD: f64 = 10.0;

// ============================================================================
// NON-NEWTONIAN FLUID PARAMETERS
// ============================================================================

/// Typical yield stress for Bingham plastics [Pa]
pub const TYPICAL_YIELD_STRESS: f64 = 1.0;

/// Power-law indices
pub const SHEAR_THINNING_INDEX: f64 = 0.8;
pub const NEWTONIAN_INDEX: f64 = 1.0;
pub const SHEAR_THICKENING_INDEX: f64 = 1.2;

/// Viscosity limits for numerical stability
pub const MIN_VISCOSITY: f64 = 1e-10;
pub const MAX_VISCOSITY: f64 = 1e10;
pub const SOLID_LIKE_VISCOSITY: f64 = 1e6;
pub const YIELD_STRESS_VISCOSITY: f64 = 1e10;

// ============================================================================
// MULTIPHASE FLOW CONSTANTS
// ============================================================================

/// Interface thickness for diffuse interface methods
pub const INTERFACE_THICKNESS: f64 = 1.5;

/// VOF compression factor
pub const VOF_COMPRESSION_FACTOR: f64 = 1.0;

/// Level set reinitialization parameters
pub const LEVEL_SET_REINIT_ITERATIONS: usize = 10;
pub const LEVEL_SET_BAND_WIDTH: f64 = 5.0;

/// Surface tension coefficient for water-air at 20°C [N/m]
pub const WATER_AIR_SURFACE_TENSION: f64 = 0.0728;

// ============================================================================
// MESH QUALITY METRICS
// ============================================================================

/// Minimum acceptable mesh quality
pub const MIN_MESH_QUALITY: f64 = 0.1;

/// Maximum acceptable aspect ratio
pub const MAX_MESH_ASPECT_RATIO: f64 = 100.0;

/// Maximum acceptable skewness
pub const MAX_MESH_SKEWNESS: f64 = 0.95;

/// Epsilon multiplier for fallback tolerance calculation
/// Used when DEFAULT_TOLERANCE cannot be converted to type T
pub const EPSILON_MULTIPLIER_FALLBACK: f64 = 100.0;

// ============================================================================
// HELPER FUNCTIONS FOR GENERIC TYPES
// ============================================================================

/// Get default tolerance for type T
pub fn default_tolerance<T: RealField>() -> T {
    T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::default_epsilon() * T::from_f64(EPSILON_MULTIPLIER_FALLBACK).unwrap())
}

/// Get epsilon for type T
pub fn epsilon<T: RealField>() -> T {
    T::from_f64(EPSILON).unwrap_or_else(T::default_epsilon)
}

/// Get small number for type T
pub fn small_number<T: RealField>() -> T {
    T::from_f64(SMALL_NUMBER).unwrap_or_else(T::default_epsilon)
}

/// Get a numerical constant for type T
pub fn get_constant<T: RealField>(value: f64) -> T {
    T::from_f64(value).unwrap_or_else(T::zero)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_constants_consistency() {
        // Verify derived constants
        assert_eq!(WATER_KINEMATIC_VISCOSITY, WATER_VISCOSITY / WATER_DENSITY);
        assert!((AIR_KINEMATIC_VISCOSITY - AIR_VISCOSITY / AIR_DENSITY).abs() < 1e-7);
        
        // Verify mathematical constants
        assert_eq!(HALF, 1.0 / 2.0);
        assert_eq!(ONE_THIRD, 1.0 / 3.0);
        assert_eq!(ONE_SIXTH, 1.0 / 6.0);
        
        // Verify LBM weights sum to 1
        let d2q9_sum = D2Q9_WEIGHT_CENTER + 4.0 * D2Q9_WEIGHT_FACE + 4.0 * D2Q9_WEIGHT_CORNER;
        assert!((d2q9_sum - 1.0).abs() < EPSILON);
    }
    
    #[test]
    fn test_generic_functions() {
        let tol_f64: f64 = default_tolerance();
        assert_eq!(tol_f64, DEFAULT_TOLERANCE);
        
        let eps_f64: f64 = epsilon();
        assert_eq!(eps_f64, EPSILON);
        
        let small_f64: f64 = small_number();
        assert_eq!(small_f64, SMALL_NUMBER);
        
        let const_f64: f64 = get_constant(TWO);
        assert_eq!(const_f64, 2.0);
    }
}