//! Named constants for 2D CFD algorithms
//!
//! This module provides a single source of truth for all numerical constants
//! used in 2D CFD simulations, following the SSOT principle.

use nalgebra::RealField;

// Numerical tolerances
/// Default convergence tolerance for iterative solvers
pub const DEFAULT_TOLERANCE: f64 = 1e-6;
/// Machine epsilon for floating point comparisons
pub const EPSILON: f64 = 1e-10;
/// Small number to avoid division by zero
pub const SMALL_NUMBER: f64 = 1e-14;

// Solver parameters
/// Default maximum iterations for iterative solvers
pub const DEFAULT_MAX_ITERATIONS: usize = 1000;
/// Default CFL number for stability
pub const DEFAULT_CFL_NUMBER: f64 = 0.5;
/// Default relaxation factor for SIMPLE algorithm
pub const DEFAULT_RELAXATION_FACTOR: f64 = 0.7;
/// Optimal SOR relaxation factor for Poisson equation
pub const SOR_OPTIMAL_FACTOR: f64 = 1.85;

// PISO algorithm parameters
/// Default number of correctors for PISO
pub const DEFAULT_MAX_CORRECTORS: usize = 2;
/// Default number of non-orthogonal correctors
pub const DEFAULT_NON_ORTHOGONAL_CORRECTORS: usize = 1;

// LBM parameters
/// Speed of sound squared in lattice units (c_s^2 = 1/3)
pub const LBM_CS_SQUARED: f64 = 1.0 / 3.0;
/// Lattice speed
pub const LBM_LATTICE_SPEED: f64 = 1.0;

// Numerical scheme factors
/// Factor for gradient calculations (central difference)
pub const GRADIENT_FACTOR: f64 = 2.0;
/// Factor for Laplacian stencil center coefficient
pub const LAPLACIAN_CENTER_FACTOR: f64 = -4.0;
/// Factor for second-order accuracy
pub const SECOND_ORDER_FACTOR: f64 = 2.0;

// Physical constants
/// Reference density for water at 20°C [kg/m³]
pub const WATER_DENSITY: f64 = 998.2;
/// Reference dynamic viscosity for water at 20°C [Pa·s]
pub const WATER_VISCOSITY: f64 = 1.002e-3;
/// Reference density for air at 20°C, 1 atm [kg/m³]
pub const AIR_DENSITY: f64 = 1.204;
/// Reference dynamic viscosity for air at 20°C [Pa·s]
pub const AIR_VISCOSITY: f64 = 1.81e-5;

/// Get default tolerance for type T
pub fn default_tolerance<T: RealField>() -> T {
    T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::default_epsilon() * T::from_f64(100.0).unwrap())
}

/// Get epsilon for type T
pub fn epsilon<T: RealField>() -> T {
    T::from_f64(EPSILON).unwrap_or_else(T::default_epsilon)
}

/// Get small number for type T
pub fn small_number<T: RealField>() -> T {
    T::from_f64(SMALL_NUMBER).unwrap_or_else(T::default_epsilon)
}

/// Get gradient factor for type T
pub fn gradient_factor<T: RealField>() -> T {
    T::from_f64(GRADIENT_FACTOR).unwrap()
}

/// Get Laplacian center factor for type T
pub fn laplacian_center_factor<T: RealField>() -> T {
    T::from_f64(LAPLACIAN_CENTER_FACTOR).unwrap()
}

/// Get LBM speed of sound squared for type T
pub fn lbm_cs_squared<T: RealField>() -> T {
    T::from_f64(LBM_CS_SQUARED).unwrap()
}