//! Centralized solver and numerical constants
//!
//! All numerical constants used in solvers, convergence checks, and numerical methods
//! are defined here to maintain SSOT and enable literature validation.

// ============================================================================
// CONVERGENCE TOLERANCES
// ============================================================================

/// Default iterative solver tolerance
/// Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems, p. 40
pub const DEFAULT_SOLVER_TOLERANCE: f64 = 1e-6;

/// Tight tolerance for pressure solvers in incompressible flow
/// Reference: Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow, p. 126
pub const PRESSURE_SOLVER_TOLERANCE: f64 = 1e-8;

/// Machine epsilon safety factor for zero comparisons
/// Reference: Goldberg, D. (1991). What Every Computer Scientist Should Know About Floating-Point
pub const EPSILON_TOLERANCE: f64 = 1e-10;

/// Ultra-tight tolerance for spectral methods
/// Reference: Boyd, J.P. (2001). Chebyshev and Fourier Spectral Methods, p. 98
pub const SPECTRAL_TOLERANCE: f64 = 1e-12;

// ============================================================================
// ITERATION LIMITS
// ============================================================================

/// Default maximum iterations for iterative solvers
/// Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems, p. 272
pub const DEFAULT_MAX_ITERATIONS: usize = 1000;

/// Maximum iterations for pressure correction in PISO/SIMPLE
/// Reference: Issa, R.I. (1986). Solution of implicitly discretised fluid flow equations
pub const PRESSURE_CORRECTION_MAX_ITERATIONS: usize = 100;

/// Maximum outer iterations for non-linear problems
/// Reference: Ferziger & Peric (2002). Computational Methods for Fluid Dynamics, p. 178
pub const NONLINEAR_MAX_ITERATIONS: usize = 50;

// ============================================================================
// BOUNDARY CONDITION ENFORCEMENT
// ============================================================================

/// Large number for Dirichlet BC enforcement using penalty method
/// Reference: Babuška, I. (1973). The finite element method with penalty
/// Note: Should be >> 1/tolerance but not cause overflow
pub const DIRICHLET_PENALTY_FACTOR: f64 = 1.0e20;

// ============================================================================
// RELAXATION FACTORS
// ============================================================================

/// Under-relaxation factor for velocity in SIMPLE algorithm
/// Reference: Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow, p. 127
pub const VELOCITY_RELAXATION_FACTOR: f64 = 0.7;

/// Under-relaxation factor for pressure in SIMPLE algorithm
/// Reference: Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow, p. 127
pub const PRESSURE_RELAXATION_FACTOR: f64 = 0.3;

/// Under-relaxation for turbulence quantities
/// Reference: Menter, F.R. (1994). Two-equation eddy-viscosity turbulence models
pub const TURBULENCE_RELAXATION_FACTOR: f64 = 0.8;

// ============================================================================
// TIME STEPPING
// ============================================================================

/// Default CFL number for explicit time stepping
/// Reference: Courant, Friedrichs, Lewy (1928). Über die partiellen Differenzengleichungen
pub const DEFAULT_CFL_NUMBER: f64 = 0.5;

/// Maximum allowable CFL for stability
/// Reference: LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems
pub const MAX_CFL_NUMBER: f64 = 1.0;

/// Minimum time step to prevent numerical issues
pub const MIN_TIME_STEP: f64 = 1e-10;

// ============================================================================
// GRID QUALITY
// ============================================================================

/// Maximum acceptable aspect ratio for grid cells
/// Reference: Mavriplis, D.J. (1997). Unstructured grid techniques
pub const MAX_ASPECT_RATIO: f64 = 100.0;

/// Minimum acceptable cell volume (relative to domain)
pub const MIN_CELL_VOLUME_RATIO: f64 = 1e-12;

// ============================================================================
// PHYSICAL THRESHOLDS
// ============================================================================

/// Threshold for considering flow as incompressible (Mach number)
/// Reference: Anderson, J.D. (1995). Computational Fluid Dynamics, p. 87
pub const INCOMPRESSIBLE_MACH_LIMIT: f64 = 0.3;

/// Minimum Reynolds number for turbulent flow
/// Reference: Schlichting, H. (1979). Boundary Layer Theory, p. 416
pub const TURBULENT_REYNOLDS_THRESHOLD: f64 = 2300.0;

// ============================================================================
// NUMERICAL STABILITY
// ============================================================================

/// Small number to prevent division by zero
/// Note: Should be > machine epsilon but small enough not to affect results
pub const DIVISION_EPSILON: f64 = 1e-30;

/// Threshold for detecting stagnation in iterative methods
pub const STAGNATION_RATIO: f64 = 0.99;

/// Factor for detecting divergence in iterative methods
pub const DIVERGENCE_FACTOR: f64 = 1e10;