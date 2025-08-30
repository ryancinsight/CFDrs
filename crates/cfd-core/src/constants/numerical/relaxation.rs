//! Relaxation factors for iterative solvers
//!
//! Under-relaxation and over-relaxation parameters for various solution variables.

/// Default under-relaxation factor
pub const DEFAULT_RELAXATION: f64 = 1.0;

/// Pressure under-relaxation factor (SIMPLE/PISO)
pub const PRESSURE_RELAXATION: f64 = 0.3;

/// Velocity under-relaxation factor (SIMPLE/PISO)
pub const VELOCITY_RELAXATION: f64 = 0.7;

/// Density under-relaxation factor
pub const DENSITY_RELAXATION: f64 = 0.5;

/// Temperature under-relaxation factor
pub const TEMPERATURE_RELAXATION: f64 = 0.8;

/// Turbulent kinetic energy under-relaxation
pub const TKE_RELAXATION: f64 = 0.5;

/// Turbulent dissipation rate under-relaxation
pub const EPSILON_RELAXATION: f64 = 0.5;

/// Species concentration under-relaxation
pub const SPECIES_RELAXATION: f64 = 0.8;

/// SOR (Successive Over-Relaxation) optimal factor range
/// Minimum SOR relaxation factor (1.0 = Gauss-Seidel)
pub const SOR_OPTIMAL_MIN: f64 = 1.0;
/// Maximum SOR relaxation factor (typically < 2.0 for convergence)
pub const SOR_OPTIMAL_MAX: f64 = 2.0;
/// Default SOR relaxation factor for Poisson equation
pub const SOR_OMEGA_DEFAULT: f64 = 1.8;

/// Gauss-Seidel relaxation factor
pub const GAUSS_SEIDEL_RELAXATION: f64 = 1.0;
