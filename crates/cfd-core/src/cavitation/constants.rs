//! Physical constants for cavitation.

/// Surface tension of water at 20°C (N/m)
pub const SURFACE_TENSION_WATER: f64 = 0.0728;

/// Vapor pressure of water at 20°C (Pa)
pub const VAPOR_PRESSURE_WATER_20C: f64 = 2339.0;

/// Saturation pressure ratio threshold for inception
pub const CAVITATION_INCEPTION_THRESHOLD: f64 = 0.9;

/// Blake critical radius coefficient
pub const BLAKE_CRITICAL_COEFFICIENT: f64 = 0.85;

/// Minimum bubble radius for numerical stability (m)
pub const MIN_BUBBLE_RADIUS: f64 = 1e-9;

/// Maximum void fraction before numerical issues
pub const MAX_VOID_FRACTION: f64 = 0.999;

/// Nucleation site density in clean water (#/m³)
pub const NUCLEATION_DENSITY_CLEAN: f64 = 1e6;

/// Nucleation site density in technical water (#/m³)
pub const NUCLEATION_DENSITY_TECHNICAL: f64 = 1e8;
