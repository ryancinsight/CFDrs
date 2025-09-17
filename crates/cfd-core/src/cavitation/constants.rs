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

/// Nurick correlation coefficient for cavity length
pub const NURICK_K_COEFFICIENT: f64 = 0.88;

/// Nurick correlation exponent
pub const NURICK_EXPONENT: f64 = 0.5;

/// Incipient cavitation number for typical venturi
pub const SIGMA_INCIPIENT: f64 = 1.2;

/// Critical cavitation number for choked flow
pub const SIGMA_CRITICAL: f64 = 1.0;

/// Plesset-Chapman erosion model constants
/// Material erosion coefficient for steel (m³/N²·Hz·s)
pub const EROSION_COEFFICIENT_STEEL: f64 = 5e-13;

/// Pressure exponent for erosion (typical range 2.0-2.5)
pub const EROSION_PRESSURE_EXPONENT: f64 = 2.25;

/// Basquin's law constants for fatigue
/// Fatigue strength coefficient ratio (`σ_f`'/UTS)
pub const FATIGUE_STRENGTH_RATIO: f64 = 0.9;

/// Basquin exponent for steels (typical -0.085 to -0.12)
pub const BASQUIN_EXPONENT_STEEL: f64 = -0.087;
