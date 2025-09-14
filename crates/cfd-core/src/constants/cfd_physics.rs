//! Comprehensive CFD physics constants
//!
//! All physical constants used in CFD simulations to eliminate magic numbers
//! Based on standard literature values

use std::f64::consts;

// ============================================================================
// FLUID PROPERTIES
// ============================================================================

/// Water properties at 20°C
pub mod water {
    /// Density [kg/m³] - White (2011) Fluid Mechanics
    pub const DENSITY: f64 = 998.2;

    /// Dynamic viscosity [Pa·s] - White (2011)
    pub const DYNAMIC_VISCOSITY: f64 = 1.002e-3;

    /// Kinematic viscosity [m²/s]
    pub const KINEMATIC_VISCOSITY: f64 = 1.004e-6;

    /// Bulk modulus [Pa]
    pub const BULK_MODULUS: f64 = 2.2e9;

    /// Speed of sound [m/s]
    pub const SPEED_OF_SOUND: f64 = 1483.0;

    /// Surface tension with air [N/m]
    pub const SURFACE_TENSION: f64 = 0.0728;

    /// Thermal conductivity [W/(m·K)]
    pub const THERMAL_CONDUCTIVITY: f64 = 0.598;

    /// Specific heat capacity [J/(kg·K)]
    pub const SPECIFIC_HEAT: f64 = 4182.0;
}

/// Air properties at 20°C, 1 atm
pub mod air {
    /// Density [kg/m³] - Anderson (2011) Fundamentals of Aerodynamics
    pub const DENSITY: f64 = 1.204;

    /// Dynamic viscosity [Pa·s]
    pub const DYNAMIC_VISCOSITY: f64 = 1.82e-5;

    /// Kinematic viscosity [m²/s]
    pub const KINEMATIC_VISCOSITY: f64 = 1.51e-5;

    /// Specific gas constant [J/(kg·K)]
    pub const GAS_CONSTANT: f64 = 287.05;

    /// Specific heat ratio (γ = Cp/Cv)
    pub const HEAT_RATIO: f64 = 1.4;

    /// Speed of sound [m/s]
    pub const SPEED_OF_SOUND: f64 = 343.0;
}

// ============================================================================
// DIMENSIONLESS NUMBERS
// ============================================================================

/// Reynolds number ranges for different flow regimes
pub mod reynolds {
    /// Laminar flow upper limit for pipes - White (2011)
    pub const PIPE_LAMINAR_MAX: f64 = 2300.0;

    /// Turbulent flow lower limit for pipes
    pub const PIPE_TURBULENT_MIN: f64 = 4000.0;

    /// Transition zone midpoint
    pub const PIPE_TRANSITION: f64 = 3150.0;

    /// Critical Reynolds for flat plate - Schlichting (1979)
    pub const FLAT_PLATE_CRITICAL: f64 = 5e5;
}

/// Prandtl number values for common fluids
pub mod prandtl {
    /// Water at 20°C
    pub const WATER: f64 = 7.0;

    /// Air at 20°C
    pub const AIR: f64 = 0.71;

    /// Liquid metals typical value
    pub const LIQUID_METAL: f64 = 0.01;
}

// ============================================================================
// TURBULENCE CONSTANTS
// ============================================================================

/// Turbulence modeling constants from fluid mechanics literature
pub mod turbulence {
    /// von Kármán constant - Pope (2000)
    pub const VON_KARMAN: f64 = 0.41;

    /// Log-law constant B - Pope (2000)
    pub const LOG_LAW_B: f64 = 5.2;

    /// Wall y+ values
    pub mod y_plus {
        /// Viscous sublayer limit
        pub const VISCOUS_SUBLAYER: f64 = 5.0;

        /// Buffer layer start
        pub const BUFFER_START: f64 = 5.0;

        /// Buffer layer end
        pub const BUFFER_END: f64 = 30.0;

        /// Log layer start
        pub const LOG_LAYER: f64 = 30.0;

        /// Recommended first cell y+
        pub const FIRST_CELL_TARGET: f64 = 1.0;
    }

    /// k-ε turbulence model constants from Launder & Spalding (1974)
    pub mod k_epsilon {
        /// Model constant Cμ in turbulent viscosity calculation
        pub const C_MU: f64 = 0.09;
        /// Model constant C₁ε in ε transport equation
        pub const C_1E: f64 = 1.44;
        /// Model constant C₂ε in ε transport equation  
        pub const C_2E: f64 = 1.92;
        /// Turbulent Prandtl number for k equation
        pub const SIGMA_K: f64 = 1.0;
        /// Turbulent Prandtl number for ε equation
        pub const SIGMA_E: f64 = 1.3;
    }

    /// k-ω SST turbulence model constants from Menter (1994)
    pub mod k_omega_sst {
        /// Model constant β* in specific dissipation rate calculation
        pub const BETA_STAR: f64 = 0.09;
        /// Turbulent Prandtl number for k equation (inner region)
        pub const SIGMA_K1: f64 = 0.85;
        /// Turbulent Prandtl number for k equation (outer region)
        pub const SIGMA_K2: f64 = 1.0;
        /// Turbulent Prandtl number for ω equation (inner region)
        pub const SIGMA_W1: f64 = 0.5;
        /// Turbulent Prandtl number for ω equation (outer region)
        pub const SIGMA_W2: f64 = 0.856;
        /// Model constant β₁ for inner region
        pub const BETA_1: f64 = 0.075;
        /// Model constant β₂ for outer region
        pub const BETA_2: f64 = 0.0828;
        /// Model constant α₁ for transition
        pub const ALPHA_1: f64 = 0.31;
        /// Model constant γ₁ for production term
        pub const GAMMA_1: f64 = 0.553;
        /// Model constant γ₂ for production term
        pub const GAMMA_2: f64 = 0.440;
    }
}

// ============================================================================
// NUMERICAL CONSTANTS
// ============================================================================

/// Numerical constants for computational methods
pub mod numerical {
    /// CFL number limits
    pub mod cfl {
        /// Explicit schemes maximum
        pub const EXPLICIT_MAX: f64 = 1.0;

        /// Implicit schemes typical
        pub const IMPLICIT_TYPICAL: f64 = 10.0;

        /// Conservative value for stability
        pub const CONSERVATIVE: f64 = 0.5;
    }

    /// Convergence criteria
    pub mod convergence {
        /// Default residual tolerance
        pub const RESIDUAL_TOLERANCE: f64 = 1e-6;

        /// Machine epsilon safety factor
        pub const EPSILON_SAFETY: f64 = 1e-14;

        /// Maximum iterations default
        pub const MAX_ITERATIONS_DEFAULT: usize = 1000;
    }

    /// Relaxation factors - Patankar (1980)
    pub mod relaxation {
        /// Pressure equation
        pub const PRESSURE: f64 = 0.3;

        /// Momentum equations
        pub const MOMENTUM: f64 = 0.7;

        /// Turbulence equations
        pub const TURBULENCE: f64 = 0.7;

        /// Temperature equation
        pub const TEMPERATURE: f64 = 0.9;
    }
}

// ============================================================================
// LATTICE BOLTZMANN CONSTANTS
// ============================================================================

/// Lattice Boltzmann method constants
pub mod lattice_boltzmann {
    /// D2Q9 lattice sound speed squared (cs²)
    pub const LATTICE_SOUND_SPEED_SQUARED: f64 = 1.0 / 3.0;

    /// D2Q9 lattice sound speed
    pub const LATTICE_SOUND_SPEED: f64 = 0.577_350_269; // 1/√3

    /// BGK relaxation time offset
    pub const RELAXATION_TIME_OFFSET: f64 = 0.5;

    /// Reference temperature for thermal LBM
    pub const REFERENCE_TEMPERATURE: f64 = 1.0;
}

// ============================================================================
// CAVITATION CONSTANTS
// ============================================================================

/// Cavitation modeling constants
pub mod cavitation {
    /// Vapor pressure of water at 20°C [Pa]
    pub const WATER_VAPOR_PRESSURE: f64 = 2339.0;

    /// Critical pressure ratio for cavitation inception
    pub const CRITICAL_PRESSURE_RATIO: f64 = 0.9;

    /// Rayleigh-Plesset constants
    pub mod rayleigh_plesset {
        /// Polytropic exponent for gas
        pub const POLYTROPIC_EXPONENT: f64 = 1.4;

        /// Initial gas content fraction
        pub const INITIAL_GAS_FRACTION: f64 = 1e-6;
    }
}

// ============================================================================
// PHYSICAL CONSTANTS
// ============================================================================

/// Universal constants
pub mod universal {
    /// Standard gravity [m/s²]
    pub const GRAVITY: f64 = 9.80665;

    /// Standard atmospheric pressure [Pa]
    pub const ATMOSPHERIC_PRESSURE: f64 = 101_325.0;

    /// Universal gas constant [J/(mol·K)]
    pub const GAS_CONSTANT: f64 = 8.314462618;

    /// Stefan-Boltzmann constant [W/(m²·K⁴)]
    pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;
}

// ============================================================================
// GEOMETRIC CONSTANTS
// ============================================================================

/// Geometric constants for CFD calculations
pub mod geometry {
    use super::consts;

    /// Pi
    pub const PI: f64 = consts::PI;

    /// Two pi
    pub const TWO_PI: f64 = 2.0 * consts::PI;

    /// Pi over two
    pub const HALF_PI: f64 = consts::PI / 2.0;

    /// Pi over four
    pub const QUARTER_PI: f64 = consts::PI / 4.0;

    /// Degrees to radians conversion
    pub const DEG_TO_RAD: f64 = consts::PI / 180.0;

    /// Radians to degrees conversion
    pub const RAD_TO_DEG: f64 = 180.0 / consts::PI;
}
