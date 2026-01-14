//! Centralized constants for CFD simulations.
//!
//! This module provides a single source of truth for all mathematical, physical,
//! and numerical constants used throughout the CFDrs suite.

/// Mathematical constants
pub mod mathematical {
    use std::f64::consts;

    /// Pi (π)
    pub const PI: f64 = consts::PI;
    /// 2π
    pub const TWO_PI: f64 = 2.0 * consts::PI;
    /// π/2
    pub const PI_OVER_TWO: f64 = consts::PI / 2.0;
    /// π/4
    pub const PI_OVER_FOUR: f64 = consts::PI / 4.0;
    /// Euler's number (e)
    pub const E: f64 = consts::E;
    /// √2
    pub const SQRT_TWO: f64 = consts::SQRT_2;
    /// ln(2)
    pub const LN_TWO: f64 = consts::LN_2;

    /// Common fractions and integers as f64
    pub mod numeric {
        /// Zero (0.0)
        pub const ZERO: f64 = 0.0;
        /// One (1.0)
        pub const ONE: f64 = 1.0;
        /// Two (2.0)
        pub const TWO: f64 = 2.0;
        /// Three (3.0)
        pub const THREE: f64 = 3.0;
        /// Four (4.0)
        pub const FOUR: f64 = 4.0;
        /// Five (5.0)
        pub const FIVE: f64 = 5.0;
        /// Six (6.0)
        pub const SIX: f64 = 6.0;
        /// Eight (8.0)
        pub const EIGHT: f64 = 8.0;
        /// Ten (10.0)
        pub const TEN: f64 = 10.0;
        /// One half (0.5)
        pub const ONE_HALF: f64 = 0.5;
        /// One third (1/3)
        pub const ONE_THIRD: f64 = 1.0 / 3.0;
        /// Two thirds (2/3)
        pub const TWO_THIRDS: f64 = 2.0 / 3.0;
        /// One quarter (0.25)
        pub const ONE_QUARTER: f64 = 0.25;
        /// Three quarters (0.75)
        pub const THREE_QUARTERS: f64 = 0.75;
        /// One sixth (1/6)
        pub const ONE_SIXTH: f64 = 1.0 / 6.0;
        /// Five sixths (5/6)
        pub const FIVE_SIXTHS: f64 = 5.0 / 6.0;
        /// One twelfth (1/12)
        pub const ONE_TWELFTH: f64 = 1.0 / 12.0;
    }
}

/// Physical constants and fluid properties
pub mod physics {
    /// Universal physical constants
    pub mod universal {
        /// Universal gas constant [J/(mol·K)]
        pub const GAS_CONSTANT: f64 = 8.314_462_618;
        /// Standard acceleration due to gravity [m/s²]
        pub const GRAVITY: f64 = 9.806_65;
        /// Stefan-Boltzmann constant [W/(m²·K⁴)]
        pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
    }

    /// Thermodynamic constants
    pub mod thermo {
        /// Specific gas constant for air [J/(kg·K)]
        pub const R_AIR: f64 = 287.052_874;
        /// Ratio of specific heats for air (γ)
        pub const GAMMA_AIR: f64 = 1.4;
        /// Standard atmospheric pressure [Pa]
        pub const P_ATM: f64 = 101_325.0;
        /// Standard temperature [K]
        pub const T_STANDARD: f64 = 288.15;
        /// Celsius to Kelvin conversion offset
        pub const CELSIUS_TO_KELVIN: f64 = 273.15;
    }

    /// Fluid properties at 20°C and 1 atm
    pub mod fluid {
        /// Water density [kg/m³]
        pub const WATER_DENSITY: f64 = 998.2;
        /// Water dynamic viscosity [Pa·s]
        pub const WATER_VISCOSITY: f64 = 1.002e-3;
        /// Water kinematic viscosity [m²/s]
        pub const WATER_KINEMATIC_VISCOSITY: f64 = 1.004e-6;
        /// Water specific heat [J/(kg·K)]
        pub const WATER_SPECIFIC_HEAT: f64 = 4182.0;
        /// Water thermal conductivity [W/(m·K)]
        pub const WATER_THERMAL_CONDUCTIVITY: f64 = 0.606;
        /// Water speed of sound [m/s]
        pub const WATER_SPEED_OF_SOUND: f64 = 1482.0;

        /// Air density [kg/m³]
        pub const AIR_DENSITY: f64 = 1.204;
        /// Air dynamic viscosity [Pa·s]
        pub const AIR_VISCOSITY: f64 = 1.825e-5;
        /// Air kinematic viscosity [m²/s]
        pub const AIR_KINEMATIC_VISCOSITY: f64 = 1.51e-5;
        /// Air specific heat [J/(kg·K)]
        pub const AIR_SPECIFIC_HEAT: f64 = 1005.0;
        /// Air thermal conductivity [W/(m·K)]
        pub const AIR_THERMAL_CONDUCTIVITY: f64 = 0.0257;
        /// Air speed of sound [m/s]
        pub const AIR_SPEED_OF_SOUND: f64 = 343.2;

        /// Von Kármán constant
        pub const VON_KARMAN: f64 = 0.41;
    }

    /// Turbulence model constants
    pub mod turbulence {
        /// k-epsilon model constants (Launder & Spalding 1974)
        /// C_mu constant for eddy viscosity calculation
        pub const K_EPSILON_C_MU: f64 = 0.09;
        /// C1 constant for epsilon production
        pub const K_EPSILON_C1: f64 = 1.44;
        /// C2 constant for epsilon destruction
        pub const K_EPSILON_C2: f64 = 1.92;
        /// Turbulent Prandtl number for k
        pub const K_EPSILON_SIGMA_K: f64 = 1.0;
        /// Turbulent Prandtl number for epsilon
        pub const K_EPSILON_SIGMA_EPSILON: f64 = 1.3;
    }

    /// Hydraulics and pipe flow constants
    pub mod hydraulics {
        /// Blasius equation limits and coefficients
        /// Maximum Reynolds number for Blasius correlation validity
        pub const BLASIUS_MAX_RE: f64 = 1e5;
        /// Blasius friction factor coefficient (0.3164 / Re^0.25)
        pub const BLASIUS_COEFFICIENT: f64 = 0.3164;
        /// Blasius friction factor exponent
        pub const BLASIUS_EXPONENT: f64 = 0.25;

        /// Colebrook-White equation constants
        /// Constant divisor for relative roughness (3.7)
        pub const COLEBROOK_ROUGHNESS_DIVISOR: f64 = 3.7;
        /// Reynolds term numerator (2.51)
        pub const COLEBROOK_REYNOLDS_NUMERATOR: f64 = 2.51;

        /// Haaland equation constants
        /// Logarithmic coefficient (-1.8)
        pub const HAALAND_LOG_COEFFICIENT: f64 = -1.8;
        /// Roughness exponent (1.11)
        pub const HAALAND_ROUGHNESS_EXPONENT: f64 = 1.11;
    }

    /// Dimensionless numbers constants
    pub mod dimensionless {
        /// Reynolds number limits
        pub mod reynolds {
            /// Lower critical Reynolds number for pipe flow (laminar limit)
            pub const PIPE_LAMINAR_MAX: f64 = 2300.0;
            /// Upper critical Reynolds number for pipe flow (turbulent onset)
            pub const PIPE_TURBULENT_MIN: f64 = 4000.0;
            /// Critical Reynolds number for flow over a flat plate
            pub const PLATE_CRITICAL: f64 = 500_000.0;
        }

        /// Mach number limits and regimes
        pub mod mach {
            /// Limit for incompressible flow assumption (Ma < 0.3)
            pub const INCOMPRESSIBLE_LIMIT: f64 = 0.3;
            /// Lower limit for transonic regime (Ma ~ 0.8)
            pub const TRANSONIC_LOWER: f64 = 0.8;
            /// Upper limit for transonic regime (Ma ~ 1.2)
            pub const TRANSONIC_UPPER: f64 = 1.2;
            /// Limit for hypersonic flow (Ma > 5.0)
            pub const HYPERSONIC: f64 = 5.0;
        }
    }
}

/// Numerical and solver constants
pub mod numerical {
    /// Convergence and iteration parameters
    pub mod solver {
        /// Default relative convergence tolerance
        pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;
        /// Small value for denominator protection
        pub const EPSILON_TOLERANCE: f64 = 1e-10;
        /// Maximum outer iterations for non-linear systems
        pub const MAX_ITERATIONS_OUTER: usize = 100;
        /// Maximum inner iterations for linear solvers
        pub const MAX_ITERATIONS_INNER: usize = 1000;
    }

    /// Time integration parameters
    pub mod time {
        /// Default Courant-Friedrichs-Lewy number
        pub const DEFAULT_CFL: f64 = 0.5;
        /// Maximum stable CFL number
        pub const MAX_CFL: f64 = 1.0;
        /// Default time step [s]
        pub const DEFAULT_TIME_STEP: f64 = 0.001;
        /// Minimum allowable time step [s]
        pub const MIN_TIME_STEP: f64 = 1e-12;
    }

    /// Relaxation factors
    pub mod relaxation {
        /// Under-relaxation for pressure correction
        pub const PRESSURE_RELAXATION: f64 = 0.3;
        /// Under-relaxation for momentum equations
        pub const VELOCITY_RELAXATION: f64 = 0.7;
    }

    /// Machine epsilon and small values
    pub mod precision {
        /// Small value for numerical stability
        pub const SMALL: f64 = 1e-10;
        /// Large value for numerical stability
        pub const LARGE: f64 = 1e10;
        /// Machine epsilon for f64
        pub const EPSILON: f64 = f64::EPSILON;
    }
}
