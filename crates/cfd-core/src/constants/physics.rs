//! Physics constants for CFD computations
//!
//! All physical constants used throughout the codebase, organized by domain.
//! References: NIST, ISO standards, and standard CFD literature.

/// Fluid mechanics constants
pub mod fluid {
    /// Von Kármán constant for wall functions
    pub const VON_KARMAN: f64 = 0.41;

    /// Wall function constant E
    pub const WALL_FUNCTION_E: f64 = 9.8;

    /// Y+ threshold for laminar sublayer
    pub const Y_PLUS_LAMINAR: f64 = 11.63;

    /// Molecular Prandtl number for air
    pub const PRANDTL_AIR: f64 = 0.72;

    /// Molecular Prandtl number for water
    pub const PRANDTL_WATER: f64 = 7.0;

    /// Sutherland's law reference temperature [K]
    pub const SUTHERLAND_T_REF: f64 = 273.15;

    /// Sutherland's law constant for air [K]
    pub const SUTHERLAND_S_AIR: f64 = 110.4;

    /// Parabolic velocity profile coefficient for channel flow
    pub const CHANNEL_FLOW_COEFFICIENT: f64 = 4.0;

    /// Pipe flow shape factor (circular cross-section)
    pub const PIPE_SHAPE_FACTOR: f64 = 2.0;
}

/// Thermodynamics constants
pub mod thermo {
    /// Universal gas constant [J/(mol·K)]
    pub const R_UNIVERSAL: f64 = 8.314462618;

    /// Specific gas constant for air [J/(kg·K)]
    pub const R_AIR: f64 = 287.052874;

    /// Ratio of specific heats for air (γ)
    pub const GAMMA_AIR: f64 = 1.4;

    /// Standard atmospheric pressure [Pa]
    pub const P_ATM: f64 = 101325.0;

    /// Standard temperature [K]
    pub const T_STANDARD: f64 = 288.15;

    /// Celsius to Kelvin conversion offset
    pub const CELSIUS_TO_KELVIN: f64 = 273.15;

    /// Stefan-Boltzmann constant [W/(m²·K⁴)]
    pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;
}

/// Turbulence model constants
pub mod turbulence {
    /// k-ε model constants (Launder & Spalding 1974)
    /// C_μ constant for eddy viscosity calculation
    pub const K_EPSILON_C_MU: f64 = 0.09;
    /// C₁ε constant for production term
    pub const K_EPSILON_C1: f64 = 1.44;
    /// C₂ε constant for dissipation term
    pub const K_EPSILON_C2: f64 = 1.92;
    /// σ_k turbulent Prandtl number for k
    pub const K_EPSILON_SIGMA_K: f64 = 1.0;
    /// σ_ε turbulent Prandtl number for ε
    pub const K_EPSILON_SIGMA_EPSILON: f64 = 1.3;

    /// k-ω SST model constants (Menter 1994)
    /// α₁ constant for SST limiter
    pub const K_OMEGA_ALPHA1: f64 = 0.31;
    /// β₁ constant for inner layer
    pub const K_OMEGA_BETA1: f64 = 0.075;
    /// β₂ constant for outer layer
    pub const K_OMEGA_BETA2: f64 = 0.0828;
    /// σ_k1 diffusion constant for k (inner)
    pub const K_OMEGA_SIGMA_K1: f64 = 0.85;
    /// σ_k2 diffusion constant for k (outer)
    pub const K_OMEGA_SIGMA_K2: f64 = 1.0;
    /// σ_ω1 diffusion constant for ω (inner)
    pub const K_OMEGA_SIGMA_W1: f64 = 0.5;
    /// σ_ω2 diffusion constant for ω (outer)
    pub const K_OMEGA_SIGMA_W2: f64 = 0.856;

    /// Spalart-Allmaras model constants
    /// cb1 production constant
    pub const SA_CB1: f64 = 0.1355;
    /// cb2 destruction constant
    pub const SA_CB2: f64 = 0.622;
    /// cv1 viscosity constant
    pub const SA_CV1: f64 = 7.1;
    /// cw1 wall destruction constant
    pub const SA_CW1: f64 = 3.24;
    /// cw2 wall destruction constant
    pub const SA_CW2: f64 = 0.3;
    /// cw3 wall destruction constant
    pub const SA_CW3: f64 = 2.0;
    /// Spalart-Allmaras model sigma constant
    pub const SA_SIGMA: f64 = 2.0 / 3.0;
}

/// Dimensionless numbers thresholds
pub mod dimensionless {
    /// Critical Reynolds numbers for different geometries
    pub mod reynolds {
        /// Pipe flow transition start
        pub const PIPE_CRITICAL_LOWER: f64 = 2300.0;

        /// Pipe flow fully turbulent
        pub const PIPE_CRITICAL_UPPER: f64 = 4000.0;

        /// Flat plate boundary layer transition
        pub const PLATE_CRITICAL: f64 = 5e5;

        /// Sphere drag crisis
        pub const SPHERE_CRITICAL: f64 = 3e5;

        /// Cylinder vortex shedding onset
        pub const CYLINDER_VORTEX_ONSET: f64 = 40.0;

        /// Cylinder drag crisis
        pub const CYLINDER_CRITICAL: f64 = 2e5;
    }

    /// Critical Mach numbers
    pub mod mach {
        /// Incompressible flow limit
        pub const INCOMPRESSIBLE_LIMIT: f64 = 0.3;

        /// Transonic regime start
        pub const TRANSONIC_LOWER: f64 = 0.8;

        /// Transonic regime end
        pub const TRANSONIC_UPPER: f64 = 1.2;

        /// Hypersonic regime start
        pub const HYPERSONIC: f64 = 5.0;
    }

    /// Critical Froude numbers
    pub mod froude {
        /// Subcritical flow
        pub const SUBCRITICAL: f64 = 1.0;

        /// Hydraulic jump threshold
        pub const JUMP_THRESHOLD: f64 = 1.7;
    }
}

/// Mathematical constants used in physics
pub mod math {
    use std::f64::consts;

    /// π
    pub const PI: f64 = consts::PI;

    /// 2π
    pub const TWO_PI: f64 = 2.0 * consts::PI;

    /// π/2
    pub const HALF_PI: f64 = consts::PI / 2.0;

    /// √2
    pub const SQRT_2: f64 = consts::SQRT_2;

    /// Natural logarithm base e
    pub const E: f64 = consts::E;

    /// Golden ratio
    pub const PHI: f64 = 1.618033988749895;
}
