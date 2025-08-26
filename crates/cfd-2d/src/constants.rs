//! Physical and numerical constants for 2D CFD calculations

/// Numerical constants
pub mod numerical {
    /// Half value for averaging
    pub const HALF: f64 = 0.5;

    /// Double value for central differences
    pub const TWO: f64 = 2.0;

    /// Triple value for higher-order schemes
    pub const THREE: f64 = 3.0;

    /// Machine epsilon for convergence checks
    pub const EPSILON_MACHINE: f64 = 1e-14;

    /// Convergence tolerance
    pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;

    /// Small value for division safety
    pub const EPSILON_DIV: f64 = 1e-10;
}

/// Relaxation factors for iterative methods
pub mod relaxation {
    /// Velocity under-relaxation (SIMPLE algorithm)
    pub const VELOCITY: f64 = 0.7;

    /// Pressure under-relaxation (SIMPLE algorithm)
    pub const PRESSURE: f64 = 0.3;

    /// Temperature under-relaxation
    pub const TEMPERATURE: f64 = 0.8;

    /// Turbulence under-relaxation
    pub const TURBULENCE: f64 = 0.8;
}

/// Turbulence model constants
pub mod turbulence {
    /// von Karman constant
    pub const KAPPA: f64 = cfd_core::constants::fluid::VON_KARMAN;

    /// Wall function constant
    pub const E_WALL: f64 = cfd_core::constants::fluid::WALL_FUNCTION_E;

    /// k-epsilon model constant `C_mu`
    pub const C_MU: f64 = 0.09;

    /// k-epsilon model constant `C_1e`
    pub const C1_EPSILON: f64 = 1.44;

    /// k-epsilon model constant `C_2e`
    pub const C2_EPSILON: f64 = 1.92;

    /// Turbulent Prandtl number for k
    pub const SIGMA_K: f64 = 1.0;

    /// Turbulent Prandtl number for epsilon
    pub const SIGMA_EPSILON: f64 = 1.3;

    /// Minimum turbulent kinetic energy
    pub const K_MIN: f64 = 1e-10;

    /// Minimum dissipation rate
    pub const EPSILON_MIN: f64 = 1e-10;
}

/// CFL conditions for stability
pub mod cfl {
    /// Maximum CFL number for explicit schemes
    pub const EXPLICIT_MAX: f64 = 0.5;

    /// Maximum CFL number for implicit schemes
    pub const IMPLICIT_MAX: f64 = 1.0;

    /// Safety factor for adaptive time stepping
    pub const SAFETY_FACTOR: f64 = 0.9;
}
