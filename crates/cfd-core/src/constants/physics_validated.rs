//! Physics constants with literature validation
//!
//! All constants are cross-referenced against authoritative sources.
//! Values use SI units unless otherwise specified.

/// Fluid dynamics constants
pub mod fluid_dynamics {
    /// Kinematic viscosity of water at 20°C [m²/s]
    /// Source: White, F.M. (2016). Fluid Mechanics, 8th ed., Table 1.4
    pub const WATER_KINEMATIC_VISCOSITY_20C: f64 = 1.004e-6;

    /// Density of water at 20°C, 1 atm [kg/m³]
    /// Source: NIST Chemistry `WebBook`, SRD 69
    pub const WATER_DENSITY_20C: f64 = 998.2071;

    /// Dynamic viscosity of water at 20°C [Pa·s]
    /// Source: Kestin et al. (1978). J. Phys. Chem. Ref. Data, 7(3), 941-948
    pub const WATER_DYNAMIC_VISCOSITY_20C: f64 = 1.0016e-3;

    /// Density of air at 20°C, 1 atm [kg/m³]
    /// Source: ISO 2533:1975 Standard Atmosphere
    pub const AIR_DENSITY_20C: f64 = 1.2041;

    /// Dynamic viscosity of air at 20°C [Pa·s]
    /// Source: Sutherland, W. (1893). Phil. Mag., 36, 507-531
    pub const AIR_DYNAMIC_VISCOSITY_20C: f64 = 1.8205e-5;
}

/// Reynolds number thresholds
pub mod reynolds {
    /// Critical Reynolds number for pipe flow transition (lower bound)
    /// Source: Schlichting & Gersten (2017). Boundary-Layer Theory, 9th ed.
    pub const PIPE_TRANSITION_LOWER: f64 = 2300.0;

    /// Critical Reynolds number for pipe flow transition (upper bound)
    /// Source: Schlichting & Gersten (2017). Boundary-Layer Theory, 9th ed.
    pub const PIPE_TRANSITION_UPPER: f64 = 4000.0;

    /// Critical Reynolds number for flat plate boundary layer
    /// Source: Anderson, J.D. (2017). Fundamentals of Aerodynamics, 6th ed.
    pub const FLAT_PLATE_CRITICAL: f64 = 5e5;
}

/// Thermodynamic properties
pub mod thermodynamics {
    /// Specific heat of water at 20°C [J/(kg·K)]
    /// Source: NIST Chemistry `WebBook`, SRD 69
    pub const WATER_SPECIFIC_HEAT_20C: f64 = 4181.3;

    /// Thermal conductivity of water at 20°C [W/(m·K)]
    /// Source: Ramires et al. (1995). J. Phys. Chem. Ref. Data, 24(3), 1377-1381
    pub const WATER_THERMAL_CONDUCTIVITY_20C: f64 = 0.5984;

    /// Specific heat of air at 20°C, constant pressure [J/(kg·K)]
    /// Source: Lemmon et al. (2000). J. Phys. Chem. Ref. Data, 29(3), 331-385
    pub const AIR_SPECIFIC_HEAT_CP_20C: f64 = 1005.0;

    /// Thermal conductivity of air at 20°C [W/(m·K)]
    /// Source: Kadoya et al. (1985). J. Phys. Chem. Ref. Data, 14(4), 947-970
    pub const AIR_THERMAL_CONDUCTIVITY_20C: f64 = 0.02514;
}

/// Universal physical constants
pub mod universal {
    /// Universal gas constant [J/(mol·K)]
    /// Source: CODATA 2018, NIST SP 961
    pub const GAS_CONSTANT: f64 = 8.314_462_618;

    /// Standard acceleration due to gravity [m/s²]
    /// Source: ISO 80000-3:2019
    pub const GRAVITY_STANDARD: f64 = 9.80665;

    /// Stefan-Boltzmann constant [W/(m²·K⁴)]
    /// Source: CODATA 2018, NIST SP 961
    pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
}

/// Dimensionless numbers for validation
pub mod validation {
    /// Prandtl number for water at 20°C
    /// Pr = μ·Cp/k = 1.0016e-3 * 4181.3 / 0.5984 = 6.99
    /// Source: Calculated from validated properties above
    pub const WATER_PRANDTL_20C: f64 = 6.99;

    /// Prandtl number for air at 20°C
    /// Pr = μ·Cp/k = 1.8205e-5 * 1005.0 / 0.02514 = 0.728
    /// Source: Calculated from validated properties above
    pub const AIR_PRANDTL_20C: f64 = 0.728;
}
