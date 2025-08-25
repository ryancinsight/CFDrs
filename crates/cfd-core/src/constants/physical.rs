//! Physical constants for CFD simulations
//!
//! This module provides standardized physical constants to ensure
//! consistency across the codebase (SSOT principle).

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Temperature conversion constants
pub mod temperature {
    /// Celsius to Kelvin offset
    pub const CELSIUS_TO_KELVIN_OFFSET: f64 = 273.15;

    /// Fahrenheit to Celsius conversion factor (5/9)
    pub const FAHRENHEIT_TO_CELSIUS_FACTOR: f64 = 5.0 / 9.0;

    /// Fahrenheit zero point offset
    pub const FAHRENHEIT_ZERO_OFFSET: f64 = 32.0;
}

/// Pressure conversion constants
pub mod pressure {
    /// Standard atmospheric pressure in Pascals
    pub const STANDARD_ATMOSPHERE_PA: f64 = 101_325.0;

    /// Bar to Pascal conversion factor
    pub const BAR_TO_PASCAL: f64 = 100_000.0;

    /// PSI to Pascal conversion factor
    pub const PSI_TO_PASCAL: f64 = 6_894.757;
}

/// Fluid properties at standard conditions
pub mod fluid {
    /// Water density at standard conditions (kg/m³)
    pub const WATER_DENSITY_STD: f64 = 998.2;

    /// Air density at standard conditions (kg/m³)
    pub const AIR_DENSITY_STD: f64 = 1.225;

    /// Water dynamic viscosity at 20°C (Pa·s)
    pub const WATER_VISCOSITY_20C: f64 = 1.002e-3;

    /// Air dynamic viscosity at 20°C (Pa·s)
    pub const AIR_VISCOSITY_20C: f64 = 1.81e-5;

    /// Water kinematic viscosity at 20°C (m²/s)
    pub const WATER_KINEMATIC_VISCOSITY_20C: f64 = 1.004e-6;

    /// Air kinematic viscosity at 20°C (m²/s)
    pub const AIR_KINEMATIC_VISCOSITY_20C: f64 = 1.48e-5;
}

/// Universal physical constants
pub mod universal {
    /// Universal gas constant (J/(mol·K))
    pub const GAS_CONSTANT: f64 = 8.314_462_618;

    /// Standard gravity (m/s²)
    pub const STANDARD_GRAVITY: f64 = 9.806_65;

    /// Stefan-Boltzmann constant (W/(m²·K⁴))
    pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
}

/// Get a physical constant as a generic RealField type
pub fn get_constant<T: RealField + FromPrimitive>(value: f64) -> T {
    T::from_f64(value).unwrap_or_else(T::one)
}
