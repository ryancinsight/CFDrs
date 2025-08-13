//! Mathematical and numerical constants used throughout the CFD suite
//!
//! This module provides named constants to eliminate magic numbers and ensure
//! consistency across the codebase. All constants follow the SSOT principle.

use nalgebra::RealField;

/// Common numerical factors
pub mod factors {
    /// Factor of 2.0
    pub const TWO: f64 = 2.0;
    
    /// Factor of 3.0
    pub const THREE: f64 = 3.0;
    
    /// Factor of 4.0
    pub const FOUR: f64 = 4.0;
    
    /// Factor of 0.5 (half)
    pub const HALF: f64 = 0.5;
    
    /// Factor of 0.25 (quarter)
    pub const QUARTER: f64 = 0.25;
    
    /// Factor of 10.0
    pub const TEN: f64 = 10.0;
    
    /// Factor of 20.0
    pub const TWENTY: f64 = 20.0;
}

/// Finite element method constants
pub mod fem {
    /// Linear shape function value at centroid for tetrahedral element
    pub const LINEAR_SHAPE_CENTROID_WEIGHT: f64 = 0.25;
    
    /// Consistent mass matrix factor for tetrahedral element
    pub const CONSISTENT_MASS_FACTOR: f64 = 0.05; // 1/20
    
    /// Diagonal mass matrix factor
    pub const DIAGONAL_MASS_FACTOR: f64 = 2.0;
    
    /// Stability reduction factor for stiffness matrix
    pub const STABILITY_FACTOR: f64 = 0.1;
}

/// Channel flow constants
pub mod channel {
    /// Entrance length coefficient for laminar flow (Re < 10)
    pub const ENTRANCE_LENGTH_LAMINAR_LOW: f64 = 0.06;
    
    /// Entrance length coefficient for transitional flow
    pub const ENTRANCE_LENGTH_TRANSITIONAL: f64 = 4.4;
    
    /// Entrance length exponent for transitional flow
    pub const ENTRANCE_LENGTH_EXPONENT: f64 = 1.0 / 6.0;
    
    /// Entrance length coefficient for turbulent flow
    pub const ENTRANCE_LENGTH_TURBULENT: f64 = 0.1;
    
    /// Ellipse perimeter approximation constants
    pub const ELLIPSE_PERIMETER_FACTOR_1: f64 = 3.0;
    pub const ELLIPSE_PERIMETER_FACTOR_2: f64 = 10.0;
    pub const ELLIPSE_PERIMETER_FACTOR_3: f64 = 4.0;
}

/// Physical property constants
pub mod physical {
    /// Default water density at 20°C [kg/m³]
    pub const WATER_DENSITY_20C: f64 = 998.2;
    
    /// Default water viscosity at 20°C [Pa·s]
    pub const WATER_VISCOSITY_20C: f64 = 0.001;
    
    /// Default air density at 20°C, 1 atm [kg/m³]
    pub const AIR_DENSITY_20C: f64 = 1.204;
    
    /// Default air viscosity at 20°C [Pa·s]
    pub const AIR_VISCOSITY_20C: f64 = 1.81e-5;
    
    /// Standard temperature in Kelvin (20°C)
    pub const STANDARD_TEMPERATURE_K: f64 = 293.15;
}

/// Helper trait to convert constants to generic RealField types
pub trait ConstantConvert<T: RealField> {
    /// Convert the constant to the target type
    fn to_field(value: f64) -> T;
}

impl<T: RealField> ConstantConvert<T> for f64 {
    fn to_field(value: f64) -> T {
        T::from_f64(value).unwrap_or_else(T::one)
    }
}

/// Macro to easily convert constants to RealField types
#[macro_export]
macro_rules! const_field {
    ($const_path:path, $T:ty) => {
        <$T>::from_f64($const_path).unwrap_or_else(<$T>::one)
    };
}