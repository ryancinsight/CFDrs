//! Temperature value object

use crate::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::fmt;

// Temperature conversion constants
const KELVIN_TO_CELSIUS_OFFSET: f64 = 273.15;
const FAHRENHEIT_SCALE_FACTOR: f64 = 9.0 / 5.0;
const FAHRENHEIT_OFFSET: f64 = 32.0;
const RANKINE_TO_KELVIN_FACTOR: f64 = 5.0 / 9.0;

/// Temperature value with unit conversions
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Temperature<T: RealField + Copy> {
    /// Value in Kelvin (SI unit)
    kelvin: T,
}

impl<T: RealField + Copy + FromPrimitive> Temperature<T> {
    /// Create temperature in Kelvin
    pub fn from_kelvin(value: T) -> Result<Self> {
        if value < T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Temperature in Kelvin cannot be negative".into(),
            ));
        }
        Ok(Self { kelvin: value })
    }

    /// Create temperature in Celsius
    pub fn from_celsius(value: T) -> Result<Self> {
        let offset = T::from_f64(KELVIN_TO_CELSIUS_OFFSET).unwrap_or_else(T::zero);
        let kelvin = value + offset;
        if kelvin < T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Temperature below absolute zero".into(),
            ));
        }
        Ok(Self { kelvin })
    }

    /// Create temperature in Fahrenheit
    pub fn from_fahrenheit(value: T) -> Result<Self> {
        let scale = T::from_f64(FAHRENHEIT_SCALE_FACTOR).unwrap_or_else(T::one);
        let f_offset = T::from_f64(FAHRENHEIT_OFFSET).unwrap_or_else(T::zero);
        let k_offset = T::from_f64(KELVIN_TO_CELSIUS_OFFSET).unwrap_or_else(T::zero);

        let celsius = (value - f_offset) / scale;
        let kelvin = celsius + k_offset;

        if kelvin < T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Temperature below absolute zero".into(),
            ));
        }
        Ok(Self { kelvin })
    }

    /// Get temperature in Kelvin
    pub fn kelvin(&self) -> T {
        self.kelvin
    }

    /// Get temperature in Celsius
    pub fn celsius(&self) -> T {
        self.kelvin - T::from_f64(KELVIN_TO_CELSIUS_OFFSET).unwrap_or_else(T::zero)
    }

    /// Get temperature in Fahrenheit
    pub fn fahrenheit(&self) -> T {
        let scale = T::from_f64(FAHRENHEIT_SCALE_FACTOR).unwrap_or_else(T::one);
        let f_offset = T::from_f64(FAHRENHEIT_OFFSET).unwrap_or_else(T::zero);
        self.celsius() * scale + f_offset
    }

    /// Get temperature in Rankine
    pub fn rankine(&self) -> T {
        let factor = T::from_f64(RANKINE_TO_KELVIN_FACTOR).unwrap_or_else(T::one);
        self.kelvin / factor
    }

    /// Check if temperature is at standard conditions (20°C)
    pub fn is_standard(&self) -> bool {
        let standard = T::from_f64(293.15).unwrap_or_else(T::zero); // 20°C in Kelvin
        let tolerance = T::from_f64(0.1).unwrap_or_else(T::zero);
        (self.kelvin - standard).abs() < tolerance
    }

    /// Check if temperature is cryogenic (< 120K)
    pub fn is_cryogenic(&self) -> bool {
        let limit = T::from_f64(120.0).unwrap_or_else(T::zero);
        self.kelvin < limit
    }

    /// Check if temperature is high (> 1000K)
    pub fn is_high_temperature(&self) -> bool {
        let limit = T::from_f64(1000.0).unwrap_or_else(T::zero);
        self.kelvin > limit
    }
}

impl<T: RealField + Copy + fmt::Display> fmt::Display for Temperature<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} K", self.kelvin)
    }
}

/// Temperature units
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TemperatureUnit {
    /// Kelvin (K)
    Kelvin,
    /// Celsius (°C)
    Celsius,
    /// Fahrenheit (°F)
    Fahrenheit,
    /// Rankine (°R)
    Rankine,
}
