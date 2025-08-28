//! Pressure value object

use crate::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::fmt;

// Pressure conversion constants
const PA_TO_BAR: f64 = 1e-5;
const PA_TO_PSI: f64 = 0.000145038;
const PA_TO_ATM: f64 = 9.86923e-6;
const PA_TO_MMHG: f64 = 0.00750062;

/// Pressure value with unit conversions
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Pressure<T: RealField + Copy> {
    /// Value in Pascals (SI unit)
    pascals: T,
}

impl<T: RealField + Copy + FromPrimitive> Pressure<T> {
    /// Create pressure in Pascals
    pub fn from_pascals(value: T) -> Result<Self> {
        Ok(Self { pascals: value })
    }

    /// Create pressure in bar
    pub fn from_bar(value: T) -> Result<Self> {
        let pa_to_bar = T::from_f64(PA_TO_BAR).unwrap_or_else(T::zero);
        if pa_to_bar > T::zero() {
            Ok(Self {
                pascals: value / pa_to_bar,
            })
        } else {
            Err(crate::error::Error::InvalidConfiguration(
                "Invalid conversion factor".into(),
            ))
        }
    }

    /// Create pressure in PSI
    pub fn from_psi(value: T) -> Result<Self> {
        let pa_to_psi = T::from_f64(PA_TO_PSI).unwrap_or_else(T::zero);
        if pa_to_psi > T::zero() {
            Ok(Self {
                pascals: value / pa_to_psi,
            })
        } else {
            Err(crate::error::Error::InvalidConfiguration(
                "Invalid conversion factor".into(),
            ))
        }
    }

    /// Get pressure in Pascals
    pub fn pascals(&self) -> T {
        self.pascals
    }

    /// Get pressure in bar
    pub fn bar(&self) -> T {
        self.pascals * T::from_f64(PA_TO_BAR).unwrap_or_else(T::zero)
    }

    /// Get pressure in PSI
    pub fn psi(&self) -> T {
        self.pascals * T::from_f64(PA_TO_PSI).unwrap_or_else(T::zero)
    }

    /// Get pressure in atmospheres
    pub fn atmospheres(&self) -> T {
        self.pascals * T::from_f64(PA_TO_ATM).unwrap_or_else(T::zero)
    }

    /// Get pressure in mmHg
    pub fn mmhg(&self) -> T {
        self.pascals * T::from_f64(PA_TO_MMHG).unwrap_or_else(T::zero)
    }

    /// Check if pressure is gauge (relative to atmospheric)
    pub fn is_gauge(&self) -> bool {
        self.pascals < T::from_f64(101325.0).unwrap_or_else(T::zero)
    }
}

impl<T: RealField + Copy + fmt::Display> fmt::Display for Pressure<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} Pa", self.pascals)
    }
}

/// Pressure units
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum PressureUnit {
    /// Pascals (Pa)
    Pascal,
    /// Bar
    Bar,
    /// Pounds per square inch (PSI)
    Psi,
    /// Atmospheres
    Atmosphere,
    /// Millimeters of mercury
    MmHg,
    /// Kilopascals
    KiloPascal,
    /// Megapascals
    MegaPascal,
}
