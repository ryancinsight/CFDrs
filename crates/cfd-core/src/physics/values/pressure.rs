//! Pressure value object

use crate::error::Result;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};
use std::fmt;

// Pressure conversion constants
const PA_TO_BAR: f64 = 1e-5;
const PA_TO_PSI: f64 = 0.000_145_038;
const PA_TO_ATM: f64 = 9.86923e-6;
const PA_TO_MMHG: f64 = 0.007_500_62;

/// Pressure value with unit conversions
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Pressure<T: FloatElement + Copy> {
    /// Value in Pascals (SI unit)
    pascals: T,
}

impl<T: FloatElement + Copy> Pressure<T> {
    /// Create zero pressure
    #[must_use]
    pub fn zero() -> Self {
        Self {
            pascals: <T as NumericElement>::ZERO,
        }
    }

    /// Create pressure in Pascals
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn from_pascals(value: T) -> Result<Self> {
        Ok(Self { pascals: value })
    }

    /// Create pressure in bar
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn from_bar(value: T) -> Result<Self> {
        let pa_to_bar = <T as FloatElement>::from_f64(PA_TO_BAR);
        if pa_to_bar > <T as NumericElement>::ZERO {
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
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn from_psi(value: T) -> Result<Self> {
        let pa_to_psi = <T as FloatElement>::from_f64(PA_TO_PSI);
        if pa_to_psi > <T as NumericElement>::ZERO {
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
        self.pascals * <T as FloatElement>::from_f64(PA_TO_BAR)
    }

    /// Get pressure in PSI
    pub fn psi(&self) -> T {
        self.pascals * <T as FloatElement>::from_f64(PA_TO_PSI)
    }

    /// Get pressure in atmospheres
    pub fn atmospheres(&self) -> T {
        self.pascals * <T as FloatElement>::from_f64(PA_TO_ATM)
    }

    /// Get pressure in mmHg
    pub fn mmhg(&self) -> T {
        self.pascals * <T as FloatElement>::from_f64(PA_TO_MMHG)
    }

    /// Check if pressure is gauge (relative to atmospheric)
    pub fn is_gauge(&self) -> bool {
        use crate::physics::constants::physics::thermo::P_ATM;
        self.pascals < <T as FloatElement>::from_f64(P_ATM)
    }
}

impl<T: FloatElement + Copy + fmt::Display> fmt::Display for Pressure<T> {
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
