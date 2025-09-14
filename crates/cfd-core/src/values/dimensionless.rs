//! Dimensionless numbers for CFD

use crate::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Type of dimensionless number
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum DimensionlessType {
    /// Reynolds number (Re = ρvL/μ)
    Reynolds,
    /// Prandtl number (Pr = μCp/k)
    Prandtl,
    /// Nusselt number (Nu = hL/k)
    Nusselt,
    /// Peclet number (Pe = Re*Pr)
    Peclet,
    /// Froude number (Fr = v/√(gL))
    Froude,
    /// Weber number (We = ρv²L/σ)
    Weber,
    /// Mach number (Ma = v/a)
    Mach,
    /// Strouhal number (St = fL/v)
    Strouhal,
    /// Grashof number (Gr = gβΔTL³/ν²)
    Grashof,
    /// Rayleigh number (Ra = Gr*Pr)
    Rayleigh,
    /// Schmidt number (Sc = ν/D)
    Schmidt,
    /// Lewis number (Le = α/D)
    Lewis,
    /// Damköhler number (Da = reaction rate / flow rate)
    Damkohler,
    /// Biot number (Bi = hL/k)
    Biot,
    /// Euler number (Eu = Δp/(ρv²))
    Euler,
}

/// Generic dimensionless number
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct DimensionlessNumber<T: RealField + Copy> {
    value: T,
    number_type: DimensionlessType,
}

impl<T: RealField + Copy + FromPrimitive> DimensionlessNumber<T> {
    /// Create a new dimensionless number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn new(value: T, number_type: DimensionlessType) -> Result<Self> {
        // Some dimensionless numbers must be non-negative
        match number_type {
            DimensionlessType::Reynolds
            | DimensionlessType::Prandtl
            | DimensionlessType::Schmidt
            | DimensionlessType::Lewis => {
                if value < T::zero() {
                    return Err(crate::error::Error::InvalidConfiguration(format!(
                        "{:?} must be non-negative",
                        number_type
                    )));
                }
            }
            _ => {}
        }
        Ok(Self { value, number_type })
    }

    /// Get the value
    pub fn value(&self) -> T {
        self.value
    }

    /// Get the type
    pub fn number_type(&self) -> DimensionlessType {
        self.number_type
    }

    /// Create Reynolds number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn reynolds(density: T, velocity: T, length: T, viscosity: T) -> Result<Self> {
        if viscosity <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Viscosity must be positive".into(),
            ));
        }
        let re = density * velocity * length / viscosity;
        Self::new(re, DimensionlessType::Reynolds)
    }

    /// Create Prandtl number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn prandtl(viscosity: T, specific_heat: T, thermal_conductivity: T) -> Result<Self> {
        if thermal_conductivity <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Thermal conductivity must be positive".into(),
            ));
        }
        let pr = viscosity * specific_heat / thermal_conductivity;
        Self::new(pr, DimensionlessType::Prandtl)
    }

    /// Create Nusselt number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn nusselt(heat_transfer_coeff: T, length: T, thermal_conductivity: T) -> Result<Self> {
        if thermal_conductivity <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Thermal conductivity must be positive".into(),
            ));
        }
        let nu = heat_transfer_coeff * length / thermal_conductivity;
        Self::new(nu, DimensionlessType::Nusselt)
    }

    /// Create Mach number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn mach(velocity: T, speed_of_sound: T) -> Result<Self> {
        if speed_of_sound <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Speed of sound must be positive".into(),
            ));
        }
        let ma = velocity / speed_of_sound;
        Self::new(ma, DimensionlessType::Mach)
    }

    /// Create Froude number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn froude(velocity: T, gravity: T, length: T) -> Result<Self> {
        if gravity <= T::zero() || length <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Gravity and length must be positive".into(),
            ));
        }
        let fr = velocity / (gravity * length).sqrt();
        Self::new(fr, DimensionlessType::Froude)
    }

    /// Create Weber number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn weber(density: T, velocity: T, length: T, surface_tension: T) -> Result<Self> {
        if surface_tension <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Surface tension must be positive".into(),
            ));
        }
        let we = density * velocity * velocity * length / surface_tension;
        Self::new(we, DimensionlessType::Weber)
    }
}

impl<T: RealField + Copy + fmt::Display> fmt::Display for DimensionlessNumber<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} = {}", self.number_type, self.value)
    }
}
