//! Basic fluid properties and calculations

use crate::error::Error;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// A block of computed fluid properties at a single state point
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FluidProperties<T: RealField + Copy> {
    /// Density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub dynamic_viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
}

impl<T: RealField + Copy> FluidProperties<T> {
    /// Create new fluid properties
    pub fn new(
        density: T,
        dynamic_viscosity: T,
        specific_heat: T,
        thermal_conductivity: T,
    ) -> Self {
        Self {
            density,
            dynamic_viscosity,
            specific_heat,
            thermal_conductivity,
        }
    }

    /// Calculate kinematic viscosity [m²/s] from base properties
    ///
    /// # Errors
    /// Returns an error if density is non-positive
    pub fn kinematic_viscosity(&self) -> Result<T, Error> {
        if self.density <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        Ok(self.dynamic_viscosity / self.density)
    }

    /// Calculate Prandtl number from base properties
    ///
    /// # Errors
    /// Returns an error if thermal conductivity is non-positive
    pub fn prandtl_number(&self) -> Result<T, Error> {
        if self.thermal_conductivity <= T::zero() {
            return Err(Error::InvalidInput(
                "Thermal conductivity must be positive".to_string(),
            ));
        }
        Ok(self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity)
    }

    /// Calculate Reynolds number for given flow conditions
    ///
    /// # Errors
    /// Returns an error if the dynamic viscosity is zero or negative
    pub fn reynolds_number(&self, velocity: T, characteristic_length: T) -> Result<T, Error> {
        if self.dynamic_viscosity <= T::zero() {
            return Err(Error::InvalidInput(
                "Viscosity must be positive".to_string(),
            ));
        }
        Ok(self.density * velocity * characteristic_length / self.dynamic_viscosity)
    }

    /// Calculate Peclet number for given flow conditions
    ///
    /// # Errors
    /// Returns an error if viscosity, thermal conductivity, density, or specific heat are invalid
    pub fn peclet_number(&self, velocity: T, characteristic_length: T) -> Result<T, Error> {
        let re = self.reynolds_number(velocity, characteristic_length)?;
        let pr = self.prandtl_number()?;
        Ok(re * pr)
    }

    /// Calculate thermal diffusivity [m²/s]
    ///
    /// # Errors
    /// Returns an error if density or specific heat are non-positive
    pub fn thermal_diffusivity(&self) -> Result<T, Error> {
        if self.density <= T::zero() || self.specific_heat <= T::zero() {
            return Err(Error::InvalidInput(
                "Density and specific heat must be positive".to_string(),
            ));
        }
        Ok(self.thermal_conductivity / (self.density * self.specific_heat))
    }

    /// Calculate speed of sound for ideal gas [m/s]
    /// Requires ratio of specific heats (gamma) as parameter
    ///
    /// # Errors
    /// Returns an error if any input parameter is non-positive
    pub fn speed_of_sound(&self, gamma: T, temperature: T, gas_constant: T) -> Result<T, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        Ok((gamma * gas_constant * temperature).sqrt())
    }

    /// Calculate Mach number for given flow velocity
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidInput` if sound speed is non-positive.
    pub fn mach_number(&self, velocity: T, sound_speed: T) -> Result<T, Error> {
        if sound_speed <= T::zero() {
            return Err(Error::InvalidInput(
                "Sound speed must be positive".to_string(),
            ));
        }
        Ok(velocity / sound_speed)
    }
}
