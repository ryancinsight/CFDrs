//! Temperature-dependent fluid models

use super::traits::{Fluid as FluidTrait, FluidState};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Temperature-dependent viscosity model using polynomial fit
///
/// μ(T) = Σ(a_i * T^i) for i = 0 to n
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolynomialViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Base density [kg/m³]
    pub density_ref: T,
    /// Temperature coefficient for density [1/K]
    pub thermal_expansion: T,
    /// Polynomial coefficients for viscosity [Pa·s/K^i]
    pub viscosity_coeffs: Vec<T>,
    /// Reference temperature [K]
    pub t_ref: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
}

impl<T: RealField + FromPrimitive + Copy> PolynomialViscosity<T> {
    /// Calculate viscosity at given temperature
    pub fn calculate_viscosity(&self, temperature: T) -> T {
        let mut viscosity = T::zero();
        let mut t_power = T::one();

        for coeff in &self.viscosity_coeffs {
            viscosity = viscosity + *coeff * t_power;
            t_power = t_power * temperature;
        }

        viscosity
    }

    /// Calculate density with thermal expansion
    pub fn calculate_density(&self, temperature: T) -> T {
        let dt = temperature - self.t_ref;
        self.density_ref * (T::one() - self.thermal_expansion * dt)
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for PolynomialViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        let density = self.calculate_density(temperature);
        let viscosity = self.calculate_viscosity(temperature);

        Ok(FluidState {
            density,
            dynamic_viscosity: viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }

    fn reference_temperature(&self) -> Option<T> {
        Some(self.t_ref)
    }
}

/// Arrhenius viscosity model for liquids
///
/// μ(T) = A * exp(B/T)
/// Common for liquids where viscosity decreases with temperature
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArrheniusViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Pre-exponential factor A [Pa·s]
    pub a_factor: T,
    /// Activation energy parameter B [K]
    pub b_factor: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
}

impl<T: RealField + FromPrimitive + Copy> ArrheniusViscosity<T> {
    /// Calculate viscosity using Arrhenius model
    pub fn calculate_viscosity(&self, temperature: T) -> T {
        self.a_factor * (self.b_factor / temperature).exp()
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for ArrheniusViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        let viscosity = self.calculate_viscosity(temperature);

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }
}

/// Andrade viscosity model
///
/// μ(T) = A * exp(B/(T - C))
/// Extended version of Arrhenius with additional parameter
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AndradeViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Pre-exponential factor A [Pa·s]
    pub a_factor: T,
    /// Temperature coefficient B [K]
    pub b_factor: T,
    /// Temperature offset C [K]
    pub c_factor: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
}

impl<T: RealField + FromPrimitive + Copy> AndradeViscosity<T> {
    /// Calculate viscosity using Andrade model
    pub fn calculate_viscosity(&self, temperature: T) -> Result<T, Error> {
        let denominator = temperature - self.c_factor;
        if denominator <= T::zero() {
            return Err(Error::InvalidInput(format!(
                "Temperature must be greater than {}",
                self.c_factor
            )));
        }

        Ok(self.a_factor * (self.b_factor / denominator).exp())
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for AndradeViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let viscosity = self.calculate_viscosity(temperature)?;

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }
}
