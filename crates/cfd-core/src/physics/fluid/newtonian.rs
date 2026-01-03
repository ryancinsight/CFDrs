//! Newtonian fluid models with constant and variable properties

use super::traits::{ConstantFluid, Fluid as FluidTrait, FluidState};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Constant property fluid model (incompressible, Newtonian)
///
/// This model assumes fluid properties are independent of temperature and pressure.
/// Suitable for isothermal, incompressible flow simulations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConstantPropertyFluid<T: RealField + Copy> {
    /// Descriptive name of the fluid
    pub name: String,
    /// Constant density [kg/m³]
    pub density: T,
    /// Constant dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Constant specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Constant thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Constant speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + Copy> ConstantPropertyFluid<T> {
    /// Create a new constant property fluid
    pub fn new(
        name: String,
        density: T,
        viscosity: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
    ) -> Self {
        Self {
            name,
            density,
            viscosity,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
        }
    }

    /// Validate that all properties are physically reasonable
    ///
    /// # Errors
    /// Returns an error if any physical property is non-positive or invalid
    pub fn validate(&self) -> Result<(), Error> {
        if self.density <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        if self.viscosity <= T::zero() {
            return Err(Error::InvalidInput(
                "Viscosity must be positive".to_string(),
            ));
        }
        if self.specific_heat <= T::zero() {
            return Err(Error::InvalidInput(
                "Specific heat must be positive".to_string(),
            ));
        }
        if self.thermal_conductivity <= T::zero() {
            return Err(Error::InvalidInput(
                "Thermal conductivity must be positive".to_string(),
            ));
        }
        if self.speed_of_sound <= T::zero() {
            return Err(Error::InvalidInput(
                "Speed of sound must be positive".to_string(),
            ));
        }
        Ok(())
    }

    /// Water at 20°C and standard pressure
    /// Properties from NIST webbook
    ///
    /// # Errors
    /// Returns an error if numeric conversion from f64 fails for the target type T
    pub fn water_20c() -> Result<Self, Error>
    where
        T: FromPrimitive,
    {
        let fluid = Self::new(
            "Water (20°C)".to_string(),
            T::from_f64(998.2).unwrap(),     // density [kg/m³]
            T::from_f64(0.001_002).unwrap(), // viscosity [Pa·s]
            T::from_f64(4186.0).unwrap(),    // specific heat [J/(kg·K)]
            T::from_f64(0.599).unwrap(),     // thermal conductivity [W/(m·K)]
            T::from_f64(1482.0).unwrap(),    // speed of sound [m/s]
        );
        fluid.validate()?;
        Ok(fluid)
    }

    /// Air at 20°C and standard pressure
    /// Properties from NIST webbook
    ///
    /// # Errors
    /// Returns an error if numeric conversion from f64 fails for the target type T
    pub fn air_20c() -> Result<Self, Error>
    where
        T: FromPrimitive,
    {
        let fluid = Self::new(
            "Air (20°C)".to_string(),
            T::from_f64(1.204).unwrap(),    // density [kg/m³]
            T::from_f64(1.825e-5).unwrap(), // viscosity [Pa·s]
            T::from_f64(1005.0).unwrap(),   // specific heat [J/(kg·K)]
            T::from_f64(0.02538).unwrap(),  // thermal conductivity [W/(m·K)]
            T::from_f64(343.2).unwrap(),    // speed of sound [m/s]
        );
        fluid.validate()?;
        Ok(fluid)
    }
}

impl<T: RealField + Copy> FluidTrait<T> for ConstantPropertyFluid<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: self.viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

impl<T: RealField + Copy> ConstantFluid<T> for ConstantPropertyFluid<T> {
    fn density(&self) -> T {
        self.density
    }

    fn dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn speed_of_sound(&self) -> T {
        self.speed_of_sound
    }
}

/// Ideal gas model with constant specific heats
///
/// Properties vary with temperature and pressure according to ideal gas law
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IdealGas<T: RealField + Copy> {
    /// Gas name
    pub name: String,
    /// Specific gas constant [J/(kg·K)]
    pub gas_constant: T,
    /// Specific heat at constant pressure [J/(kg·K)]
    pub cp: T,
    /// Reference viscosity [Pa·s] at reference temperature
    pub mu_ref: T,
    /// Reference temperature [K] for viscosity
    pub t_ref: T,
    /// Sutherland constant [K] for viscosity model
    pub sutherland_constant: T,
    /// Thermal conductivity coefficient
    pub k_coeff: T,
}

impl<T: RealField + FromPrimitive + Copy> IdealGas<T> {
    /// Create a new ideal gas
    pub fn new(
        name: String,
        gas_constant: T,
        cp: T,
        mu_ref: T,
        t_ref: T,
        sutherland_constant: T,
    ) -> Self {
        // Thermal conductivity from Prandtl number assumption (Pr ≈ 0.7 for air)
        let pr = T::from_f64(0.7).unwrap_or_else(T::one);
        let k_coeff = cp / pr;

        Self {
            name,
            gas_constant,
            cp,
            mu_ref,
            t_ref,
            sutherland_constant,
            k_coeff,
        }
    }

    /// Calculate density from ideal gas law
    fn calculate_density(&self, temperature: T, pressure: T) -> T {
        pressure / (self.gas_constant * temperature)
    }

    /// Calculate viscosity using Sutherland's law
    fn calculate_viscosity(&self, temperature: T) -> T {
        let t_ratio = temperature / self.t_ref;
        let numerator = self.t_ref + self.sutherland_constant;
        let denominator = temperature + self.sutherland_constant;

        self.mu_ref * t_ratio.powf(T::from_f64(1.5).unwrap_or_else(T::one)) * numerator
            / denominator
    }

    /// Calculate thermal conductivity using kinetic theory
    fn calculate_thermal_conductivity(&self, viscosity: T) -> T {
        viscosity * self.k_coeff
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for IdealGas<T> {
    fn properties_at(&self, temperature: T, pressure: T) -> Result<FluidState<T>, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        if pressure <= T::zero() {
            return Err(Error::InvalidInput("Pressure must be positive".to_string()));
        }

        let density = self.calculate_density(temperature, pressure);
        let viscosity = self.calculate_viscosity(temperature);
        let thermal_conductivity = self.calculate_thermal_conductivity(viscosity);

        // Speed of sound for ideal gas: c = sqrt(gamma * R * T)
        // gamma = cp / cv = cp / (cp - R)
        let cv = self.cp - self.gas_constant;
        let gamma = self.cp / cv;
        let speed_of_sound = (gamma * self.gas_constant * temperature).sqrt();

        Ok(FluidState {
            density,
            dynamic_viscosity: viscosity,
            specific_heat: self.cp,
            thermal_conductivity,
            speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }

    fn is_pressure_dependent(&self) -> bool {
        true
    }
}
