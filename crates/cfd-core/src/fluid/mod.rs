//! Fluid properties and models with trait-based extensibility
//!
//! This module provides a comprehensive framework for fluid modeling in CFD simulations,
//! supporting Newtonian, non-Newtonian, and temperature-dependent fluid behaviors.

use crate::error::Error;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

pub mod database;
pub mod newtonian;
pub mod non_newtonian;
pub mod properties;
pub mod temperature;
pub mod validation;

// Re-export core types
pub use database::{air_20c, water_20c};
pub use newtonian::ConstantPropertyFluid;
pub use properties::FluidProperties;

/// Type alias for simplified naming when using constant property fluids
pub type Fluid<T> = ConstantPropertyFluid<T>;

/// Trait defining the interface for fluid models
///
/// This trait allows for different fluid models to be implemented
/// while providing a consistent interface to the solvers.
pub trait FluidModel<T: RealField + Copy>: Send + Sync {
    /// Calculate a block of all fluid properties at the given conditions
    fn properties(&self, temperature: T, pressure: T) -> Result<FluidProperties<T>, Error>;

    /// Get a descriptive name for this fluid model
    fn name(&self) -> &str;

    /// Get fluid density [kg/m³] at given conditions
    fn density(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.density)
    }

    /// Get dynamic viscosity [Pa·s] at given conditions
    fn dynamic_viscosity(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.dynamic_viscosity)
    }

    /// Get kinematic viscosity [m²/s] at given conditions
    fn kinematic_viscosity(&self, temperature: T, pressure: T) -> Result<T, Error> {
        self.properties(temperature, pressure)?
            .kinematic_viscosity()
    }

    /// Get specific heat capacity [J/(kg·K)] at given conditions
    fn specific_heat(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.specific_heat)
    }

    /// Get thermal conductivity [W/(m·K)] at given conditions
    fn thermal_conductivity(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.thermal_conductivity)
    }

    /// Get Prandtl number (dimensionless) at given conditions
    fn prandtl_number(&self, temperature: T, pressure: T) -> Result<T, Error> {
        self.properties(temperature, pressure)?.prandtl_number()
    }

    /// Calculate Reynolds number for given flow conditions
    fn reynolds_number(
        &self,
        velocity: T,
        characteristic_length: T,
        temperature: T,
        pressure: T,
    ) -> Result<T, Error> {
        self.properties(temperature, pressure)?
            .reynolds_number(velocity, characteristic_length)
    }
}
