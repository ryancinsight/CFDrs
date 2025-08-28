//! Fluid properties and models.

use crate::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Fluid properties representation
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Fluid<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: Option<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: Option<T>,
}

impl<T: RealField + Copy> Fluid<T> {
    /// Create a fluid with required properties
    pub fn create(name: String, density: T, viscosity: T) -> Self {
        Self {
            name,
            density,
            viscosity,
            specific_heat: None,
            thermal_conductivity: None,
        }
    }

    /// Set specific heat
    pub fn with_specific_heat(mut self, cp: T) -> Self {
        self.specific_heat = Some(cp);
        self
    }

    /// Set thermal conductivity
    pub fn with_thermal_conductivity(mut self, k: T) -> Self {
        self.thermal_conductivity = Some(k);
        self
    }

    /// Get kinematic viscosity
    pub fn kinematic_viscosity(&self) -> T {
        self.viscosity / self.density
    }

    /// Get characteristic viscosity (same as dynamic viscosity for standard fluids)
    pub fn characteristic_viscosity(&self) -> T {
        self.viscosity
    }

    /// Get dynamic viscosity
    pub fn dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    /// Get Prandtl number (if thermal properties are set)
    pub fn prandtl_number(&self) -> Option<T> {
        match (self.specific_heat, self.thermal_conductivity) {
            (Some(cp), Some(k)) => Some(self.viscosity * cp / k),
            _ => None,
        }
    }

    /// Create water at 20°C
    pub fn water_20c() -> Self
    where
        T: FromPrimitive,
    {
        Self {
            name: "Water at 20°C".to_string(),
            density: T::from_f64(998.2).unwrap_or_else(T::zero),
            viscosity: T::from_f64(1.002e-3).unwrap_or_else(T::zero),
            specific_heat: Some(T::from_f64(4182.0).unwrap_or_else(T::zero)),
            thermal_conductivity: Some(T::from_f64(0.598).unwrap_or_else(T::zero)),
        }
    }

    /// Create air at 20°C
    pub fn air_20c() -> Self
    where
        T: FromPrimitive,
    {
        Self {
            name: "Air at 20°C".to_string(),
            density: T::from_f64(1.204).unwrap_or_else(T::zero),
            viscosity: T::from_f64(1.82e-5).unwrap_or_else(T::zero),
            specific_heat: Some(T::from_f64(1005.0).unwrap_or_else(T::zero)),
            thermal_conductivity: Some(T::from_f64(0.0257).unwrap_or_else(T::zero)),
        }
    }
}