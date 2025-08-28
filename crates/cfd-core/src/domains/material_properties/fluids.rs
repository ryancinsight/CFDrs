//! Newtonian fluid implementations

use super::traits::FluidProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Newtonian fluid with constant properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NewtonianFluid<T: RealField + Copy> {
    /// Fluid density
    pub density: T,
    /// Dynamic viscosity
    pub dynamic_viscosity: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
}

impl<T: RealField + Copy> FluidProperties<T> for NewtonianFluid<T> {
    fn density(&self) -> T {
        self.density
    }

    fn dynamic_viscosity(&self) -> T {
        self.dynamic_viscosity
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }
}

impl<T: RealField + Copy> NewtonianFluid<T> {
    /// Create water properties at 20°C
    pub fn water() -> Self
    where
        T: From<f64>,
    {
        Self {
            density: T::from(998.2),
            dynamic_viscosity: T::from(1.002e-3),
            thermal_conductivity: T::from(0.598),
            specific_heat: T::from(4182.0),
        }
    }

    /// Create air properties at 20°C, 1 atm
    pub fn air() -> Self
    where
        T: From<f64>,
    {
        Self {
            density: T::from(1.204),
            dynamic_viscosity: T::from(1.825e-5),
            thermal_conductivity: T::from(0.0257),
            specific_heat: T::from(1005.0),
        }
    }

    /// Create oil properties (typical engine oil)
    pub fn oil() -> Self
    where
        T: From<f64>,
    {
        Self {
            density: T::from(890.0),
            dynamic_viscosity: T::from(0.1),
            thermal_conductivity: T::from(0.15),
            specific_heat: T::from(2000.0),
        }
    }
}
