//! Casson fluid model.
//!
//! √μ = √μ_inf + √(τ_y / γ̇)

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Casson fluid model
///
/// Another common model for blood.
/// √μ = √μ_inf + √(τ_y / γ̇)
/// Or μ = (√μ_inf + √(τ_y / γ̇))²
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Casson<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Yield stress τ_y [Pa]
    pub yield_stress: T,
    /// Plastic viscosity μ_p (or μ_inf) [Pa·s]
    pub plastic_viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + FromPrimitive + Copy> Casson<T> {
    /// Create a new Casson fluid
    pub fn new(
        name: String,
        density: T,
        yield_stress: T,
        plastic_viscosity: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
    ) -> Self {
        Self {
            name,
            density,
            yield_stress,
            plastic_viscosity,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
        }
    }

    /// Calculate apparent viscosity at given shear rate
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= T::from_f64(1e-6).unwrap_or(T::zero()) {
            return T::from_f64(100.0).unwrap_or_else(T::one);
        }

        let sqrt_mu = self.plastic_viscosity.sqrt() + (self.yield_stress / shear_rate).sqrt();
        sqrt_mu * sqrt_mu
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for Casson<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: self.plastic_viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn viscosity_at_shear(&self, shear_rate: T, _temperature: T, _pressure: T) -> Result<T, Error> {
        Ok(self.apparent_viscosity(shear_rate))
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for Casson<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        Casson::apparent_viscosity(self, shear_rate)
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<T> {
        Some(self.yield_stress)
    }
}
