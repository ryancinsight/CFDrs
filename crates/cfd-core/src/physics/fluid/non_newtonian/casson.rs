//! Casson fluid model.
//!
//! √μ = √μ_inf + √(τ_y / γ̇)

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, Pressure, SpecificHeatCapacity, ThermalConductivity, Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
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
    pub density: MassDensity<T>,
    /// Yield stress τ_y \[Pa]
    pub yield_stress: Pressure<T>,
    /// Plastic viscosity μ_p (or μ_inf) [Pa·s]
    pub plastic_viscosity: DynamicViscosity<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
}

impl<T: RealField + FloatElement + Copy> Casson<T> {
    /// Create a new Casson fluid
    pub fn new(
        name: String,
        density: MassDensity<T>,
        yield_stress: Pressure<T>,
        plastic_viscosity: DynamicViscosity<T>,
        specific_heat: SpecificHeatCapacity<T>,
        thermal_conductivity: ThermalConductivity<T>,
        speed_of_sound: Velocity<T>,
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
        if shear_rate <= <T as FloatElement>::from_f64(1e-6) {
            return <T as FloatElement>::from_f64(100.0);
        }

        let sqrt_mu = <T as NumericElement>::sqrt(self.plastic_viscosity.into_base())
            + <T as NumericElement>::sqrt(self.yield_stress.into_base() / shear_rate);
        sqrt_mu * sqrt_mu
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for Casson<T> {
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

    fn viscosity_at_shear(
        &self,
        shear_rate: T,
        _temperature: T,
        _pressure: T,
    ) -> Result<DynamicViscosity<T>, Error> {
        Ok(DynamicViscosity::from_base(
            self.apparent_viscosity(shear_rate),
        ))
    }
}

impl<T: RealField + FloatElement + Copy> NonNewtonianFluid<T> for Casson<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(Casson::apparent_viscosity(self, shear_rate))
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<Pressure<T>> {
        Some(self.yield_stress)
    }
}
