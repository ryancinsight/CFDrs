//! Bingham plastic fluid model.
//!
//! τ = τ₀ + μ_p · γ̇  for τ > τ₀

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, Pressure, ReciprocalTime, SpecificHeatCapacity,
    ThermalConductivity, Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Bingham plastic fluid model
///
/// τ = τ₀ + `μ_p` * γ̇  for τ > τ₀
/// γ̇ = 0            for τ ≤ τ₀
/// where τ₀ is yield stress, `μ_p` is plastic viscosity
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinghamPlastic<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: MassDensity<T>,
    /// Yield stress τ₀ \[Pa]
    pub yield_stress: Pressure<T>,
    /// Plastic viscosity `μ_p` [Pa·s]
    pub plastic_viscosity: DynamicViscosity<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
    /// Reference shear rate for viscosity calculation [1/s]
    pub reference_shear_rate: ReciprocalTime<T>,
}

impl<T: RealField + FloatElement + Copy> BinghamPlastic<T> {
    /// Create a new Bingham plastic fluid
    pub fn new(
        name: String,
        density: MassDensity<T>,
        yield_stress: Pressure<T>,
        plastic_viscosity: DynamicViscosity<T>,
        specific_heat: SpecificHeatCapacity<T>,
        thermal_conductivity: ThermalConductivity<T>,
        speed_of_sound: Velocity<T>,
        reference_shear_rate: ReciprocalTime<T>,
    ) -> Self {
        Self {
            name,
            density,
            yield_stress,
            plastic_viscosity,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
            reference_shear_rate,
        }
    }

    /// Calculate apparent viscosity at given shear rate
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= <T as NumericElement>::ZERO {
            return <T as FloatElement>::from_f64(1e6);
        }

        self.plastic_viscosity.into_base() + self.yield_stress.into_base() / shear_rate
    }

    /// Check if fluid is yielded at given shear stress
    pub fn is_yielded(&self, shear_stress: T) -> bool {
        shear_stress > self.yield_stress.into_base()
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for BinghamPlastic<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let apparent_viscosity = self.apparent_viscosity(self.reference_shear_rate.into_base());

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: DynamicViscosity::from_base(apparent_viscosity),
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

impl<T: RealField + FloatElement + Copy> NonNewtonianFluid<T> for BinghamPlastic<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(BinghamPlastic::apparent_viscosity(self, shear_rate))
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<Pressure<T>> {
        Some(self.yield_stress)
    }
}
