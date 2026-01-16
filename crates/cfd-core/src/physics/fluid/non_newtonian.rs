//! Non-Newtonian fluid models
//!
//! References:
//! - Bird, R.B., Armstrong, R.C., Hassager, O. (1987) "Dynamics of Polymeric Liquids"
//! - Chhabra, R.P., Richardson, J.F. (2008) "Non-Newtonian Flow and Applied Rheology"

use super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Power-law fluid model (Ostwald-de Waele model)
///
/// τ = K * (γ̇)^n
/// where τ is shear stress, γ̇ is shear rate, K is consistency index, n is flow behavior index
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PowerLawFluid<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Consistency index K [Pa·s^n]
    pub consistency_index: T,
    /// Flow behavior index n [-]
    /// n < 1: shear-thinning (pseudoplastic)
    /// n = 1: Newtonian
    /// n > 1: shear-thickening (dilatant)
    pub flow_behavior_index: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
    /// Reference shear rate for viscosity calculation [1/s]
    pub reference_shear_rate: T,
}

impl<T: RealField + FromPrimitive + Copy> PowerLawFluid<T> {
    /// Create a new power-law fluid
    pub fn new(
        name: String,
        density: T,
        consistency_index: T,
        flow_behavior_index: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
        reference_shear_rate: T,
    ) -> Self {
        Self {
            name,
            density,
            consistency_index,
            flow_behavior_index,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
            reference_shear_rate,
        }
    }

    /// Calculate apparent viscosity at given shear rate
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= T::zero() {
            // Return high viscosity for zero shear rate
            return self.consistency_index * T::from_f64(1000.0).unwrap_or_else(T::one);
        }

        self.consistency_index * shear_rate.powf(self.flow_behavior_index - T::one())
    }

    /// Validate model parameters
    ///
    /// # Errors
    /// Returns an error if any parameter is non-positive or physically invalid
    pub fn validate(&self) -> Result<(), Error> {
        if self.density <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        if self.consistency_index <= T::zero() {
            return Err(Error::InvalidInput(
                "Consistency index must be positive".to_string(),
            ));
        }
        if self.flow_behavior_index <= T::zero() {
            return Err(Error::InvalidInput(
                "Flow behavior index must be positive".to_string(),
            ));
        }
        Ok(())
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for PowerLawFluid<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        // TODO: Use reference shear rate for property calculation
        // DEPENDENCIES: Implement temperature-dependent viscosity models for non-Newtonian fluids
        // BLOCKED BY: Limited understanding of shear rate effects on temperature-dependent properties
        // PRIORITY: High - Essential for accurate non-Newtonian CFD simulations
        let apparent_viscosity = self.apparent_viscosity(self.reference_shear_rate);

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: apparent_viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for PowerLawFluid<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        PowerLawFluid::apparent_viscosity(self, shear_rate)
    }
}

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
    pub density: T,
    /// Yield stress τ₀ [Pa]
    pub yield_stress: T,
    /// Plastic viscosity `μ_p` [Pa·s]
    pub plastic_viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
    /// Reference shear rate for viscosity calculation [1/s]
    pub reference_shear_rate: T,
}

impl<T: RealField + FromPrimitive + Copy> BinghamPlastic<T> {
    /// Create a new Bingham plastic fluid
    pub fn new(
        name: String,
        density: T,
        yield_stress: T,
        plastic_viscosity: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
        reference_shear_rate: T,
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
        if shear_rate <= T::zero() {
            // Return high viscosity for zero shear rate
            return T::from_f64(1e6).unwrap_or_else(T::one);
        }

        self.plastic_viscosity + self.yield_stress / shear_rate
    }

    /// Check if fluid is yielded at given shear stress
    pub fn is_yielded(&self, shear_stress: T) -> bool {
        shear_stress > self.yield_stress
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for BinghamPlastic<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        // TODO: Use reference shear rate for property calculation
        // DEPENDENCIES: Implement temperature-dependent yield stress and viscosity models
        // BLOCKED BY: Limited understanding of temperature effects on yield stress behavior
        // PRIORITY: High - Essential for accurate Bingham plastic CFD simulations
        let apparent_viscosity = self.apparent_viscosity(self.reference_shear_rate);

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: apparent_viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for BinghamPlastic<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        BinghamPlastic::apparent_viscosity(self, shear_rate)
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<T> {
        Some(self.yield_stress)
    }
}

/// Herschel-Bulkley fluid model (generalized Bingham plastic)
///
/// τ = τ₀ + K * (γ̇)^n  for τ > τ₀
/// γ̇ = 0               for τ ≤ τ₀
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HerschelBulkley<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Yield stress τ₀ [Pa]
    pub yield_stress: T,
    /// Consistency index K [Pa·s^n]
    pub consistency_index: T,
    /// Flow behavior index n [-]
    pub flow_behavior_index: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
    /// Reference shear rate [1/s]
    pub reference_shear_rate: T,
}

impl<T: RealField + FromPrimitive + Copy> HerschelBulkley<T> {
    /// Create a new Herschel-Bulkley fluid
    pub fn new(
        name: String,
        density: T,
        yield_stress: T,
        consistency_index: T,
        flow_behavior_index: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
        reference_shear_rate: T,
    ) -> Self {
        Self {
            name,
            density,
            yield_stress,
            consistency_index,
            flow_behavior_index,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
            reference_shear_rate,
        }
    }
    /// Calculate apparent viscosity at given shear rate
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= T::zero() {
            return T::from_f64(1e6).unwrap_or_else(T::one);
        }

        let power_law_term =
            self.consistency_index * shear_rate.powf(self.flow_behavior_index - T::one());
        self.yield_stress / shear_rate + power_law_term
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for HerschelBulkley<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        // TODO: Implement temperature-dependent Herschel-Bulkley properties
        // DEPENDENCIES: Add temperature effects to yield stress, consistency index, and flow behavior
        // BLOCKED BY: Complex coupling between temperature and non-Newtonian parameters
        // PRIORITY: High - Essential for accurate Herschel-Bulkley CFD simulations
        let apparent_viscosity = self.apparent_viscosity(self.reference_shear_rate);

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: apparent_viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for HerschelBulkley<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        HerschelBulkley::apparent_viscosity(self, shear_rate)
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<T> {
        Some(self.yield_stress)
    }
}
