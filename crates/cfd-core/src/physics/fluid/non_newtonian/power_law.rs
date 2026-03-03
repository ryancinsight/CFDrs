//! Power-law fluid model (Ostwald–de Waele model).
//!
//! τ = K · (γ̇)^n

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use crate::physics::constants::physics::universal::GAS_CONSTANT;
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
    /// Reference temperature for properties [K]
    #[serde(default)]
    pub reference_temperature: Option<T>,
    /// Activation energy for consistency index [J/mol]
    #[serde(default)]
    pub activation_energy_k: Option<T>,
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
            reference_temperature: None,
            activation_energy_k: None,
        }
    }

    /// Set temperature dependence parameters
    pub fn with_temperature_dependence(
        mut self,
        reference_temperature: T,
        activation_energy_k: Option<T>,
    ) -> Self {
        self.reference_temperature = Some(reference_temperature);
        self.activation_energy_k = activation_energy_k;
        self
    }

    /// Calculate apparent viscosity at given shear rate
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= T::zero() {
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
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let k_adj = if let Some(t_ref) = self.reference_temperature {
            let r_gas = T::from_f64(GAS_CONSTANT).unwrap_or_else(T::one);
            let inv_t = T::one() / temperature;
            let inv_t_ref = T::one() / t_ref;
            let diff_inv_t = inv_t - inv_t_ref;

            if let Some(ea_k) = self.activation_energy_k {
                let arg = (ea_k / r_gas) * diff_inv_t;
                self.consistency_index * arg.exp()
            } else {
                self.consistency_index
            }
        } else {
            self.consistency_index
        };

        let shear_rate = self.reference_shear_rate;
        let apparent_viscosity = if shear_rate <= T::zero() {
            k_adj * T::from_f64(1000.0).unwrap_or_else(T::one)
        } else {
            k_adj * shear_rate.powf(self.flow_behavior_index - T::one())
        };

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

    fn is_temperature_dependent(&self) -> bool {
        self.reference_temperature.is_some()
    }

    fn reference_temperature(&self) -> Option<T> {
        self.reference_temperature
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for PowerLawFluid<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        PowerLawFluid::apparent_viscosity(self, shear_rate)
    }
}
