//! Herschel–Bulkley fluid model (generalised Bingham plastic).
//!
//! τ = τ₀ + K · (γ̇)^n  for τ > τ₀

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use crate::physics::constants::physics::universal::GAS_CONSTANT;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

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
    /// Reference temperature for properties [K]
    #[serde(default)]
    pub reference_temperature: Option<T>,
    /// Activation energy for consistency index [J/mol]
    #[serde(default)]
    pub activation_energy_k: Option<T>,
    /// Activation energy for yield stress [J/mol]
    #[serde(default)]
    pub activation_energy_tau: Option<T>,
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
            reference_temperature: None,
            activation_energy_k: None,
            activation_energy_tau: None,
        }
    }

    /// Set temperature dependence parameters
    pub fn with_temperature_dependence(
        mut self,
        reference_temperature: T,
        activation_energy_k: Option<T>,
        activation_energy_tau: Option<T>,
    ) -> Self {
        self.reference_temperature = Some(reference_temperature);
        self.activation_energy_k = activation_energy_k;
        self.activation_energy_tau = activation_energy_tau;
        self
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
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let (k_adj, tau_adj) = if let Some(t_ref) = self.reference_temperature {
            let r_gas = T::from_f64(GAS_CONSTANT).unwrap_or_else(T::one);
            let inv_t = T::one() / temperature;
            let inv_t_ref = T::one() / t_ref;
            let diff_inv_t = inv_t - inv_t_ref;

            let k = if let Some(ea_k) = self.activation_energy_k {
                let arg = (ea_k / r_gas) * diff_inv_t;
                self.consistency_index * arg.exp()
            } else {
                self.consistency_index
            };

            let tau = if let Some(ea_tau) = self.activation_energy_tau {
                let arg = (ea_tau / r_gas) * diff_inv_t;
                self.yield_stress * arg.exp()
            } else {
                self.yield_stress
            };

            (k, tau)
        } else {
            (self.consistency_index, self.yield_stress)
        };

        let shear_rate = self.reference_shear_rate;
        let apparent_viscosity = if shear_rate <= T::zero() {
            T::from_f64(1e6).unwrap_or_else(T::one)
        } else {
            let power_law_term = k_adj * shear_rate.powf(self.flow_behavior_index - T::one());
            tau_adj / shear_rate + power_law_term
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
