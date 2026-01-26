//! Non-Newtonian fluid models
//!
//! References:
//! - Bird, R.B., Armstrong, R.C., Hassager, O. (1987) "Dynamics of Polymeric Liquids"
//! - Chhabra, R.P., Richardson, J.F. (2008) "Non-Newtonian Flow and Applied Rheology"

use super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
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
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        // Calculate adjusted consistency index based on temperature if reference temperature is set
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

        // Calculate apparent viscosity using adjusted properties
        let shear_rate = self.reference_shear_rate;
        let apparent_viscosity = if shear_rate <= T::zero() {
            // Return high viscosity for zero shear rate
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
        // Calculate adjusted properties based on temperature if reference temperature is set
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

        // Calculate apparent viscosity using adjusted properties
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_herschel_bulkley_constant() {
        let fluid = HerschelBulkley::<f64>::new(
            "Test Fluid".to_string(),
            1000.0,
            10.0, // yield stress
            5.0,  // consistency index
            0.5,  // flow behavior index
            4000.0,
            0.6,
            1500.0,
            10.0, // reference shear rate
        );

        let props = fluid.properties_at(300.0, 101325.0).unwrap();

        // Manual calc:
        // shear_rate = 10.0
        // power_law_term = 5.0 * 10.0^(0.5 - 1.0) = 5.0 * 10.0^(-0.5) = 5.0 * 0.316227766 = 1.58113883
        // yield_term = 10.0 / 10.0 = 1.0
        // viscosity = 1.0 + 1.58113883 = 2.58113883

        assert_relative_eq!(props.dynamic_viscosity, 2.5811388300841898);
        assert!(!fluid.is_temperature_dependent());
    }

    #[test]
    fn test_herschel_bulkley_temperature_dependent() {
        let t_ref = 300.0;
        let ea_k = 5000.0; // J/mol

        let fluid = HerschelBulkley::<f64>::new(
            "Test Fluid".to_string(),
            1000.0,
            10.0, // yield stress
            5.0,  // consistency index
            0.5,  // flow behavior index
            4000.0,
            0.6,
            1500.0,
            10.0, // reference shear rate
        ).with_temperature_dependence(t_ref, Some(ea_k), None);

        assert!(fluid.is_temperature_dependent());
        assert_eq!(fluid.reference_temperature(), Some(t_ref));

        // Test at T_ref
        let props_ref = fluid.properties_at(t_ref, 101325.0).unwrap();
        assert_relative_eq!(props_ref.dynamic_viscosity, 2.5811388300841898);

        // Test at higher temperature
        // Formula: K(T) = K_ref * exp( (Ea/R) * (1/T - 1/T_ref) )
        // If T > T_ref, 1/T < 1/T_ref, so (1/T - 1/T_ref) is negative.
        // exp(negative) < 1. So K(T) < K_ref. Viscosity should decrease.

        let t_high = 350.0;
        let props_high = fluid.properties_at(t_high, 101325.0).unwrap();

        let r = 8.314462618;
        let arg = (ea_k / r) * (1.0/t_high - 1.0/t_ref);
        let k_high = 5.0 * arg.exp();

        let shear_rate = 10.0_f64;
        let power_law_term = k_high * shear_rate.powf(0.5 - 1.0);
        let yield_term = 10.0 / shear_rate;
        let expected_viscosity = yield_term + power_law_term;

        assert_relative_eq!(props_high.dynamic_viscosity, expected_viscosity);
        assert!(props_high.dynamic_viscosity < props_ref.dynamic_viscosity);
    }

    #[test]
    fn test_power_law_constant() {
        let fluid = PowerLawFluid::<f64>::new(
            "Test Fluid".to_string(),
            1000.0,
            5.0, // consistency index
            0.5, // flow behavior index
            4000.0,
            0.6,
            1500.0,
            10.0, // reference shear rate
        );

        let props = fluid.properties_at(300.0, 101325.0).unwrap();

        // Manual calc:
        // shear_rate = 10.0
        // power_law_term = 5.0 * 10.0^(0.5 - 1.0) = 5.0 * 10.0^(-0.5) = 5.0 * 0.316227766 = 1.58113883

        assert_relative_eq!(props.dynamic_viscosity, 1.5811388300841898);
        assert!(!fluid.is_temperature_dependent());
    }

    #[test]
    fn test_power_law_temperature_dependent() {
        let t_ref = 300.0;
        let ea_k = 5000.0; // J/mol

        let fluid = PowerLawFluid::<f64>::new(
            "Test Fluid".to_string(),
            1000.0,
            5.0, // consistency index
            0.5, // flow behavior index
            4000.0,
            0.6,
            1500.0,
            10.0, // reference shear rate
        )
        .with_temperature_dependence(t_ref, Some(ea_k));

        assert!(fluid.is_temperature_dependent());
        assert_eq!(fluid.reference_temperature(), Some(t_ref));

        // Test at T_ref
        let props_ref = fluid.properties_at(t_ref, 101325.0).unwrap();
        assert_relative_eq!(props_ref.dynamic_viscosity, 1.5811388300841898);

        // Test at higher temperature
        let t_high = 350.0;
        let props_high = fluid.properties_at(t_high, 101325.0).unwrap();

        let r = 8.314462618;
        let arg = (ea_k / r) * (1.0 / t_high - 1.0 / t_ref);
        let k_high = 5.0 * arg.exp();

        let shear_rate = 10.0_f64;
        let expected_viscosity = k_high * shear_rate.powf(0.5 - 1.0);

        assert_relative_eq!(props_high.dynamic_viscosity, expected_viscosity);
        assert!(props_high.dynamic_viscosity < props_ref.dynamic_viscosity);
    }
}
