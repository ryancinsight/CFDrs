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
///
/// Temperature dependence (optional):
/// K(T) = K_ref * exp(Ea_k * (1/T - 1/T_ref))
/// τ₀(T) = τ₀_ref * exp(Ea_y * (1/T - 1/T_ref))
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

    // Temperature dependence parameters
    /// Reference temperature for temperature-dependent properties [K]
    #[serde(default)]
    pub reference_temperature: Option<T>,
    /// Activation energy coefficient for consistency index (Ea_k/R) [K]
    #[serde(default)]
    pub consistency_activation_energy_coeff: Option<T>,
    /// Activation energy coefficient for yield stress (Ea_y/R) [K]
    #[serde(default)]
    pub yield_stress_activation_energy_coeff: Option<T>,
    /// Activation energy coefficient for flow behavior index (Ea_n/R) [K]
    #[serde(default)]
    pub flow_behavior_activation_energy_coeff: Option<T>,
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
            consistency_activation_energy_coeff: None,
            yield_stress_activation_energy_coeff: None,
            flow_behavior_activation_energy_coeff: None,
        }
    }

    /// Set temperature dependence parameters
    pub fn with_temperature_dependence(
        mut self,
        reference_temperature: T,
        consistency_activation_energy_coeff: Option<T>,
        yield_stress_activation_energy_coeff: Option<T>,
        flow_behavior_activation_energy_coeff: Option<T>,
    ) -> Self {
        self.reference_temperature = Some(reference_temperature);
        self.consistency_activation_energy_coeff = consistency_activation_energy_coeff;
        self.yield_stress_activation_energy_coeff = yield_stress_activation_energy_coeff;
        self.flow_behavior_activation_energy_coeff = flow_behavior_activation_energy_coeff;
        self
    }

    /// Helper method to calculate Arrhenius-type temperature dependence
    fn calculate_arrhenius_term(
        &self,
        temperature: T,
        base_value: T,
        activation_energy_coeff: Option<T>,
    ) -> T {
        match (self.reference_temperature, activation_energy_coeff) {
            (Some(t_ref), Some(ea)) if temperature > T::zero() => {
                // Arrhenius-type dependence: Val(T) = Val_ref * exp(Ea * (1/T - 1/T_ref))
                // Note: ea here is actually (Ea/R)
                let exponent = ea * (T::one() / temperature - T::one() / t_ref);
                base_value * exponent.exp()
            }
            _ => base_value,
        }
    }

    /// Calculate temperature-dependent Consistency Index K(T)
    pub fn consistency_index_at(&self, temperature: T) -> T {
        self.calculate_arrhenius_term(
            temperature,
            self.consistency_index,
            self.consistency_activation_energy_coeff,
        )
    }

    /// Calculate temperature-dependent Yield Stress τ₀(T)
    pub fn yield_stress_at(&self, temperature: T) -> T {
        self.calculate_arrhenius_term(
            temperature,
            self.yield_stress,
            self.yield_stress_activation_energy_coeff,
        )
    }

    /// Calculate temperature-dependent Flow Behavior Index n(T)
    pub fn flow_behavior_index_at(&self, temperature: T) -> T {
        self.calculate_arrhenius_term(
            temperature,
            self.flow_behavior_index,
            self.flow_behavior_activation_energy_coeff,
        )
    }

    /// Calculate apparent viscosity at given shear rate and temperature
    ///
    /// If temperature is None, uses reference properties.
    pub fn apparent_viscosity_at_temp(&self, shear_rate: T, temperature: Option<T>) -> T {
        if shear_rate <= T::zero() {
            return T::from_f64(1e6).unwrap_or_else(T::one);
        }

        let k = if let Some(t) = temperature {
            self.consistency_index_at(t)
        } else {
            self.consistency_index
        };

        let y_stress = if let Some(t) = temperature {
            self.yield_stress_at(t)
        } else {
            self.yield_stress
        };

        let n = if let Some(t) = temperature {
            self.flow_behavior_index_at(t)
        } else {
            self.flow_behavior_index
        };

        let power_law_term = k * shear_rate.powf(n - T::one());
        y_stress / shear_rate + power_law_term
    }

    /// Calculate apparent viscosity at given shear rate (using stored reference properties)
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        self.apparent_viscosity_at_temp(shear_rate, None)
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for HerschelBulkley<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let apparent_viscosity = self.apparent_viscosity_at_temp(self.reference_shear_rate, Some(temperature));

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
            && (self.consistency_activation_energy_coeff.is_some()
                || self.yield_stress_activation_energy_coeff.is_some()
                || self.flow_behavior_activation_energy_coeff.is_some())
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

    #[test]
    fn test_power_law_fluid_viscosity() {
        let fluid = PowerLawFluid::<f64> {
            name: "Blood".to_string(),
            density: 1060.0,
            consistency_index: 0.035, // Pa·s^n
            flow_behavior_index: 0.6, // Shear-thinning
            specific_heat: 3600.0,
            thermal_conductivity: 0.5,
            speed_of_sound: 1540.0,
            reference_shear_rate: 100.0,
        };

        // At shear rate 100 s^-1:
        // mu = K * gamma^(n-1) = 0.035 * 100^(0.6-1) = 0.035 * 100^-0.4
        // 100^-0.4 = (10^2)^-0.4 = 10^-0.8 = 0.158489
        // mu = 0.035 * 0.158489 = 0.005547
        let mu = fluid.apparent_viscosity(100.0);
        assert!((mu - 0.005547).abs() < 1e-6);
    }

    #[test]
    fn test_bingham_plastic_yield() {
        let fluid = BinghamPlastic::<f64> {
            name: "Toothpaste".to_string(),
            density: 1600.0,
            yield_stress: 200.0,      // Pa
            plastic_viscosity: 10.0,  // Pa·s
            specific_heat: 2000.0,
            thermal_conductivity: 0.5,
            speed_of_sound: 1500.0,
            reference_shear_rate: 1.0,
        };

        assert!(!fluid.is_yielded(150.0));
        assert!(fluid.is_yielded(250.0));

        // Viscosity at shear rate 1.0
        // mu = mu_p + tau_0/gamma = 10 + 200/1 = 210
        assert!((fluid.apparent_viscosity(1.0) - 210.0).abs() < 1e-10);
    }

    #[test]
    fn test_herschel_bulkley_viscosity() {
        let fluid = HerschelBulkley::<f64> {
            name: "Drilling Mud".to_string(),
            density: 1200.0,
            yield_stress: 5.0,
            consistency_index: 0.8,
            flow_behavior_index: 0.7,
            specific_heat: 2000.0,
            thermal_conductivity: 0.6,
            speed_of_sound: 1400.0,
            reference_shear_rate: 10.0,
            reference_temperature: None,
            consistency_activation_energy_coeff: None,
            yield_stress_activation_energy_coeff: None,
            flow_behavior_activation_energy_coeff: None,
        };

        // At shear rate 10.0:
        // mu = tau_0/gamma + K * gamma^(n-1)
        // mu = 5.0/10.0 + 0.8 * 10.0^(0.7-1)
        // mu = 0.5 + 0.8 * 10^-0.3
        // 10^-0.3 = 0.501187
        // mu = 0.5 + 0.8 * 0.501187 = 0.5 + 0.40095 = 0.90095
        let mu = fluid.apparent_viscosity(10.0);
        assert!((mu - 0.90095).abs() < 1e-5);
    }

    #[test]
    fn test_herschel_bulkley_temperature_dependence() {
        let fluid = HerschelBulkley::<f64> {
            name: "Temp Dependent Mud".to_string(),
            density: 1200.0,
            yield_stress: 10.0,
            consistency_index: 1.0,
            flow_behavior_index: 1.0, // Bingham-like for simplicity of calc
            specific_heat: 2000.0,
            thermal_conductivity: 0.6,
            speed_of_sound: 1400.0,
            reference_shear_rate: 1.0,
            reference_temperature: None,
            consistency_activation_energy_coeff: None,
            yield_stress_activation_energy_coeff: None,
            flow_behavior_activation_energy_coeff: None,
        }
        .with_temperature_dependence(300.0, Some(1000.0), Some(500.0), Some(200.0)); // Ea/R values

        // Test at reference temperature (300K)
        // Properties should be equal to base values
        let k_ref = fluid.consistency_index_at(300.0);
        let y_ref = fluid.yield_stress_at(300.0);
        let n_ref = fluid.flow_behavior_index_at(300.0);
        assert!((k_ref - 1.0).abs() < 1e-10);
        assert!((y_ref - 10.0).abs() < 1e-10);
        assert!((n_ref - 1.0).abs() < 1e-10);

        // Test at higher temperature (350K)
        // K(350) = 1.0 * exp(1000 * (1/350 - 1/300))
        // 1/350 - 1/300 = 0.002857 - 0.003333 = -0.000476
        // exp(1000 * -0.000476) = exp(-0.47619) = 0.6211
        let k_hot = fluid.consistency_index_at(350.0);
        let expected_k = 1.0 * (1000.0f64 * (1.0/350.0 - 1.0/300.0)).exp();
        assert!((k_hot - expected_k).abs() < 1e-10);
        assert!(k_hot < 1.0); // Viscosity should decrease

        // Yield stress with different activation energy
        let y_hot = fluid.yield_stress_at(350.0);
        let expected_y = 10.0 * (500.0f64 * (1.0/350.0 - 1.0/300.0)).exp();
        assert!((y_hot - expected_y).abs() < 1e-10);
        assert!(y_hot < 10.0);

        // Flow behavior index
        let n_hot = fluid.flow_behavior_index_at(350.0);
        let expected_n = 1.0 * (200.0f64 * (1.0/350.0 - 1.0/300.0)).exp();
        assert!((n_hot - expected_n).abs() < 1e-10);

        // Apparent viscosity via FluidTrait
        // At shear rate 1.0 (ref), T=350
        // mu = y_hot / 1.0 + k_hot * 1.0^(n_hot - 1)
        // 1.0^(n-1) is always 1.0 regardless of n
        let state = fluid.properties_at(350.0, 101325.0).unwrap();
        assert!((state.dynamic_viscosity - (y_hot + k_hot)).abs() < 1e-10);
    }
}
