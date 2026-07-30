//! Herschel–Bulkley fluid model (generalised Bingham plastic).
//!
//! τ = τ₀ + K · (γ̇)^n  for τ > τ₀

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use crate::physics::constants::physics::universal::GAS_CONSTANT;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, MassDensity, MolarEnergy, Pressure, ReciprocalTime,
    SpecificHeatCapacity, ThermalConductivity, ThermodynamicTemperature, Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
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
    pub density: MassDensity<T>,
    /// Yield stress τ₀ \[Pa]
    pub yield_stress: Pressure<T>,
    /// Consistency index K; its units are `Pa·s^n` and depend on the runtime
    /// flow exponent, so it remains scalar formula data.
    pub consistency_index: T,
    /// Flow behavior index n [-]
    pub flow_behavior_index: Dimensionless<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
    /// Reference shear rate [1/s]
    pub reference_shear_rate: ReciprocalTime<T>,
    /// Reference temperature for properties \[K]
    #[serde(default)]
    pub reference_temperature: Option<ThermodynamicTemperature<T>>,
    /// Activation energy for consistency index [J/mol]
    #[serde(default)]
    pub activation_energy_k: Option<MolarEnergy<T>>,
    /// Activation energy for yield stress [J/mol]
    #[serde(default)]
    pub activation_energy_tau: Option<MolarEnergy<T>>,
}

impl<T: RealField + FloatElement + Copy> HerschelBulkley<T> {
    /// Create a new Herschel-Bulkley fluid
    pub fn new(
        name: String,
        density: MassDensity<T>,
        yield_stress: Pressure<T>,
        consistency_index: T,
        flow_behavior_index: Dimensionless<T>,
        specific_heat: SpecificHeatCapacity<T>,
        thermal_conductivity: ThermalConductivity<T>,
        speed_of_sound: Velocity<T>,
        reference_shear_rate: ReciprocalTime<T>,
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
        reference_temperature: ThermodynamicTemperature<T>,
        activation_energy_k: Option<MolarEnergy<T>>,
        activation_energy_tau: Option<MolarEnergy<T>>,
    ) -> Self {
        self.reference_temperature = Some(reference_temperature);
        self.activation_energy_k = activation_energy_k;
        self.activation_energy_tau = activation_energy_tau;
        self
    }

    /// Calculate apparent viscosity at a scalar shear rate at the formula
    /// boundary.
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= <T as NumericElement>::ZERO {
            return <T as FloatElement>::from_f64(1e6);
        }

        let power_law_term = self.consistency_index
            * <T as FloatElement>::powf(
                shear_rate,
                self.flow_behavior_index.into_base() - <T as NumericElement>::ONE,
            );
        self.yield_stress.into_base() / shear_rate + power_law_term
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for HerschelBulkley<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let (k_adj, tau_adj) = if let Some(t_ref) = self.reference_temperature {
            let r_gas = <T as FloatElement>::from_f64(GAS_CONSTANT);
            let inv_t = <T as NumericElement>::ONE / temperature;
            let inv_t_ref = <T as NumericElement>::ONE / t_ref.into_base();
            let diff_inv_t = inv_t - inv_t_ref;

            let k = if let Some(ea_k) = self.activation_energy_k {
                let arg = (ea_k.into_base() / r_gas) * diff_inv_t;
                self.consistency_index * <T as FloatElement>::exp(arg)
            } else {
                self.consistency_index
            };

            let tau = if let Some(ea_tau) = self.activation_energy_tau {
                let arg = (ea_tau.into_base() / r_gas) * diff_inv_t;
                self.yield_stress.into_base() * <T as FloatElement>::exp(arg)
            } else {
                self.yield_stress.into_base()
            };

            (k, tau)
        } else {
            (self.consistency_index, self.yield_stress.into_base())
        };

        let shear_rate = self.reference_shear_rate.into_base();
        let apparent_viscosity = if shear_rate <= <T as NumericElement>::ZERO {
            <T as FloatElement>::from_f64(1e6)
        } else {
            let power_law_term = k_adj
                * <T as FloatElement>::powf(
                    shear_rate,
                    self.flow_behavior_index.into_base() - <T as NumericElement>::ONE,
                );
            tau_adj / shear_rate + power_law_term
        };

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

    fn is_temperature_dependent(&self) -> bool {
        self.reference_temperature.is_some()
    }

    fn reference_temperature(&self) -> Option<ThermodynamicTemperature<T>> {
        self.reference_temperature
    }
}

impl<T: RealField + FloatElement + Copy> NonNewtonianFluid<T> for HerschelBulkley<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(HerschelBulkley::apparent_viscosity(self, shear_rate))
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<Pressure<T>> {
        Some(self.yield_stress)
    }
}
