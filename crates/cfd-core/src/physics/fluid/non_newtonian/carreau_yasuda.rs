//! Carreau–Yasuda fluid model.
//!
//! μ = μ_inf + (μ₀ − μ_inf) · [1 + (λ · γ̇)^a]^((n−1)/a)

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Carreau-Yasuda fluid model
///
/// Common model for blood rheology.
/// μ = μ_inf + (μ_0 - μ_inf) * [1 + (λ * γ̇)^a ]^((n-1)/a)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CarreauYasuda<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Zero-shear viscosity μ₀ [Pa·s]
    pub viscosity_zero: T,
    /// Infinite-shear viscosity μ_inf [Pa·s]
    pub viscosity_inf: T,
    /// Relaxation time λ [s]
    pub lambda: T,
    /// Power law index n [-]
    pub power_index: T,
    /// Yasuda parameter a [-] (default 2.0 for Carreau model)
    pub yasuda_index: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + FromPrimitive + Copy> CarreauYasuda<T> {
    /// Create a new Carreau-Yasuda fluid
    pub fn new(
        name: String,
        density: T,
        viscosity_zero: T,
        viscosity_inf: T,
        lambda: T,
        power_index: T,
        yasuda_index: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
    ) -> Self {
        Self {
            name,
            density,
            viscosity_zero,
            viscosity_inf,
            lambda,
            power_index,
            yasuda_index,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
        }
    }

    /// Calculate apparent viscosity at given shear rate
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        let one = T::one();
        let term1 = (self.lambda * shear_rate).powf(self.yasuda_index);
        let base = one + term1;
        let exponent = (self.power_index - one) / self.yasuda_index;

        self.viscosity_inf + (self.viscosity_zero - self.viscosity_inf) * base.powf(exponent)
    }

    /// Standard blood parameters (approximate)
    pub fn blood() -> Self {
        Self::new(
            "Blood".to_string(),
            T::from_f64(1060.0).unwrap(),
            T::from_f64(0.056).unwrap(),
            T::from_f64(0.0035).unwrap(),
            T::from_f64(3.313).unwrap(),
            T::from_f64(0.3568).unwrap(),
            T::from_f64(2.0).unwrap(),
            T::from_f64(3600.0).unwrap(),
            T::from_f64(0.5).unwrap(),
            T::from_f64(1540.0).unwrap(),
        )
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for CarreauYasuda<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: self.viscosity_zero,
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
    ) -> Result<T, Error> {
        Ok(self.apparent_viscosity(shear_rate))
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for CarreauYasuda<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        CarreauYasuda::apparent_viscosity(self, shear_rate)
    }
}
