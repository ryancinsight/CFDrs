//! Carreau–Yasuda fluid model.
//!
//! μ = μ_inf + (μ₀ − μ_inf) · [1 + (λ · γ̇)^a]^((n−1)/a)

use super::super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, MassDensity, SpecificHeatCapacity, ThermalConductivity, Time,
    Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Carreau-Yasuda fluid model
///
/// Common model for blood rheology.
/// μ = `μ_inf` + (`μ_0` - `μ_inf`) * [1 + (λ * γ̇)^a ]^((n-1)/a)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CarreauYasuda<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: MassDensity<T>,
    /// Zero-shear viscosity μ₀ [Pa·s]
    pub viscosity_zero: DynamicViscosity<T>,
    /// Infinite-shear viscosity `μ_inf` [Pa·s]
    pub viscosity_inf: DynamicViscosity<T>,
    /// Relaxation time λ \[s]
    pub lambda: Time<T>,
    /// Power law index n [-]
    pub power_index: Dimensionless<T>,
    /// Yasuda parameter a [-] (default 2.0 for Carreau model)
    pub yasuda_index: Dimensionless<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
}

impl<T: RealField + FloatElement + Copy> CarreauYasuda<T> {
    /// Create a new Carreau-Yasuda fluid
    pub fn new(
        name: String,
        density: MassDensity<T>,
        viscosity_zero: DynamicViscosity<T>,
        viscosity_inf: DynamicViscosity<T>,
        lambda: Time<T>,
        power_index: Dimensionless<T>,
        yasuda_index: Dimensionless<T>,
        specific_heat: SpecificHeatCapacity<T>,
        thermal_conductivity: ThermalConductivity<T>,
        speed_of_sound: Velocity<T>,
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
        let one = <T as NumericElement>::ONE;
        let term1 = <T as FloatElement>::powf(
            self.lambda.into_base() * shear_rate,
            self.yasuda_index.into_base(),
        );
        let base = one + term1;
        let exponent = (self.power_index.into_base() - one) / self.yasuda_index.into_base();

        self.viscosity_inf.into_base()
            + (self.viscosity_zero.into_base() - self.viscosity_inf.into_base())
                * <T as FloatElement>::powf(base, exponent)
    }

    /// Standard blood parameters (approximate)
    #[must_use]
    pub fn blood() -> Self {
        Self::new(
            "Blood".to_string(),
            MassDensity::from_base(<T as FloatElement>::from_f64(1060.0)),
            DynamicViscosity::from_base(<T as FloatElement>::from_f64(0.056)),
            DynamicViscosity::from_base(<T as FloatElement>::from_f64(0.0035)),
            Time::from_base(<T as FloatElement>::from_f64(3.313)),
            Dimensionless::from_base(<T as FloatElement>::from_f64(0.3568)),
            Dimensionless::from_base(<T as FloatElement>::from_f64(2.0)),
            SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(3600.0)),
            ThermalConductivity::from_base(<T as FloatElement>::from_f64(0.5)),
            Velocity::from_base(<T as FloatElement>::from_f64(1540.0)),
        )
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for CarreauYasuda<T> {
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
    ) -> Result<DynamicViscosity<T>, Error> {
        Ok(DynamicViscosity::from_base(
            self.apparent_viscosity(shear_rate),
        ))
    }
}

impl<T: RealField + FloatElement + Copy> NonNewtonianFluid<T> for CarreauYasuda<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(CarreauYasuda::apparent_viscosity(self, shear_rate))
    }
}
