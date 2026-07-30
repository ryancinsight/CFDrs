use super::constants;
use crate::error::Error;
use crate::physics::fluid::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, MassDensity, ReciprocalTime, SpecificHeatCapacity,
    ThermalConductivity, Time, Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Cross blood model (simpler alternative to Carreau-Yasuda)
///
/// # Constitutive Equation
/// ```text
/// μ(γ̇) = μ_∞ + (μ_0 - μ_∞) / (1 + (K·γ̇)^n)
/// ```
///
/// Computationally simpler than Carreau-Yasuda but provides good fit for blood.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossBlood<T> {
    /// Fluid name
    pub name: String,
    /// Blood density [kg/m³]
    pub density: MassDensity<T>,
    /// Zero-shear viscosity μ_0 [Pa·s]
    pub zero_shear_viscosity: DynamicViscosity<T>,
    /// Infinite-shear viscosity μ_∞ [Pa·s]
    pub infinite_shear_viscosity: DynamicViscosity<T>,
    /// Time constant K \[s]
    pub time_constant: Time<T>,
    /// Rate index n [-]
    pub rate_index: Dimensionless<T>,
    /// Hematocrit [-]
    pub hematocrit: Dimensionless<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
    /// Reference shear rate [1/s]
    pub reference_shear_rate: ReciprocalTime<T>,
}

impl<T: RealField + FloatElement + Copy> CrossBlood<T> {
    /// Create Cross blood model with default parameters
    pub fn normal_blood() -> Self {
        Self {
            name: "Normal Human Blood (Cross)".to_string(),
            density: MassDensity::from_base(<T as FloatElement>::from_f64(
                constants::BLOOD_DENSITY,
            )),
            zero_shear_viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(
                constants::ZERO_SHEAR_VISCOSITY,
            )),
            infinite_shear_viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(
                constants::INFINITE_SHEAR_VISCOSITY,
            )),
            time_constant: Time::from_base(<T as FloatElement>::from_f64(1.007)),
            rate_index: Dimensionless::from_base(<T as FloatElement>::from_f64(1.028)),
            hematocrit: Dimensionless::from_base(<T as FloatElement>::from_f64(
                constants::NORMAL_HEMATOCRIT,
            )),
            specific_heat: SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(
                constants::BLOOD_SPECIFIC_HEAT,
            )),
            thermal_conductivity: ThermalConductivity::from_base(<T as FloatElement>::from_f64(
                constants::BLOOD_THERMAL_CONDUCTIVITY,
            )),
            speed_of_sound: Velocity::from_base(<T as FloatElement>::from_f64(
                constants::BLOOD_SPEED_OF_SOUND,
            )),
            reference_shear_rate: ReciprocalTime::from_base(<T as FloatElement>::from_f64(100.0)),
        }
    }

    /// Calculate apparent viscosity at given shear rate
    ///
    /// μ(γ̇) = μ_∞ + (μ_0 - μ_∞) / (1 + (K·γ̇)^n)
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        if shear_rate <= <T as NumericElement>::ZERO {
            return self.zero_shear_viscosity.into_base();
        }

        let k_gamma = self.time_constant.into_base() * shear_rate;
        let denominator = <T as NumericElement>::ONE
            + <T as FloatElement>::powf(k_gamma, self.rate_index.into_base());

        self.infinite_shear_viscosity.into_base()
            + (self.zero_shear_viscosity.into_base() - self.infinite_shear_viscosity.into_base())
                / denominator
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for CrossBlood<T> {
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

    /// Return shear-rate-dependent viscosity via the Cross model.
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

impl<T: RealField + FloatElement + Copy> NonNewtonianFluid<T> for CrossBlood<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(CrossBlood::apparent_viscosity(self, shear_rate))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_cross_blood_limits() {
        let blood = CrossBlood::<f64>::normal_blood();

        // Zero shear → μ_0
        let mu_zero = blood.apparent_viscosity(0.0);
        assert_relative_eq!(mu_zero, constants::ZERO_SHEAR_VISCOSITY, epsilon = 1e-10);

        // High shear → μ_∞
        let mu_high = blood.apparent_viscosity(10000.0);
        assert!(
            mu_high < constants::ZERO_SHEAR_VISCOSITY / 2.0,
            "High shear viscosity should be significantly reduced"
        );
    }

    #[test]
    fn normal_blood_preserves_typed_fixed_metrics() {
        let blood = CrossBlood::<f64>::normal_blood();

        assert_eq!(blood.density.into_base(), constants::BLOOD_DENSITY);
        assert_eq!(
            blood.zero_shear_viscosity.into_base(),
            constants::ZERO_SHEAR_VISCOSITY
        );
        assert_eq!(blood.time_constant.into_base(), 1.007);
        assert_eq!(blood.rate_index.into_base(), 1.028);
        assert_eq!(blood.hematocrit.into_base(), constants::NORMAL_HEMATOCRIT);
    }
}
