use super::constants;
use crate::error::Error;
use crate::physics::fluid::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use nalgebra::RealField;
use num_traits::FromPrimitive;
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
pub struct CrossBlood<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Blood density [kg/m³]
    pub density: T,
    /// Zero-shear viscosity μ_0 [Pa·s]
    pub zero_shear_viscosity: T,
    /// Infinite-shear viscosity μ_∞ [Pa·s]
    pub infinite_shear_viscosity: T,
    /// Time constant K [s]
    pub time_constant: T,
    /// Rate index n [-]
    pub rate_index: T,
    /// Hematocrit [-]
    pub hematocrit: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
    /// Reference shear rate [1/s]
    pub reference_shear_rate: T,
}

impl<T: RealField + FromPrimitive + Copy> CrossBlood<T> {
    /// Create Cross blood model with default parameters
    pub fn normal_blood() -> Self {
        Self {
            name: "Normal Human Blood (Cross)".to_string(),
            density: T::from_f64(constants::BLOOD_DENSITY).unwrap_or_else(num_traits::Zero::zero),
            zero_shear_viscosity: T::from_f64(constants::ZERO_SHEAR_VISCOSITY).unwrap_or_else(num_traits::Zero::zero),
            infinite_shear_viscosity: T::from_f64(constants::INFINITE_SHEAR_VISCOSITY).unwrap_or_else(num_traits::Zero::zero),
            time_constant: T::from_f64(1.007).unwrap_or_else(num_traits::Zero::zero), // Fitted for blood
            rate_index: T::from_f64(1.028).unwrap_or_else(num_traits::Zero::zero),    // Fitted for blood
            hematocrit: T::from_f64(constants::NORMAL_HEMATOCRIT).unwrap_or_else(num_traits::Zero::zero),
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap_or_else(num_traits::Zero::zero),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap_or_else(num_traits::Zero::zero),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap_or_else(num_traits::Zero::zero),
            reference_shear_rate: T::from_f64(100.0).unwrap_or_else(num_traits::Zero::zero),
        }
    }

    /// Calculate apparent viscosity at given shear rate
    ///
    /// μ(γ̇) = μ_∞ + (μ_0 - μ_∞) / (1 + (K·γ̇)^n)
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        let one = T::one();

        if shear_rate <= T::zero() {
            return self.zero_shear_viscosity;
        }

        let k_gamma = self.time_constant * shear_rate;
        let denominator = one + k_gamma.powf(self.rate_index);

        self.infinite_shear_viscosity
            + (self.zero_shear_viscosity - self.infinite_shear_viscosity) / denominator
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for CrossBlood<T> {
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

    /// Return shear-rate-dependent viscosity via the Cross model.
    fn viscosity_at_shear(
        &self,
        shear_rate: T,
        _temperature: T,
        _pressure: T,
    ) -> Result<T, Error> {
        Ok(self.apparent_viscosity(shear_rate))
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for CrossBlood<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        CrossBlood::apparent_viscosity(self, shear_rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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
}
