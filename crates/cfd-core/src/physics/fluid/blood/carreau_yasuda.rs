use super::constants;
use crate::error::Error;
use crate::physics::fluid::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Carreau-Yasuda blood model for wide shear rate range
///
/// The Carreau-Yasuda model accurately describes blood viscosity across the full
/// physiological shear rate range (0.01 - 1000 s⁻¹), capturing both the low-shear
/// plateau and high-shear limiting behavior.
///
/// # Constitutive Equation
/// ```text
/// μ(γ̇) = μ_∞ + (μ_0 - μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
/// ```
///
/// # Parameters
/// - μ_0 = 0.056 Pa·s (zero-shear viscosity)
/// - μ_∞ = 0.00345 Pa·s (infinite-shear viscosity)
/// - λ = 3.313 s (relaxation time)
/// - n = 0.3568 (power-law index)
/// - a = 2.0 (transition parameter)
///
/// # Reference
/// Cho, Y.I., Kensey, K.R. (1991) "Effects of the non-Newtonian viscosity of blood
/// on flows in a diseased arterial vessel. Part 1: Steady flows"
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CarreauYasudaBlood<T: RealField + Copy> {
    /// Blood density [kg/m³]
    pub density: T,
    /// Zero-shear viscosity μ_0 [Pa·s]
    pub zero_shear_viscosity: T,
    /// Infinite-shear viscosity μ_∞ [Pa·s]
    pub infinite_shear_viscosity: T,
    /// Relaxation time λ [s]
    pub relaxation_time: T,
    /// Power-law index n [-]
    pub power_law_index: T,
    /// Transition parameter a [-]
    pub transition_parameter: T,
    /// Hematocrit (volume fraction of RBCs) [-]
    pub hematocrit: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
    /// Reference shear rate for default viscosity calculation [1/s]
    pub reference_shear_rate: T,
}

impl<T: RealField + FromPrimitive + Copy> CarreauYasudaBlood<T> {
    /// Create Carreau-Yasuda blood model with literature-validated parameters
    ///
    /// Uses parameters from Cho & Kensey (1991):
    /// - μ_0 = 0.056 Pa·s
    /// - μ_∞ = 0.00345 Pa·s
    /// - λ = 3.313 s
    /// - n = 0.3568
    /// - a = 2.0
    pub fn normal_blood() -> Self {
        Self {
            density: T::from_f64(constants::BLOOD_DENSITY).unwrap_or_else(num_traits::Zero::zero),
            zero_shear_viscosity: T::from_f64(constants::ZERO_SHEAR_VISCOSITY).unwrap_or_else(num_traits::Zero::zero),
            infinite_shear_viscosity: T::from_f64(constants::INFINITE_SHEAR_VISCOSITY).unwrap_or_else(num_traits::Zero::zero),
            relaxation_time: T::from_f64(constants::CARREAU_LAMBDA).unwrap_or_else(num_traits::Zero::zero),
            power_law_index: T::from_f64(constants::CARREAU_N).unwrap_or_else(num_traits::Zero::zero),
            transition_parameter: T::from_f64(constants::CARREAU_A).unwrap_or_else(num_traits::Zero::zero),
            hematocrit: T::from_f64(constants::NORMAL_HEMATOCRIT).unwrap_or_else(num_traits::Zero::zero),
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap_or_else(num_traits::Zero::zero),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap_or_else(num_traits::Zero::zero),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap_or_else(num_traits::Zero::zero),
            reference_shear_rate: T::from_f64(100.0).unwrap_or_else(num_traits::Zero::zero),
        }
    }

    /// Create Carreau-Yasuda blood model with custom parameters
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        density: T,
        zero_shear_viscosity: T,
        infinite_shear_viscosity: T,
        relaxation_time: T,
        power_law_index: T,
        transition_parameter: T,
        hematocrit: T,
    ) -> Self {
        Self {
            density,
            zero_shear_viscosity,
            infinite_shear_viscosity,
            relaxation_time,
            power_law_index,
            transition_parameter,
            hematocrit,
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap_or_else(num_traits::Zero::zero),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap_or_else(num_traits::Zero::zero),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap_or_else(num_traits::Zero::zero),
            reference_shear_rate: T::from_f64(100.0).unwrap_or_else(num_traits::Zero::zero),
        }
    }

    /// Calculate apparent viscosity at given shear rate
    ///
    /// # Mathematical Formula
    /// ```text
    /// μ(γ̇) = μ_∞ + (μ_0 - μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
    /// ```
    ///
    /// # Limiting Behavior
    /// - As γ̇ → 0: μ → μ_0 (zero-shear viscosity)
    /// - As γ̇ → ∞: μ → μ_∞ (infinite-shear viscosity)
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        let one = T::one();

        // Handle zero shear rate
        if shear_rate <= T::zero() {
            return self.zero_shear_viscosity;
        }

        // λγ̇
        let lambda_gamma = self.relaxation_time * shear_rate;

        // (λγ̇)^a
        let lambda_gamma_a = lambda_gamma.powf(self.transition_parameter);

        // 1 + (λγ̇)^a
        let bracketed = one + lambda_gamma_a;

        // (n-1)/a
        let exponent = (self.power_law_index - one) / self.transition_parameter;

        // [1 + (λγ̇)^a]^((n-1)/a)
        let shear_factor = bracketed.powf(exponent);

        // μ_∞ + (μ_0 - μ_∞) · shear_factor
        self.infinite_shear_viscosity
            + (self.zero_shear_viscosity - self.infinite_shear_viscosity) * shear_factor
    }

    /// Validate model parameters
    pub fn validate(&self) -> Result<(), Error> {
        if self.density <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        if self.zero_shear_viscosity <= T::zero() {
            return Err(Error::InvalidInput(
                "Zero-shear viscosity must be positive".to_string(),
            ));
        }
        if self.infinite_shear_viscosity <= T::zero() {
            return Err(Error::InvalidInput(
                "Infinite-shear viscosity must be positive".to_string(),
            ));
        }
        if self.infinite_shear_viscosity >= self.zero_shear_viscosity {
            return Err(Error::InvalidInput(
                "Infinite-shear viscosity must be less than zero-shear viscosity".to_string(),
            ));
        }
        if self.relaxation_time <= T::zero() {
            return Err(Error::InvalidInput(
                "Relaxation time must be positive".to_string(),
            ));
        }
        if self.power_law_index <= T::zero() || self.power_law_index >= T::one() {
            return Err(Error::InvalidInput(
                "Power-law index must be in (0, 1) for shear-thinning".to_string(),
            ));
        }
        if self.transition_parameter <= T::zero() {
            return Err(Error::InvalidInput(
                "Transition parameter must be positive".to_string(),
            ));
        }
        Ok(())
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for CarreauYasudaBlood<T> {
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

    fn name(&self) -> &'static str {
        "Carreau-Yasuda Blood Model"
    }

    /// Return shear-rate-dependent viscosity via the Carreau-Yasuda model.
    /// This override ensures correct non-Newtonian apparent viscosity
    /// when called through the unified `Fluid::viscosity_at_shear` interface.
    fn viscosity_at_shear(
        &self,
        shear_rate: T,
        _temperature: T,
        _pressure: T,
    ) -> Result<T, Error> {
        Ok(self.apparent_viscosity(shear_rate))
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for CarreauYasudaBlood<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        CarreauYasudaBlood::apparent_viscosity(self, shear_rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_carreau_yasuda_normal_blood() {
        let blood = CarreauYasudaBlood::<f64>::normal_blood();
        assert_eq!(blood.zero_shear_viscosity, constants::ZERO_SHEAR_VISCOSITY);
        assert_eq!(
            blood.infinite_shear_viscosity,
            constants::INFINITE_SHEAR_VISCOSITY
        );
        assert_eq!(blood.relaxation_time, constants::CARREAU_LAMBDA);
        assert_eq!(blood.power_law_index, constants::CARREAU_N);
    }

    #[test]
    fn test_carreau_yasuda_limits() {
        let blood = CarreauYasudaBlood::<f64>::normal_blood();

        // Zero shear rate → μ_0
        let mu_zero = blood.apparent_viscosity(0.0);
        assert_relative_eq!(mu_zero, constants::ZERO_SHEAR_VISCOSITY, epsilon = 1e-10);

        // Very high shear rate → μ_∞
        let mu_high = blood.apparent_viscosity(100000.0);
        assert!(
            (mu_high - constants::INFINITE_SHEAR_VISCOSITY).abs() < 0.0001,
            "High shear viscosity {} should approach μ_∞ = {}",
            mu_high,
            constants::INFINITE_SHEAR_VISCOSITY
        );
    }

    #[test]
    fn test_carreau_yasuda_shear_thinning() {
        let blood = CarreauYasudaBlood::<f64>::normal_blood();

        let shear_rates = [0.1, 1.0, 10.0, 100.0, 1000.0];
        let mut prev_viscosity = f64::MAX;

        for &gamma in &shear_rates {
            let mu = blood.apparent_viscosity(gamma);
            assert!(
                mu < prev_viscosity,
                "Viscosity should decrease: μ({}) = {} should be less than {}",
                gamma,
                mu,
                prev_viscosity
            );
            prev_viscosity = mu;
        }
    }

    #[test]
    fn test_carreau_yasuda_cho_kensey_validation() {
        // Validate against Cho & Kensey (1991) Table 1
        let blood = CarreauYasudaBlood::<f64>::normal_blood();

        // At γ̇ = 1 s⁻¹: μ ≈ 0.035 Pa·s
        let mu_1 = blood.apparent_viscosity(1.0);
        assert!(
            mu_1 > 0.02 && mu_1 < 0.05,
            "μ(1 s⁻¹) = {} Pa·s should be ~0.035 Pa·s",
            mu_1
        );

        // At γ̇ = 100 s⁻¹: μ ≈ 0.005 Pa·s
        let mu_100 = blood.apparent_viscosity(100.0);
        assert!(
            mu_100 > 0.003 && mu_100 < 0.01,
            "μ(100 s⁻¹) = {} Pa·s should be ~0.005 Pa·s",
            mu_100
        );
    }
}
