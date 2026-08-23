use super::constants;
use crate::error::Error;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, MassDensity, Pressure, ReciprocalTime, SpecificHeatCapacity,
    ThermalConductivity, ThermodynamicTemperature, Velocity,
};

// ── Temperature correction ────────────────────────────────────────────────────

/// Andrade viscosity correction factor for whole blood temperature dependence.
///
/// Blood viscosity decreases with temperature above 37 °C (e.g., post-cavitation
/// local heating) and increases below 37 °C (e.g., therapeutic hypothermia).
///
/// ```text
///   factor(T) = exp(B × (1/T − 1/T_ref))
/// ```
///
/// where B = 1500 K (empirical constant for whole blood, Merrill 1969),
/// T_ref = 310 K (37 °C reference temperature).
///
/// # Arguments
/// * `temp_c` — Blood temperature [°C]
///
/// # Returns
/// Dimensionless viscosity scaling factor (1.0 at 37 °C, > 1 below, < 1 above).
///
/// # Example
///
/// ```
/// use cfd_core::physics::fluid::blood::temperature_viscosity_factor;
/// // At 37 °C, factor = 1.0 (no correction)
/// assert!((temperature_viscosity_factor(37.0) - 1.0).abs() < 1e-10);
/// // At 42 °C (post-cavitation), viscosity is reduced by ~7%
/// assert!(temperature_viscosity_factor(42.0) < 1.0);
/// // At 30 °C (hypothermia), viscosity is increased by ~12%
/// assert!(temperature_viscosity_factor(30.0) > 1.0);
/// ```
#[inline]
#[must_use]
pub fn temperature_viscosity_factor(temp_c: f64) -> f64 {
    let b_k = super::constants::ANDRADE_B_BLOOD;
    let t_ref_k = super::constants::BODY_TEMPERATURE_K;
    let t_k = temp_c + crate::physics::constants::physics::thermo::CELSIUS_TO_KELVIN;
    (b_k * (1.0 / t_k - 1.0 / t_ref_k)).exp()
}
use crate::physics::fluid::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Casson blood model for hemodynamic simulations
///
/// The Casson model is the standard rheological model for blood flow in larger vessels
/// (arteries and veins). It captures both the yield stress behavior and shear-thinning.
///
/// # Constitutive Equation
/// ```text
/// √τ = √τ_y + √(μ_∞ · γ̇)    for τ > τ_y
/// γ̇ = 0                      for τ ≤ τ_y
/// ```
///
/// # Apparent Viscosity
/// ```text
/// μ_app(γ̇) = τ/γ̇ = (√τ_y + √(μ_∞ · γ̇))² / γ̇
///          = μ_∞ + τ_y/γ̇ + 2√(τ_y · μ_∞ / γ̇)
/// ```
///
/// # Literature Validation
/// - Yield stress τ_y ≈ 0.0056 Pa for normal blood (H_t = 45%)
/// - Infinite shear viscosity μ_∞ ≈ 0.00345 Pa·s
/// - Reference: Merrill et al. (1969), Fung (1993)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CassonBlood<T> {
    /// Blood density [kg/m³]
    pub density: MassDensity<T>,
    /// Yield stress τ_y \[Pa]
    /// Literature value: 0.0056 Pa for normal blood (Merrill 1969)
    pub yield_stress: Pressure<T>,
    /// Infinite-shear (Casson) viscosity μ_∞ [Pa·s]
    /// Literature value: 0.00345 Pa·s (Merrill 1969)
    pub infinite_shear_viscosity: DynamicViscosity<T>,
    /// Hematocrit (volume fraction of RBCs) [-]
    pub hematocrit: Dimensionless<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
    /// Reference shear rate for default viscosity calculation [1/s]
    pub reference_shear_rate: ReciprocalTime<T>,
    /// Regularization parameter for low shear rates [1/s]
    /// Prevents singularity at γ̇ → 0
    pub regularization_shear_rate: ReciprocalTime<T>,
}

impl<T: RealField + FloatElement + Copy> CassonBlood<T> {
    /// Create Casson blood model with literature-validated parameters for normal blood
    ///
    /// Uses parameters from Merrill et al. (1969) and Fung (1993):
    /// - τ_y = 0.0056 Pa
    /// - μ_∞ = 0.00345 Pa·s
    /// - H_t = 0.45
    ///
    /// # Example
    /// ```
    /// use cfd_core::physics::fluid::blood::CassonBlood;
    ///
    /// let blood = CassonBlood::<f64>::normal_blood();
    /// let viscosity = blood.apparent_viscosity(100.0); // At γ̇ = 100 s⁻¹
    /// assert!(viscosity > 0.003 && viscosity < 0.01);
    /// ```
    pub fn normal_blood() -> Self {
        Self {
            density: MassDensity::from_base(<T as FloatElement>::from_f64(
                constants::BLOOD_DENSITY,
            )),
            yield_stress: Pressure::from_base(<T as FloatElement>::from_f64(
                constants::YIELD_STRESS,
            )),
            infinite_shear_viscosity: DynamicViscosity::from_base(<T as FloatElement>::from_f64(
                constants::INFINITE_SHEAR_VISCOSITY,
            )),
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
            regularization_shear_rate: ReciprocalTime::from_base(<T as FloatElement>::from_f64(
                0.001,
            )),
        }
    }

    /// Create Casson blood model with custom parameters
    pub fn new(
        density: MassDensity<T>,
        yield_stress: Pressure<T>,
        infinite_shear_viscosity: DynamicViscosity<T>,
        hematocrit: Dimensionless<T>,
    ) -> Self {
        Self {
            density,
            yield_stress,
            infinite_shear_viscosity,
            hematocrit,
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
            regularization_shear_rate: ReciprocalTime::from_base(<T as FloatElement>::from_f64(
                0.001,
            )),
        }
    }

    /// Create blood model with hematocrit-dependent properties
    ///
    /// Yields stress and viscosity scale with hematocrit using empirical correlations:
    /// - τ_y(H_t) = τ_y,ref · (H_t / H_t,ref)^3  (Chien 1970)
    /// - μ_∞(H_t) = μ_plasma · exp(k · H_t / (1 - H_t))  (Quemada 1978)
    pub fn with_hematocrit(hematocrit: Dimensionless<T>) -> Self {
        let h_ref = <T as FloatElement>::from_f64(constants::NORMAL_HEMATOCRIT);
        let tau_ref = <T as FloatElement>::from_f64(constants::YIELD_STRESS);
        let mu_plasma = <T as FloatElement>::from_f64(constants::PLASMA_VISCOSITY_37C);
        let hematocrit_base = hematocrit.into_base();

        // Yield stress scaling with hematocrit (Chien 1970)
        let ratio = hematocrit_base / h_ref;
        let yield_stress = tau_ref * ratio * ratio * ratio;

        // Infinite-shear viscosity from Quemada model
        let k = <T as FloatElement>::from_f64(2.5); // Intrinsic viscosity coefficient
        let one = <T as NumericElement>::ONE;
        let infinite_shear_viscosity =
            mu_plasma * <T as FloatElement>::exp(k * hematocrit_base / (one - hematocrit_base));

        Self {
            density: MassDensity::from_base(<T as FloatElement>::from_f64(
                constants::BLOOD_DENSITY,
            )),
            yield_stress: Pressure::from_base(yield_stress),
            infinite_shear_viscosity: DynamicViscosity::from_base(infinite_shear_viscosity),
            hematocrit,
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
            regularization_shear_rate: ReciprocalTime::from_base(<T as FloatElement>::from_f64(
                0.001,
            )),
        }
    }

    /// Calculate apparent viscosity at given shear rate
    ///
    /// # Mathematical Derivation
    /// From the Casson constitutive equation:
    /// ```text
    /// √τ = √τ_y + √(μ_∞ · γ̇)
    /// τ = (√τ_y + √(μ_∞ · γ̇))²
    /// μ_app = τ / γ̇
    /// ```
    ///
    /// # Regularization
    /// For numerical stability at low shear rates, we use:
    /// ```text
    /// γ̇_eff = max(γ̇, γ̇_reg)
    /// ```
    ///
    /// # Arguments
    /// * `shear_rate` - Shear rate γ̇ [1/s]
    ///
    /// # Returns
    /// Apparent dynamic viscosity μ_app [Pa·s]
    pub fn apparent_viscosity(&self, shear_rate: T) -> T {
        // Apply regularization for numerical stability
        let gamma_eff = if shear_rate < self.regularization_shear_rate.into_base() {
            self.regularization_shear_rate.into_base()
        } else {
            shear_rate
        };

        // Casson model: μ_app = (√τ_y/√γ̇ + √μ_∞)²
        let sqrt_tau_y = <T as NumericElement>::sqrt(self.yield_stress.into_base());
        let sqrt_mu_inf = <T as NumericElement>::sqrt(self.infinite_shear_viscosity.into_base());
        let sqrt_gamma = <T as NumericElement>::sqrt(gamma_eff);

        let casson_sqrt = sqrt_tau_y / sqrt_gamma + sqrt_mu_inf;
        casson_sqrt * casson_sqrt
    }

    /// Calculate shear stress from shear rate
    ///
    /// τ = (√τ_y + √(μ_∞ · γ̇))²
    pub fn shear_stress(&self, shear_rate: T) -> T {
        self.apparent_viscosity(shear_rate) * shear_rate
    }

    /// Compute blood apparent viscosity at a given temperature.
    ///
    /// Applies the Andrade temperature correction to the Casson apparent viscosity:
    ///
    /// ```text
    ///   μ(T) = μ(37 °C) × exp(B × (1/T_ref − 1/T))
    /// ```
    ///
    /// where B ≈ 1500 K (empirical blood constant), T_ref = 310 K (37 °C).
    ///
    /// # Arguments
    /// * `shear_rate` — Wall shear rate [s⁻¹]
    /// * `temperature` — Blood temperature \[K]; if ≤ 0, defaults to 310 K (37 °C)
    pub fn apparent_viscosity_at_temp(
        &self,
        shear_rate: T,
        temperature: ThermodynamicTemperature<T>,
    ) -> T {
        let t_ref = <T as FloatElement>::from_f64(super::constants::BODY_TEMPERATURE_K);
        let b = <T as FloatElement>::from_f64(super::constants::ANDRADE_B_BLOOD);
        let temp_k = temperature.into_base();
        let t_eff = if temp_k <= <T as NumericElement>::ZERO {
            t_ref
        } else {
            temp_k
        };
        // Andrade factor: exp(B × (1/T − 1/T_ref))
        let factor = <T as FloatElement>::exp(
            b * (<T as NumericElement>::ONE / t_eff - <T as NumericElement>::ONE / t_ref),
        );
        self.apparent_viscosity(shear_rate) * factor
    }

    /// Validate model parameters
    ///
    /// # Errors
    /// Returns error if parameters are physically invalid
    pub fn validate(&self) -> Result<(), Error> {
        if self.density.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        if self.yield_stress.into_base() < <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Yield stress must be non-negative".to_string(),
            ));
        }
        if self.infinite_shear_viscosity.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Infinite-shear viscosity must be positive".to_string(),
            ));
        }
        if self.hematocrit.into_base() < <T as NumericElement>::ZERO
            || self.hematocrit.into_base() > <T as NumericElement>::ONE
        {
            return Err(Error::InvalidInput(
                "Hematocrit must be between 0 and 1".to_string(),
            ));
        }
        Ok(())
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for CassonBlood<T> {
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

    fn name(&self) -> &'static str {
        "Casson Blood Model"
    }

    /// Return shear-rate-dependent viscosity via the Casson constitutive law.
    /// This override ensures that any code calling `Fluid::viscosity_at_shear`
    /// (e.g. the 1D bifurcation junction solver) obtains the correct
    /// non-Newtonian apparent viscosity rather than a constant reference value.
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

impl<T: RealField + FloatElement + Copy> NonNewtonianFluid<T> for CassonBlood<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(CassonBlood::apparent_viscosity(self, shear_rate))
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<Pressure<T>> {
        Some(self.yield_stress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_casson_normal_blood_creation() {
        let blood = CassonBlood::<f64>::normal_blood();
        assert_eq!(blood.density.into_base(), constants::BLOOD_DENSITY);
        assert_eq!(blood.yield_stress.into_base(), constants::YIELD_STRESS);
        assert_eq!(
            blood.infinite_shear_viscosity.into_base(),
            constants::INFINITE_SHEAR_VISCOSITY
        );
        assert_eq!(blood.hematocrit.into_base(), constants::NORMAL_HEMATOCRIT);
    }

    #[test]
    fn test_casson_viscosity_limits() {
        let blood = CassonBlood::<f64>::normal_blood();

        // At very high shear rate, should approach μ_∞
        let mu_high = blood.apparent_viscosity(10000.0);
        assert!(
            (mu_high - constants::INFINITE_SHEAR_VISCOSITY).abs() < 0.001,
            "High shear viscosity {} should approach μ_∞ = {}",
            mu_high,
            constants::INFINITE_SHEAR_VISCOSITY
        );

        // At intermediate shear rate, viscosity should be between limits
        let mu_100 = blood.apparent_viscosity(100.0);
        assert!(mu_100 > constants::INFINITE_SHEAR_VISCOSITY);
        assert!(mu_100 < constants::ZERO_SHEAR_VISCOSITY);
    }

    #[test]
    fn test_casson_shear_thinning() {
        let blood = CassonBlood::<f64>::normal_blood();

        let mu_1 = blood.apparent_viscosity(1.0);
        let mu_10 = blood.apparent_viscosity(10.0);
        let mu_100 = blood.apparent_viscosity(100.0);
        let mu_1000 = blood.apparent_viscosity(1000.0);

        // Blood should be shear-thinning: viscosity decreases with shear rate
        assert!(mu_1 > mu_10, "Viscosity should decrease: μ(1) > μ(10)");
        assert!(mu_10 > mu_100, "Viscosity should decrease: μ(10) > μ(100)");
        assert!(
            mu_100 > mu_1000,
            "Viscosity should decrease: μ(100) > μ(1000)"
        );
    }

    #[test]
    fn test_casson_literature_validation() {
        // Validate against Merrill et al. (1969) Fig. 5
        // At γ̇ = 100 s⁻¹, μ_app ≈ 4.0 mPa·s for H_t = 45%
        let blood = CassonBlood::<f64>::normal_blood();
        let mu_100 = blood.apparent_viscosity(100.0);

        // Literature value: ~4.0 mPa·s at 100 s⁻¹
        let literature_value = 0.004; // Pa·s
        let tolerance = 0.5; // 50% tolerance for model vs literature

        assert!(
            (mu_100 - literature_value).abs() / literature_value < tolerance,
            "Casson viscosity at 100 s⁻¹ = {} Pa·s, expected ~{} Pa·s (within {}%)",
            mu_100,
            literature_value,
            tolerance * 100.0
        );
    }

    #[test]
    fn test_casson_validation_passes() {
        let blood = CassonBlood::<f64>::normal_blood();
        assert!(blood.validate().is_ok());

        // Invalid parameters
        let mut invalid = blood;
        invalid.density = MassDensity::from_base(-1.0);
        assert!(invalid.validate().is_err());
    }

    // ── temperature_viscosity_factor ─────────────────────────────────────────

    #[test]
    fn temperature_factor_at_37c_is_unity() {
        assert!((temperature_viscosity_factor(37.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn temperature_factor_below_37c_exceeds_one() {
        // Blood is more viscous at lower temperatures
        assert!(temperature_viscosity_factor(30.0) > 1.0);
        assert!(temperature_viscosity_factor(0.0) > 1.0);
    }

    #[test]
    fn temperature_factor_above_37c_less_than_one() {
        // Blood is less viscous at higher temperatures (post-cavitation heating)
        assert!(temperature_viscosity_factor(42.0) < 1.0);
        assert!(temperature_viscosity_factor(50.0) < 1.0);
    }

    #[test]
    fn temperature_factor_monotone_decreasing() {
        let f30 = temperature_viscosity_factor(30.0);
        let f37 = temperature_viscosity_factor(37.0);
        let f42 = temperature_viscosity_factor(42.0);
        assert!(f30 > f37 && f37 > f42);
    }

    // ── apparent_viscosity_at_temp ───────────────────────────────────────────

    #[test]
    fn apparent_viscosity_at_37c_matches_apparent_viscosity() {
        let blood = CassonBlood::<f64>::normal_blood();
        let mu_37 = blood.apparent_viscosity(100.0);
        let mu_temp =
            blood.apparent_viscosity_at_temp(100.0, ThermodynamicTemperature::from_base(310.15));
        assert!((mu_temp - mu_37).abs() < 1e-12);
    }

    #[test]
    fn apparent_viscosity_at_42c_is_lower() {
        let blood = CassonBlood::<f64>::normal_blood();
        let mu_37 =
            blood.apparent_viscosity_at_temp(100.0, ThermodynamicTemperature::from_base(310.15));
        let mu_42 =
            blood.apparent_viscosity_at_temp(100.0, ThermodynamicTemperature::from_base(315.15));
        assert!(mu_42 < mu_37);
    }

    #[test]
    fn apparent_viscosity_zero_temp_defaults_to_37c() {
        let blood = CassonBlood::<f64>::normal_blood();
        let mu_default =
            blood.apparent_viscosity_at_temp(100.0, ThermodynamicTemperature::from_base(0.0));
        let mu_37 =
            blood.apparent_viscosity_at_temp(100.0, ThermodynamicTemperature::from_base(310.15));
        assert!((mu_default - mu_37).abs() < 1e-12);
    }
}
