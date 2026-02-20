//! Blood rheology models for hemodynamic CFD simulations
//!
//! Blood is a shear-thinning, non-Newtonian fluid composed of plasma, red blood cells,
//! white blood cells, and platelets. This module provides mathematically rigorous
//! implementations of blood viscosity models validated against published literature.
//!
//! # Rheological Models
//!
//! ## Casson Model
//! Standard model for blood flow in larger vessels. The constitutive equation is:
//! ```text
//! √τ = √τ_y + √(μ_∞ · γ̇)
//! ```
//! Yields apparent viscosity:
//! ```text
//! μ_app = (√τ_y / √γ̇ + √μ_∞)²
//! ```
//!
//! ## Carreau-Yasuda Model
//! Accurate across full shear rate range (0.01 - 1000 s⁻¹):
//! ```text
//! μ(γ̇) = μ_∞ + (μ_0 - μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
//! ```
//!
//! ## Fåhræus-Lindqvist Effect
//! Apparent viscosity reduction in microvessels (D < 300 μm) due to RBC migration:
//! ```text
//! μ_rel = μ_app / μ_plasma = f(D, H_t)
//! ```
//!
//! # References
//! - Merrill, E.W. et al. (1969) "Pressure-flow relations of human blood in hollow fibers"
//! - Cho, Y.I., Kensey, K.R. (1991) "Effects of the non-Newtonian viscosity of blood"
//! - Pries, A.R. et al. (1992) "Blood viscosity in tube flow: dependence on diameter"
//! - Chien, S. (1970) "Shear dependence of effective cell volume as a determinant of blood viscosity"
//! - Fung, Y.C. (1993) "Biomechanics: Mechanical Properties of Living Tissues"

use super::traits::{Fluid as FluidTrait, FluidState, NonNewtonianFluid};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

// ============================================================================
// Physical Constants for Blood
// ============================================================================

/// Blood physical constants at 37°C (body temperature)
///
/// Reference: Fung (1993), Chien (1970)
pub mod constants {
    /// Plasma viscosity at 37°C [Pa·s]
    /// Reference: Merrill et al. (1969)
    pub const PLASMA_VISCOSITY_37C: f64 = 0.00122;

    /// Blood density [kg/m³]
    /// Reference: Fung (1993)
    pub const BLOOD_DENSITY: f64 = 1060.0;

    /// Zero-shear viscosity for normal blood (H_t = 45%) [Pa·s]
    /// Reference: Cho & Kensey (1991)
    pub const ZERO_SHEAR_VISCOSITY: f64 = 0.056;

    /// Infinite-shear viscosity for normal blood [Pa·s]
    /// Reference: Cho & Kensey (1991)
    pub const INFINITE_SHEAR_VISCOSITY: f64 = 0.00345;

    /// Yield stress for normal blood (H_t = 45%) [Pa]
    /// Reference: Merrill et al. (1969)
    pub const YIELD_STRESS: f64 = 0.0056;

    /// Casson viscosity parameter √μ_∞ for normal blood [√(Pa·s)]
    /// Reference: Merrill et al. (1969)
    pub const CASSON_VISCOSITY_SQRT: f64 = 0.0588; // √0.00345 ≈ 0.0587

    /// Carreau-Yasuda relaxation time λ [s]
    /// Reference: Cho & Kensey (1991)
    pub const CARREAU_LAMBDA: f64 = 3.313;

    /// Carreau-Yasuda power-law index n [-]
    /// Reference: Cho & Kensey (1991)
    pub const CARREAU_N: f64 = 0.3568;

    /// Carreau-Yasuda transition parameter a [-]
    /// Reference: Cho & Kensey (1991)
    pub const CARREAU_A: f64 = 2.0;

    /// Blood specific heat capacity [J/(kg·K)]
    /// Reference: Fung (1993)
    pub const BLOOD_SPECIFIC_HEAT: f64 = 3770.0;

    /// Blood thermal conductivity [W/(m·K)]
    /// Reference: Fung (1993)
    pub const BLOOD_THERMAL_CONDUCTIVITY: f64 = 0.52;

    /// Speed of sound in blood [m/s]
    /// Reference: Fung (1993)
    pub const BLOOD_SPEED_OF_SOUND: f64 = 1570.0;

    /// Normal hematocrit (volume fraction of RBCs)
    pub const NORMAL_HEMATOCRIT: f64 = 0.45;

    /// Critical vessel diameter for Fåhræus-Lindqvist effect [m]
    pub const FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER: f64 = 300e-6;
}

// ============================================================================
// Casson Blood Model
// ============================================================================

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
pub struct CassonBlood<T: RealField + Copy> {
    // Note: String fields removed to allow Copy trait
    // /// Fluid name (removed for Copy compatibility)
    // pub name: String,
    /// Blood density [kg/m³]
    pub density: T,
    /// Yield stress τ_y [Pa]
    /// Literature value: 0.0056 Pa for normal blood (Merrill 1969)
    pub yield_stress: T,
    /// Infinite-shear (Casson) viscosity μ_∞ [Pa·s]
    /// Literature value: 0.00345 Pa·s (Merrill 1969)
    pub infinite_shear_viscosity: T,
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
    /// Regularization parameter for low shear rates [1/s]
    /// Prevents singularity at γ̇ → 0
    pub regularization_shear_rate: T,
}

impl<T: RealField + FromPrimitive + Copy> CassonBlood<T> {
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
            // name field removed for Copy compatibility
            density: T::from_f64(constants::BLOOD_DENSITY).unwrap(),
            yield_stress: T::from_f64(constants::YIELD_STRESS).unwrap(),
            infinite_shear_viscosity: T::from_f64(constants::INFINITE_SHEAR_VISCOSITY).unwrap(),
            hematocrit: T::from_f64(constants::NORMAL_HEMATOCRIT).unwrap(),
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap(),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap(),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap(),
            reference_shear_rate: T::from_f64(100.0).unwrap(),
            regularization_shear_rate: T::from_f64(0.001).unwrap(),
        }
    }

    /// Create Casson blood model with custom parameters
    pub fn new(density: T, yield_stress: T, infinite_shear_viscosity: T, hematocrit: T) -> Self {
        Self {
            // name field removed for Copy compatibility
            density,
            yield_stress,
            infinite_shear_viscosity,
            hematocrit,
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap(),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap(),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap(),
            reference_shear_rate: T::from_f64(100.0).unwrap(),
            regularization_shear_rate: T::from_f64(0.001).unwrap(),
        }
    }

    /// Create blood model with hematocrit-dependent properties
    ///
    /// Yields stress and viscosity scale with hematocrit using empirical correlations:
    /// - τ_y(H_t) = τ_y,ref · (H_t / H_t,ref)^3  (Chien 1970)
    /// - μ_∞(H_t) = μ_plasma · exp(k · H_t / (1 - H_t))  (Quemada 1978)
    pub fn with_hematocrit(hematocrit: T) -> Self {
        let h_ref = T::from_f64(constants::NORMAL_HEMATOCRIT).unwrap();
        let tau_ref = T::from_f64(constants::YIELD_STRESS).unwrap();
        let mu_plasma = T::from_f64(constants::PLASMA_VISCOSITY_37C).unwrap();

        // Yield stress scaling with hematocrit (Chien 1970)
        let ratio = hematocrit / h_ref;
        let yield_stress = tau_ref * ratio * ratio * ratio;

        // Infinite-shear viscosity from Quemada model
        let k = T::from_f64(2.5).unwrap(); // Intrinsic viscosity coefficient
        let one = T::one();
        let infinite_shear_viscosity = mu_plasma * (k * hematocrit / (one - hematocrit)).exp();

        Self {
            // name field removed for Copy compatibility
            density: T::from_f64(constants::BLOOD_DENSITY).unwrap(),
            yield_stress,
            infinite_shear_viscosity,
            hematocrit,
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap(),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap(),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap(),
            reference_shear_rate: T::from_f64(100.0).unwrap(),
            regularization_shear_rate: T::from_f64(0.001).unwrap(),
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
        let gamma_eff = if shear_rate < self.regularization_shear_rate {
            self.regularization_shear_rate
        } else {
            shear_rate
        };

        // Casson model: μ_app = (√τ_y/√γ̇ + √μ_∞)²
        let sqrt_tau_y = self.yield_stress.sqrt();
        let sqrt_mu_inf = self.infinite_shear_viscosity.sqrt();
        let sqrt_gamma = gamma_eff.sqrt();

        let casson_sqrt = sqrt_tau_y / sqrt_gamma + sqrt_mu_inf;
        casson_sqrt * casson_sqrt
    }

    /// Calculate shear stress from shear rate
    ///
    /// τ = (√τ_y + √(μ_∞ · γ̇))²
    pub fn shear_stress(&self, shear_rate: T) -> T {
        self.apparent_viscosity(shear_rate) * shear_rate
    }

    /// Validate model parameters
    ///
    /// # Errors
    /// Returns error if parameters are physically invalid
    pub fn validate(&self) -> Result<(), Error> {
        if self.density <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        if self.yield_stress < T::zero() {
            return Err(Error::InvalidInput(
                "Yield stress must be non-negative".to_string(),
            ));
        }
        if self.infinite_shear_viscosity <= T::zero() {
            return Err(Error::InvalidInput(
                "Infinite-shear viscosity must be positive".to_string(),
            ));
        }
        if self.hematocrit < T::zero() || self.hematocrit > T::one() {
            return Err(Error::InvalidInput(
                "Hematocrit must be between 0 and 1".to_string(),
            ));
        }
        Ok(())
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for CassonBlood<T> {
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
    ) -> Result<T, Error> {
        Ok(self.apparent_viscosity(shear_rate))
    }
}

impl<T: RealField + FromPrimitive + Copy> NonNewtonianFluid<T> for CassonBlood<T> {
    fn apparent_viscosity(&self, shear_rate: T) -> T {
        CassonBlood::apparent_viscosity(self, shear_rate)
    }

    fn has_yield_stress(&self) -> bool {
        true
    }

    fn yield_stress(&self) -> Option<T> {
        Some(self.yield_stress)
    }
}

// ============================================================================
// Carreau-Yasuda Blood Model
// ============================================================================

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
    // Note: String fields removed to allow Copy trait
    // /// Fluid name (removed for Copy compatibility)
    // pub name: String,
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
            // name field removed for Copy compatibility
            density: T::from_f64(constants::BLOOD_DENSITY).unwrap(),
            zero_shear_viscosity: T::from_f64(constants::ZERO_SHEAR_VISCOSITY).unwrap(),
            infinite_shear_viscosity: T::from_f64(constants::INFINITE_SHEAR_VISCOSITY).unwrap(),
            relaxation_time: T::from_f64(constants::CARREAU_LAMBDA).unwrap(),
            power_law_index: T::from_f64(constants::CARREAU_N).unwrap(),
            transition_parameter: T::from_f64(constants::CARREAU_A).unwrap(),
            hematocrit: T::from_f64(constants::NORMAL_HEMATOCRIT).unwrap(),
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap(),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap(),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap(),
            reference_shear_rate: T::from_f64(100.0).unwrap(),
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
            // name field removed for Copy compatibility
            density,
            zero_shear_viscosity,
            infinite_shear_viscosity,
            relaxation_time,
            power_law_index,
            transition_parameter,
            hematocrit,
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap(),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap(),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap(),
            reference_shear_rate: T::from_f64(100.0).unwrap(),
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

    fn name(&self) -> &str {
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

// ============================================================================
// Cross Blood Model
// ============================================================================

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
            density: T::from_f64(constants::BLOOD_DENSITY).unwrap(),
            zero_shear_viscosity: T::from_f64(constants::ZERO_SHEAR_VISCOSITY).unwrap(),
            infinite_shear_viscosity: T::from_f64(constants::INFINITE_SHEAR_VISCOSITY).unwrap(),
            time_constant: T::from_f64(1.007).unwrap(), // Fitted for blood
            rate_index: T::from_f64(1.028).unwrap(),    // Fitted for blood
            hematocrit: T::from_f64(constants::NORMAL_HEMATOCRIT).unwrap(),
            specific_heat: T::from_f64(constants::BLOOD_SPECIFIC_HEAT).unwrap(),
            thermal_conductivity: T::from_f64(constants::BLOOD_THERMAL_CONDUCTIVITY).unwrap(),
            speed_of_sound: T::from_f64(constants::BLOOD_SPEED_OF_SOUND).unwrap(),
            reference_shear_rate: T::from_f64(100.0).unwrap(),
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

// ============================================================================
// Fåhræus-Lindqvist Effect
// ============================================================================

/// Fåhræus-Lindqvist effect for microvascular blood flow
///
/// In microvessels (D < 300 μm), blood exhibits apparent viscosity reduction due to:
/// 1. Axial migration of RBCs (cell-free layer near wall)
/// 2. Fåhræus effect (reduction in tube hematocrit vs feed hematocrit)
///
/// # Physical Behavior
/// The apparent viscosity DECREASES as vessel diameter decreases from 300 μm down to
/// approximately 7 μm (minimum viscosity), then increases again for very small vessels
/// where the RBC size becomes comparable to vessel diameter.
///
/// # Empirical Correlation
/// Uses simplified Pries et al. (1994) correlation:
/// ```text
/// μ_rel(D) = 1 + (μ_bulk - 1) · f(D)
///
/// where f(D) accounts for the diameter dependence of the Fåhræus-Lindqvist effect
/// ```
///
/// # Reference
/// Pries, A.R., Neuhaus, D., Gaehtgens, P. (1992) "Blood viscosity in tube flow:
/// dependence on diameter and hematocrit"
#[derive(Debug, Clone, Copy)]
pub struct FahraeuasLindqvist<T: RealField + Copy> {
    /// Vessel diameter [m]
    pub diameter: T,
    /// Hematocrit (volume fraction of RBCs) [-]
    pub hematocrit: T,
    /// Plasma viscosity [Pa·s]
    pub plasma_viscosity: T,
}

impl<T: RealField + FromPrimitive + Copy> FahraeuasLindqvist<T> {
    /// Create new Fåhræus-Lindqvist calculator
    pub fn new(diameter: T, hematocrit: T) -> Self {
        Self {
            diameter,
            hematocrit,
            plasma_viscosity: T::from_f64(constants::PLASMA_VISCOSITY_37C).unwrap(),
        }
    }

    /// Check if Fåhræus-Lindqvist effect is significant
    ///
    /// Effect is significant for D < 300 μm
    pub fn is_significant(&self) -> bool {
        self.diameter < T::from_f64(constants::FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER).unwrap()
    }

    /// Calculate relative apparent viscosity using Pries et al. (1992) simplified model
    ///
    /// Returns μ_rel = μ_app / μ_plasma
    ///
    /// # Physical Model
    /// For a vessel of diameter D and hematocrit H_t:
    /// - Bulk viscosity: μ_bulk/μ_plasma = 1 + 2.2·H_t (large vessel limit)
    /// - F-L reduction factor: decreases for smaller vessels down to ~10 μm
    ///
    /// The relative viscosity follows:
    /// μ_rel = 1 + (μ_bulk/μ_plasma - 1) · (1 - 1.7·exp(-D/22.5))
    ///
    /// where D is in micrometers.
    pub fn relative_viscosity(&self) -> T {
        let one = T::one();
        let d_um = self.diameter * T::from_f64(1e6).unwrap(); // Convert to μm

        // Bulk relative viscosity (large vessel asymptote)
        // μ_bulk/μ_plasma ≈ 1 + 2.2·H_t for H_t = 0.45 gives ~2.0
        let mu_bulk_rel = one + T::from_f64(2.2).unwrap() * self.hematocrit;

        // Fåhræus-Lindqvist reduction factor
        // This factor goes from 0 (very small D) to 1 (large D)
        // f(D) = 1 - 1.7·exp(-D/22.5)
        let fl_factor =
            one - T::from_f64(1.7).unwrap() * (-d_um / T::from_f64(22.5).unwrap()).exp();

        // Clamp factor to valid range [0, 1]
        let fl_factor = if fl_factor < T::zero() {
            T::zero()
        } else if fl_factor > one {
            one
        } else {
            fl_factor
        };

        // Final relative viscosity
        // μ_rel = 1 + (μ_bulk_rel - 1) · fl_factor
        one + (mu_bulk_rel - one) * fl_factor
    }

    /// Calculate apparent viscosity in microvessel [Pa·s]
    pub fn apparent_viscosity(&self) -> T {
        self.plasma_viscosity * self.relative_viscosity()
    }

    /// Calculate tube hematocrit from feed hematocrit (Fåhræus effect)
    ///
    /// H_tube / H_feed = empirical correlation
    pub fn tube_hematocrit(&self) -> T {
        let d_um = self.diameter * T::from_f64(1e6).unwrap();
        let one = T::one();

        // Empirical correlation for tube hematocrit reduction
        // Valid for D > 10 μm
        if d_um < T::from_f64(10.0).unwrap() {
            return self.hematocrit * T::from_f64(0.5).unwrap(); // Approximate minimum
        }

        // Pries et al. empirical fit
        let reduction_factor =
            one - T::from_f64(1.7).unwrap() * (-d_um / T::from_f64(40.0).unwrap()).exp();
        self.hematocrit * reduction_factor
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // ========================================================================
    // Casson Model Tests
    // ========================================================================

    #[test]
    fn test_casson_normal_blood_creation() {
        let blood = CassonBlood::<f64>::normal_blood();
        assert_eq!(blood.density, constants::BLOOD_DENSITY);
        assert_eq!(blood.yield_stress, constants::YIELD_STRESS);
        assert_eq!(
            blood.infinite_shear_viscosity,
            constants::INFINITE_SHEAR_VISCOSITY
        );
        assert_eq!(blood.hematocrit, constants::NORMAL_HEMATOCRIT);
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
        let mut invalid = blood.clone();
        invalid.density = -1.0;
        assert!(invalid.validate().is_err());
    }

    // ========================================================================
    // Carreau-Yasuda Model Tests
    // ========================================================================

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

    // ========================================================================
    // Cross Model Tests
    // ========================================================================

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

    // ========================================================================
    // Fåhræus-Lindqvist Tests
    // ========================================================================

    #[test]
    fn test_fahraeus_lindqvist_significance() {
        // Large vessel - effect not significant
        let large = FahraeuasLindqvist::<f64>::new(1e-3, 0.45); // 1 mm
        assert!(!large.is_significant());

        // Small vessel - effect significant
        let small = FahraeuasLindqvist::<f64>::new(50e-6, 0.45); // 50 μm
        assert!(small.is_significant());
    }

    #[test]
    fn test_fahraeus_lindqvist_viscosity_reduction() {
        // In smaller vessels, apparent viscosity should be lower
        let d_100 = FahraeuasLindqvist::<f64>::new(100e-6, 0.45);
        let d_50 = FahraeuasLindqvist::<f64>::new(50e-6, 0.45);

        let mu_100 = d_100.apparent_viscosity();
        let mu_50 = d_50.apparent_viscosity();

        // Smaller vessel should have lower apparent viscosity
        assert!(
            mu_50 < mu_100,
            "μ(50 μm) = {} should be < μ(100 μm) = {}",
            mu_50,
            mu_100
        );
    }

    #[test]
    fn test_fahraeus_lindqvist_tube_hematocrit() {
        let fl = FahraeuasLindqvist::<f64>::new(50e-6, 0.45);
        let h_tube = fl.tube_hematocrit();

        // Tube hematocrit should be less than feed due to Fåhræus effect
        assert!(
            h_tube < 0.45,
            "Tube hematocrit {} should be less than feed {}",
            h_tube,
            0.45
        );
    }

    // ========================================================================
    // Trait Implementation Tests
    // ========================================================================

    #[test]
    fn test_fluid_trait_implementations() {
        let casson = CassonBlood::<f64>::normal_blood();
        let carreau = CarreauYasudaBlood::<f64>::normal_blood();
        let cross = CrossBlood::<f64>::normal_blood();

        // All should implement FluidTrait
        let casson_state = casson.properties_at(310.0, 101325.0).unwrap();
        let carreau_state = carreau.properties_at(310.0, 101325.0).unwrap();
        let cross_state = cross.properties_at(310.0, 101325.0).unwrap();

        assert_eq!(casson_state.density, constants::BLOOD_DENSITY);
        assert_eq!(carreau_state.density, constants::BLOOD_DENSITY);
        assert_eq!(cross_state.density, constants::BLOOD_DENSITY);

        // All have positive viscosity
        assert!(casson_state.dynamic_viscosity > 0.0);
        assert!(carreau_state.dynamic_viscosity > 0.0);
        assert!(cross_state.dynamic_viscosity > 0.0);
    }

    #[test]
    fn test_non_newtonian_trait() {
        let casson = CassonBlood::<f64>::normal_blood();
        let carreau = CarreauYasudaBlood::<f64>::normal_blood();

        // Casson has yield stress
        assert!(casson.has_yield_stress());
        assert!(casson.yield_stress().is_some());

        // Carreau-Yasuda does not have yield stress
        assert!(!carreau.has_yield_stress());
        assert!(carreau.yield_stress().is_none());
    }

    // ========================================================================
    // Comparison Tests
    // ========================================================================

    #[test]
    fn test_model_comparison_at_intermediate_shear() {
        let casson = CassonBlood::<f64>::normal_blood();
        let carreau = CarreauYasudaBlood::<f64>::normal_blood();
        let cross = CrossBlood::<f64>::normal_blood();

        // At γ̇ = 100 s⁻¹ (typical arterial condition), all models should agree within factor of 2
        let gamma = 100.0;
        let mu_casson = casson.apparent_viscosity(gamma);
        let mu_carreau = carreau.apparent_viscosity(gamma);
        let mu_cross = cross.apparent_viscosity(gamma);

        // All should be in reasonable range (3-10 mPa·s)
        for (name, mu) in [
            ("Casson", mu_casson),
            ("Carreau", mu_carreau),
            ("Cross", mu_cross),
        ] {
            assert!(
                mu > 0.003 && mu < 0.010,
                "{} viscosity at 100 s⁻¹ = {} Pa·s out of expected range",
                name,
                mu
            );
        }
    }
}

// ============================================================================
// Blood Model Dispatch Enum
// ============================================================================

/// Dispatch enum for blood rheology models.
///
/// The single authoritative selector over Casson, Carreau-Yasuda, and Newtonian
/// (constant μ) blood models. All solver crates MUST import this type rather than
/// defining their own dispatch enum.
#[derive(Debug, Clone)]
pub enum BloodModel<T: RealField + Copy> {
    /// Casson model with yield stress (suitable for large vessels, D > 500 µm)
    Casson(CassonBlood<T>),
    /// Carreau-Yasuda model (full shear-rate range, 0.01–1000 s⁻¹)
    CarreauYasuda(CarreauYasudaBlood<T>),
    /// Newtonian approximation with constant dynamic viscosity [Pa·s]
    Newtonian(T),
}

impl<T: RealField + Copy + num_traits::FromPrimitive> BloodModel<T> {
    /// Compute apparent dynamic viscosity at `shear_rate` [s⁻¹].
    #[must_use]
    pub fn viscosity(&self, shear_rate: T) -> T {
        match self {
            BloodModel::Casson(m) => m.apparent_viscosity(shear_rate),
            BloodModel::CarreauYasuda(m) => m.apparent_viscosity(shear_rate),
            BloodModel::Newtonian(mu) => *mu,
        }
    }

    /// Returns `true` for Newtonian blood (constant viscosity).
    pub fn is_newtonian(&self) -> bool {
        matches!(self, BloodModel::Newtonian(_))
    }
}

