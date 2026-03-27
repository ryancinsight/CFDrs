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

/// Carreau-Yasuda blood rheology model for wide shear rate range.
pub mod carreau_yasuda;
/// Casson blood rheology model with yield stress.
pub mod casson;
/// Blood physical constants at 37°C (body temperature).
pub mod constants;
/// Cross blood model (simpler alternative to Carreau-Yasuda).
pub mod cross;
/// Fåhræus-Lindqvist effect for microvascular blood flow.
pub mod fahraeus_lindqvist;

pub use carreau_yasuda::CarreauYasudaBlood;
pub use casson::{temperature_viscosity_factor, CassonBlood};
pub use cross::CrossBlood;
pub use fahraeus_lindqvist::FahraeuasLindqvist;

use nalgebra::RealField;
use num_traits::FromPrimitive;

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

impl<T: RealField + Copy + FromPrimitive> BloodModel<T> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::fluid::traits::{Fluid as FluidTrait, NonNewtonianFluid};

    #[test]
    fn test_fluid_trait_implementations() {
        let casson = CassonBlood::<f64>::normal_blood();
        let carreau = CarreauYasudaBlood::<f64>::normal_blood();
        let cross_blood = CrossBlood::<f64>::normal_blood();

        // All should implement FluidTrait
        let casson_state = casson.properties_at(310.0, 101325.0).unwrap();
        let carreau_state = carreau.properties_at(310.0, 101325.0).unwrap();
        let cross_state = cross_blood.properties_at(310.0, 101325.0).unwrap();

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

    #[test]
    fn test_model_comparison_at_intermediate_shear() {
        let casson = CassonBlood::<f64>::normal_blood();
        let carreau = CarreauYasudaBlood::<f64>::normal_blood();
        let cross_blood = CrossBlood::<f64>::normal_blood();

        // At γ̇ = 100 s⁻¹ (typical arterial condition), all models should agree within factor of 2
        let gamma = 100.0;
        let mu_casson = casson.apparent_viscosity(gamma);
        let mu_carreau = carreau.apparent_viscosity(gamma);
        let mu_cross = cross_blood.apparent_viscosity(gamma);

        // All should be in reasonable range (3-10 mPa·s)
        for (name, mu) in [
            ("Casson", mu_casson),
            ("Carreau", mu_carreau),
            ("Cross", mu_cross),
        ] {
            assert!(
                mu > 0.003 && mu < 0.010,
                "{name} viscosity at 100 s⁻¹ = {mu} Pa·s out of expected range"
            );
        }
    }
}
