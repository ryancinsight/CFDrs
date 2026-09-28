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
pub use casson::{CassonBlood, temperature_viscosity_factor};
pub use cross::CrossBlood;
pub use fahraeus_lindqvist::FahraeuasLindqvist;

use eunomia::FloatElement;
use eunomia::RealField;

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

impl<T: RealField + Copy + FloatElement> BloodModel<T> {
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
mod tests;
