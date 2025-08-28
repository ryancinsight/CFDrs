//! Turbulence model trait definitions

use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};

/// Trait for turbulence models
pub trait TurbulenceModel<T: RealField + Copy> {
    /// Calculate turbulent viscosity
    fn turbulent_viscosity(&self, k: T, epsilon_or_omega: T, density: T) -> T;

    /// Calculate production term
    fn production_term(&self, velocity_gradient: &[[T; 2]; 2], turbulent_viscosity: T) -> T;

    /// Calculate dissipation/destruction term
    fn dissipation_term(&self, k: T, epsilon_or_omega: T) -> T;

    /// Update turbulence quantities
    fn update(
        &mut self,
        k: &mut [T],
        epsilon_or_omega: &mut [T],
        velocity: &[Vector2<T>],
        density: T,
        molecular_viscosity: T,
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()>;

    /// Get model name
    fn name(&self) -> &str;

    /// Check if model is valid for given Reynolds number
    fn is_valid_for_reynolds(&self, reynolds: T) -> bool;
}
