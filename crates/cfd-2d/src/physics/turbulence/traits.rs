//! Turbulence model trait definitions

use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};

/// Trait for turbulence models
pub trait TurbulenceModel<T: RealField + Copy> {
    /// Calculate turbulent viscosity
    /// 
    /// # Arguments
    /// * `k` - Turbulent kinetic energy
    /// * `epsilon_or_omega` - Dissipation rate (ε) or specific dissipation rate (ω)
    /// * `density` - Fluid density
    fn turbulent_viscosity(&self, k: T, epsilon_or_omega: T, density: T) -> T;

    /// Calculate turbulent viscosity with strain rate limiter (for SST model)
    /// 
    /// This method implements the Bradshaw assumption per Menter (1994):
    /// νt = a1*k / max(a1*ω, S*F2) where S is strain rate magnitude
    /// 
    /// # Arguments
    /// * `k` - Turbulent kinetic energy
    /// * `omega` - Specific dissipation rate
    /// * `density` - Fluid density
    /// * `strain_rate_magnitude` - Magnitude of strain rate tensor S = sqrt(2*Sij*Sij)
    /// * `f2` - Blending function F2 value (1 near walls, 0 in freestream)
    fn turbulent_viscosity_with_limiter(
        &self,
        k: T,
        omega: T,
        density: T,
        strain_rate_magnitude: T,
        f2: T,
    ) -> T {
        // Default implementation without limiter (for k-ε and SA models)
        self.turbulent_viscosity(k, omega, density)
    }

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
