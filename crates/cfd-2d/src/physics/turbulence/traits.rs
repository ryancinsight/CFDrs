//! Turbulence model trait definitions

use nalgebra::{DMatrix, RealField, Vector2};
use num_traits::ToPrimitive;

/// Trait for RANS turbulence models (k-ε, k-ω, SA)
pub trait TurbulenceModel<T: RealField + Copy + ToPrimitive> {
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
        // Parameters intentionally unused in default implementation
        let _ = (strain_rate_magnitude, f2);
        self.turbulent_viscosity(k, omega, density)
    }

    /// Calculate production term
    ///
    /// # Arguments
    /// * `velocity_gradient` - Gradient of the velocity field
    /// * `turbulent_viscosity` - Eddy viscosity (νt)
    /// * `turbulence_variable` - Primary turbulence variable (e.g., k for k-ε/k-ω, ν̃ for SA)
    /// * `wall_distance` - Distance to nearest wall (needed for SA)
    /// * `molecular_viscosity` - Molecular viscosity (needed for SA)
    fn production_term(
        &self,
        velocity_gradient: &[[T; 2]; 2],
        turbulent_viscosity: T,
        turbulence_variable: T,
        wall_distance: T,
        molecular_viscosity: T,
    ) -> T;

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
    ) -> cfd_core::error::Result<()>;

    /// Get model name
    fn name(&self) -> &str;

    /// Check if model is valid for given Reynolds number
    fn is_valid_for_reynolds(&self, reynolds: T) -> bool;
}

/// Trait for LES/DES turbulence models (Smagorinsky LES, Detached Eddy Simulation)
pub trait LESTurbulenceModel {
    /// Update turbulence model with flow field data
    ///
    /// # Arguments
    /// * `velocity_u` - x-velocity field
    /// * `velocity_v` - y-velocity field
    /// * `pressure` - pressure field
    /// * `density` - fluid density
    /// * `viscosity` - molecular viscosity
    /// * `dt` - time step
    /// * `dx`, `dy` - grid spacing
    fn update(
        &mut self,
        velocity_u: &DMatrix<f64>,
        velocity_v: &DMatrix<f64>,
        pressure: &DMatrix<f64>,
        density: f64,
        viscosity: f64,
        dt: f64,
        dx: f64,
        dy: f64,
    ) -> cfd_core::error::Result<()>;

    /// Get total viscosity (molecular + turbulent) at a grid point
    fn get_viscosity(&self, i: usize, j: usize) -> f64;

    /// Get turbulent viscosity field
    fn get_turbulent_viscosity_field(&self) -> &DMatrix<f64>;

    /// Get turbulent kinetic energy at a grid point (may return 0 for LES)
    fn get_turbulent_kinetic_energy(&self, i: usize, j: usize) -> f64;

    /// Get dissipation rate at a grid point
    fn get_dissipation_rate(&self, i: usize, j: usize) -> f64;

    /// Update boundary conditions
    fn boundary_condition_update(
        &mut self,
        boundary_manager: &super::boundary_conditions::TurbulenceBoundaryManager<f64>,
    ) -> cfd_core::error::Result<()>;

    /// Get model name
    fn get_model_name(&self) -> &str;

    /// Get model constants as (name, value) pairs
    fn get_model_constants(&self) -> Vec<(&str, f64)>;
}
