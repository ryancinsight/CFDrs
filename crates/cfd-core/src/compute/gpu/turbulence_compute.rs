//! GPU turbulence operations through the Hephaestus provider.

use super::kernels::turbulence::{TurbulenceGrid, TurbulenceKernels};
use super::GpuContext;
use crate::error::Result;
use std::sync::Arc;

/// GPU turbulence operation facade sharing one provider context.
pub struct GpuTurbulenceCompute {
    context: Arc<GpuContext>,
    kernels: TurbulenceKernels,
}

impl GpuTurbulenceCompute {
    /// Compile every turbulence operation through Hephaestus.
    ///
    /// # Errors
    /// Returns a typed error when provider acquisition or kernel compilation
    /// fails.
    pub fn new() -> Result<Self> {
        let context = Arc::new(GpuContext::create()?);
        let kernels = TurbulenceKernels::new(context.clone())?;
        Ok(Self { context, kernels })
    }

    /// Compute Smagorinsky subgrid viscosity into caller-owned storage.
    ///
    /// # Errors
    /// Returns a typed contract or provider error.
    pub fn compute_smagorinsky_sgs(
        &self,
        velocity_x: &[f32],
        velocity_y: &[f32],
        grid: TurbulenceGrid,
        constant: f32,
        output: &mut [f32],
    ) -> Result<()> {
        self.kernels
            .smagorinsky(velocity_x, velocity_y, grid, constant, output)
    }

    /// Compute the grid-cutoff branch `constant * sqrt(dx * dy)` of the DES
    /// length scale into caller-owned storage.
    ///
    /// # Errors
    /// Returns a typed contract or provider error.
    pub fn compute_des_length_scale(
        &self,
        grid: TurbulenceGrid,
        constant: f32,
        output: &mut [f32],
    ) -> Result<()> {
        self.kernels.des_grid_scale(grid, constant, output)
    }

    /// Compute exact distance to the nearest rectangular domain wall.
    ///
    /// # Errors
    /// Returns a typed contract or provider error.
    pub fn compute_wall_distance(&self, grid: TurbulenceGrid, output: &mut [f32]) -> Result<()> {
        self.kernels.wall_distance(grid, output)
    }

    /// Provider context used by this facade.
    #[must_use]
    pub fn context(&self) -> &GpuContext {
        &self.context
    }
}

impl std::fmt::Debug for GpuTurbulenceCompute {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("GpuTurbulenceCompute")
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests;
