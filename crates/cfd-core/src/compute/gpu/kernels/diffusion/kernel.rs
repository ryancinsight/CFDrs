//! Validated diffusion configuration and provider dispatch.

use crate::compute::gpu::kernels::{validate_field_len, validate_finite_field};
use crate::compute::gpu::GpuContext;
use crate::error::{Error, Result};
use bytemuck::{Pod, Zeroable};
use hephaestus_wgpu::{
    ComputeDevice, DispatchGrid, MultiStorageKernel, WgslMultiStorageKernel, WgslStorageBinding,
    WgslStorageBindingLayout,
};
use std::sync::Arc;

const WORKGROUP: [usize; 3] = [8, 8, 1];
const EXPLICIT_STABILITY_LIMIT: f32 = 0.5;

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct DiffusionParams {
    dimensions: [u32; 4],
    inverse_spacing_squared_step: [f32; 4],
}

/// Validated grid and coefficients for three-dimensional explicit diffusion.
#[derive(Debug, Clone, Copy)]
pub struct DiffusionConfig {
    params: DiffusionParams,
    len: usize,
}

impl DiffusionConfig {
    /// Construct a validated diffusion configuration.
    ///
    /// Every axis must contain at least three points for the centered stencil.
    /// Spacings must be finite and positive; the timestep and diffusivity must
    /// be finite and nonnegative. The explicit update must satisfy the
    /// three-dimensional von Neumann stability bound
    /// `diffusivity * dt * sum(1 / spacing_i^2) <= 1/2`.
    ///
    /// # Errors
    /// Returns [`Error::InvalidConfiguration`] when an invariant is violated,
    /// the grid element count overflows `usize`, or an axis exceeds `u32::MAX`.
    pub fn new(
        dimensions: [usize; 3],
        spacing: [f32; 3],
        dt: f32,
        diffusivity: f32,
    ) -> Result<Self> {
        if dimensions.iter().any(|extent| *extent < 3) {
            return Err(Error::InvalidConfiguration(format!(
                "diffusion grid requires every axis to contain at least three points; got {dimensions:?}"
            )));
        }
        if spacing
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(Error::InvalidConfiguration(format!(
                "diffusion spacing must be finite and positive; got {spacing:?}"
            )));
        }
        if !dt.is_finite() || dt < 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "diffusion timestep must be finite and nonnegative; got {dt}"
            )));
        }
        if !diffusivity.is_finite() || diffusivity < 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "diffusivity must be finite and nonnegative; got {diffusivity}"
            )));
        }

        let inverse_spacing_squared = spacing.map(|value| value.recip().powi(2));
        let diffusion_step = dt * diffusivity;
        let stability_number =
            diffusion_step * inverse_spacing_squared.iter().copied().sum::<f32>();
        if !stability_number.is_finite() || stability_number > EXPLICIT_STABILITY_LIMIT {
            return Err(Error::InvalidConfiguration(format!(
                "diffusion stability number {stability_number} exceeds the explicit limit {EXPLICIT_STABILITY_LIMIT}"
            )));
        }

        let [nx, ny, nz] = dimensions;
        let len = nx
            .checked_mul(ny)
            .and_then(|plane| plane.checked_mul(nz))
            .ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "diffusion grid element count overflows usize: {dimensions:?}"
                ))
            })?;
        let axis = |value| {
            u32::try_from(value).map_err(|_| {
                Error::InvalidConfiguration(format!("diffusion grid axis {value} exceeds u32::MAX"))
            })
        };
        let [nx, ny, nz] = [axis(nx)?, axis(ny)?, axis(nz)?];

        Ok(Self {
            params: DiffusionParams {
                dimensions: [nx, ny, nz, 0],
                inverse_spacing_squared_step: [
                    inverse_spacing_squared[0],
                    inverse_spacing_squared[1],
                    inverse_spacing_squared[2],
                    diffusion_step,
                ],
            },
            len,
        })
    }

    /// Number of scalar values required by this grid.
    #[must_use]
    pub const fn element_count(self) -> usize {
        self.len
    }
}

/// Hephaestus-backed explicit central-difference diffusion kernel.
pub struct GpuDiffusionKernel {
    context: Arc<GpuContext>,
    kernel: WgslMultiStorageKernel,
}

impl GpuDiffusionKernel {
    /// Compile the diffusion shader through Hephaestus.
    ///
    /// # Errors
    /// Returns [`Error::GpuProvider`] when the shader or binding contract is
    /// rejected by the provider.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        let kernel = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Central Diffusion",
            include_str!("kernel.wgsl"),
            "diffusion_central",
            &[
                WgslStorageBindingLayout::read_only(1),
                WgslStorageBindingLayout::read_write(2),
            ],
            0,
        )?;
        Ok(Self { context, kernel })
    }

    /// Advance `input` by one validated explicit diffusion step.
    ///
    /// # Errors
    /// Returns a typed dimension error for mismatched fields, a physical
    /// invariant error for non-finite input, or a provider error for
    /// transfer/dispatch failure.
    pub fn execute(
        &self,
        input: &[f32],
        config: DiffusionConfig,
        output: &mut [f32],
    ) -> Result<()> {
        validate_field_len(config.len, input.len())?;
        validate_field_len(config.len, output.len())?;
        validate_finite_field("diffusion input", input)?;

        let provider = self.context.provider();
        let input = provider.upload(input)?;
        let result = provider.alloc_zeroed(config.len)?;
        let [nx, ny, nz, _] = config.params.dimensions;
        let grid = DispatchGrid::covering_domain(
            [
                usize::try_from(nx).expect("invariant: u32 fits usize"),
                usize::try_from(ny).expect("invariant: u32 fits usize"),
                usize::try_from(nz).expect("invariant: u32 fits usize"),
            ],
            WORKGROUP,
        )?;
        self.kernel.dispatch(
            provider,
            [
                WgslStorageBinding::new(1, &input),
                WgslStorageBinding::new(2, &result),
            ],
            &config.params,
            grid,
        )?;
        provider.download(&result, output)?;
        Ok(())
    }
}

impl std::fmt::Debug for GpuDiffusionKernel {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("GpuDiffusionKernel")
            .finish_non_exhaustive()
    }
}
