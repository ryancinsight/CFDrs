//! Validated advection configuration and provider dispatch.

use crate::compute::gpu::kernels::validate_field_len;
use crate::compute::gpu::GpuContext;
use crate::error::{Error, Result};
use bytemuck::{Pod, Zeroable};
use hephaestus_wgpu::{
    ComputeDevice, DispatchGrid, MultiStorageKernel, WgslMultiStorageKernel, WgslStorageBinding,
    WgslStorageBindingLayout,
};
use std::sync::Arc;

const WORKGROUP: [usize; 3] = [8, 8, 1];

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct AdvectionParams {
    dimensions: [u32; 4],
    spacing_step: [f32; 4],
}

/// Validated grid and timestep for two-dimensional upwind advection on z-planes.
#[derive(Debug, Clone, Copy)]
pub struct AdvectionConfig {
    params: AdvectionParams,
    len: usize,
}

impl AdvectionConfig {
    /// Construct a validated advection configuration.
    ///
    /// `nx` and `ny` must contain at least two points, `nz` must be nonzero,
    /// spacings must be finite and positive, and `dt` must be finite and
    /// nonnegative.
    ///
    /// # Errors
    /// Returns [`Error::InvalidConfiguration`] when any invariant is violated
    /// or the grid element count overflows `usize`/the dimensions exceed
    /// `u32::MAX`.
    pub fn new(dimensions: [usize; 3], spacing: [f32; 3], dt: f32) -> Result<Self> {
        let [nx, ny, nz] = dimensions;
        if nx < 2 || ny < 2 || nz == 0 {
            return Err(Error::InvalidConfiguration(format!(
                "advection grid requires nx>=2, ny>=2, nz>=1; got {dimensions:?}"
            )));
        }
        if spacing
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(Error::InvalidConfiguration(format!(
                "advection spacing must be finite and positive; got {spacing:?}"
            )));
        }
        if !dt.is_finite() || dt < 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "advection timestep must be finite and nonnegative; got {dt}"
            )));
        }
        let len = nx
            .checked_mul(ny)
            .and_then(|plane| plane.checked_mul(nz))
            .ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "advection grid element count overflows usize: {dimensions:?}"
                ))
            })?;
        let axis = |value| {
            u32::try_from(value).map_err(|_| {
                Error::InvalidConfiguration(format!("advection grid axis {value} exceeds u32::MAX"))
            })
        };
        let [nx, ny, nz] = [axis(nx)?, axis(ny)?, axis(nz)?];

        Ok(Self {
            params: AdvectionParams {
                dimensions: [nx, ny, nz, 0],
                spacing_step: [spacing[0], spacing[1], spacing[2], dt],
            },
            len,
        })
    }

    /// Number of scalar/velocity values required by this grid.
    #[must_use]
    pub const fn element_count(self) -> usize {
        self.len
    }
}

/// Hephaestus-backed first-order upwind advection kernel.
pub struct GpuAdvectionKernel {
    context: Arc<GpuContext>,
    kernel: WgslMultiStorageKernel,
}

impl GpuAdvectionKernel {
    /// Compile the advection shader through Hephaestus.
    ///
    /// # Errors
    /// Returns [`Error::GpuProvider`] when the shader or binding contract is
    /// rejected by the provider.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        let kernel = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Upwind Advection",
            include_str!("kernel.wgsl"),
            "advection_upwind",
            &[
                WgslStorageBindingLayout::read_only(1),
                WgslStorageBindingLayout::read_only(2),
                WgslStorageBindingLayout::read_only(3),
                WgslStorageBindingLayout::read_write(4),
            ],
            0,
        )?;
        Ok(Self { context, kernel })
    }

    /// Advance `scalar` by one validated upwind step.
    ///
    /// # Errors
    /// Returns a typed dimension error for mismatched fields, a physical
    /// invariant error for non-finite values or CFL violation, or a provider
    /// error for transfer/dispatch failure.
    pub fn execute(
        &self,
        scalar: &[f32],
        velocity_x: &[f32],
        velocity_y: &[f32],
        config: AdvectionConfig,
        output: &mut [f32],
    ) -> Result<()> {
        validate_field_len(config.len, scalar.len())?;
        validate_field_len(config.len, velocity_x.len())?;
        validate_field_len(config.len, velocity_y.len())?;
        validate_field_len(config.len, output.len())?;
        validate_values(scalar, velocity_x, velocity_y, config.params.spacing_step)?;

        let provider = self.context.provider();
        let scalar = provider.upload(scalar)?;
        let velocity_x = provider.upload(velocity_x)?;
        let velocity_y = provider.upload(velocity_y)?;
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
                WgslStorageBinding::new(1, &scalar),
                WgslStorageBinding::new(2, &velocity_x),
                WgslStorageBinding::new(3, &velocity_y),
                WgslStorageBinding::new(4, &result),
            ],
            &config.params,
            grid,
        )?;
        provider.download(&result, output)?;
        Ok(())
    }
}

fn validate_values(
    scalar: &[f32],
    velocity_x: &[f32],
    velocity_y: &[f32],
    spacing_step: [f32; 4],
) -> Result<()> {
    if let Some(index) = scalar.iter().position(|value| !value.is_finite()) {
        return Err(Error::PhysicsViolation(format!(
            "advection scalar field is non-finite at index {index}"
        )));
    }
    let [dx, dy, _, dt] = spacing_step;
    let mut max_cfl = 0.0_f32;
    for (index, (&vx, &vy)) in velocity_x.iter().zip(velocity_y).enumerate() {
        if !vx.is_finite() || !vy.is_finite() {
            return Err(Error::PhysicsViolation(format!(
                "advection velocity is non-finite at index {index}: ({vx}, {vy})"
            )));
        }
        max_cfl = max_cfl.max(dt * (vx.abs() / dx + vy.abs() / dy));
    }
    if max_cfl > 1.0 {
        return Err(Error::PhysicsViolation(format!(
            "advection CFL {max_cfl} exceeds the explicit upwind limit 1"
        )));
    }
    Ok(())
}

impl std::fmt::Debug for GpuAdvectionKernel {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("GpuAdvectionKernel")
            .finish_non_exhaustive()
    }
}
