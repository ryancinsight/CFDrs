//! Validated SIMPLE velocity operations and provider dispatch.

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

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct VelocityParams {
    dimensions: [u32; 4],
    gradient_scale: [f32; 4],
}

/// Validated grid and physical coefficients for SIMPLE velocity operations.
#[derive(Debug, Clone, Copy)]
pub struct VelocityConfig {
    params: VelocityParams,
    len: usize,
}

impl VelocityConfig {
    /// Construct a validated velocity-operation configuration.
    ///
    /// Every axis must contain at least three points for centered differences.
    /// Spacings, timestep, and density must be finite and strictly positive.
    ///
    /// # Errors
    /// Returns [`Error::InvalidConfiguration`] when an invariant is violated,
    /// the derived correction factor is not representable, the element count
    /// overflows `usize`, or an axis exceeds `u32::MAX`.
    pub fn new(dimensions: [usize; 3], spacing: [f32; 3], dt: f32, density: f32) -> Result<Self> {
        if dimensions.iter().any(|extent| *extent < 3) {
            return Err(Error::InvalidConfiguration(format!(
                "velocity grid requires every axis to contain at least three points; got {dimensions:?}"
            )));
        }
        if spacing
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(Error::InvalidConfiguration(format!(
                "velocity spacing must be finite and positive; got {spacing:?}"
            )));
        }
        if !dt.is_finite() || dt <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "velocity timestep must be finite and positive; got {dt}"
            )));
        }
        if !density.is_finite() || density <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "density must be finite and positive; got {density}"
            )));
        }

        let correction_factor = dt / density;
        let gradient_scale = spacing.map(|value| 0.5 / value);
        if !correction_factor.is_finite() || gradient_scale.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidConfiguration(
                "velocity coefficients are not representable as finite f32 values".to_string(),
            ));
        }

        let [nx, ny, nz] = dimensions;
        let len = nx
            .checked_mul(ny)
            .and_then(|plane| plane.checked_mul(nz))
            .ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "velocity grid element count overflows usize: {dimensions:?}"
                ))
            })?;
        let axis = |value| {
            u32::try_from(value).map_err(|_| {
                Error::InvalidConfiguration(format!("velocity grid axis {value} exceeds u32::MAX"))
            })
        };
        let [nx, ny, nz] = [axis(nx)?, axis(ny)?, axis(nz)?];

        Ok(Self {
            params: VelocityParams {
                dimensions: [nx, ny, nz, 0],
                gradient_scale: [
                    gradient_scale[0],
                    gradient_scale[1],
                    gradient_scale[2],
                    correction_factor,
                ],
            },
            len,
        })
    }

    /// Number of values required for every scalar or vector component field.
    #[must_use]
    pub const fn element_count(self) -> usize {
        self.len
    }
}

/// Hephaestus-backed SIMPLE velocity correction and divergence-source kernels.
pub struct GpuVelocityKernel {
    context: Arc<GpuContext>,
    correction: WgslMultiStorageKernel,
    divergence_source: WgslMultiStorageKernel,
}

impl GpuVelocityKernel {
    /// Compile the velocity shaders through Hephaestus.
    ///
    /// # Errors
    /// Returns [`Error::GpuProvider`] when either shader or binding contract is
    /// rejected by the provider.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        let correction = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD SIMPLE Velocity Correction",
            include_str!("correction.wgsl"),
            "velocity_correction",
            &[
                WgslStorageBindingLayout::read_only(1),
                WgslStorageBindingLayout::read_only(2),
                WgslStorageBindingLayout::read_only(3),
                WgslStorageBindingLayout::read_only(4),
                WgslStorageBindingLayout::read_write(5),
                WgslStorageBindingLayout::read_write(6),
                WgslStorageBindingLayout::read_write(7),
            ],
            0,
        )?;
        let divergence_source = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD SIMPLE Divergence Source",
            include_str!("divergence_source.wgsl"),
            "divergence_source",
            &[
                WgslStorageBindingLayout::read_only(1),
                WgslStorageBindingLayout::read_only(2),
                WgslStorageBindingLayout::read_only(3),
                WgslStorageBindingLayout::read_write(4),
            ],
            0,
        )?;
        Ok(Self {
            context,
            correction,
            divergence_source,
        })
    }

    /// Correct an intermediate velocity using the centered pressure gradient.
    ///
    /// Outer cells are set to zero to impose the kernel's no-slip boundary
    /// contract.
    ///
    /// # Errors
    /// Returns a typed dimension error for mismatched fields, a physical
    /// invariant error for non-finite input, or a provider error for
    /// transfer/dispatch failure.
    pub fn correct(
        &self,
        velocity: [&[f32]; 3],
        pressure: &[f32],
        config: VelocityConfig,
        output: [&mut [f32]; 3],
    ) -> Result<()> {
        let [velocity_x, velocity_y, velocity_z] = velocity;
        let [output_x, output_y, output_z] = output;
        for actual in [
            velocity_x.len(),
            velocity_y.len(),
            velocity_z.len(),
            pressure.len(),
            output_x.len(),
            output_y.len(),
            output_z.len(),
        ] {
            validate_field_len(config.len, actual)?;
        }
        validate_finite_field("velocity x", velocity_x)?;
        validate_finite_field("velocity y", velocity_y)?;
        validate_finite_field("velocity z", velocity_z)?;
        validate_finite_field("pressure", pressure)?;

        let provider = self.context.provider();
        let velocity_x = provider.upload(velocity_x)?;
        let velocity_y = provider.upload(velocity_y)?;
        let velocity_z = provider.upload(velocity_z)?;
        let pressure = provider.upload(pressure)?;
        let result_x = provider.alloc_zeroed(config.len)?;
        let result_y = provider.alloc_zeroed(config.len)?;
        let result_z = provider.alloc_zeroed(config.len)?;
        self.correction.dispatch(
            provider,
            [
                WgslStorageBinding::new(1, &velocity_x),
                WgslStorageBinding::new(2, &velocity_y),
                WgslStorageBinding::new(3, &velocity_z),
                WgslStorageBinding::new(4, &pressure),
                WgslStorageBinding::new(5, &result_x),
                WgslStorageBinding::new(6, &result_y),
                WgslStorageBinding::new(7, &result_z),
            ],
            &config.params,
            dispatch_grid(config.params.dimensions)?,
        )?;
        provider.download(&result_x, output_x)?;
        provider.download(&result_y, output_y)?;
        provider.download(&result_z, output_z)?;
        Ok(())
    }

    /// Assemble `(density / dt) * divergence(velocity)` for a pressure solve.
    ///
    /// Outer cells are set to zero.
    ///
    /// # Errors
    /// Returns a typed dimension error for mismatched fields, a physical
    /// invariant error for non-finite input, or a provider error for
    /// transfer/dispatch failure.
    pub fn divergence_source(
        &self,
        velocity_x: &[f32],
        velocity_y: &[f32],
        velocity_z: &[f32],
        config: VelocityConfig,
        output: &mut [f32],
    ) -> Result<()> {
        for actual in [
            velocity_x.len(),
            velocity_y.len(),
            velocity_z.len(),
            output.len(),
        ] {
            validate_field_len(config.len, actual)?;
        }
        validate_finite_field("velocity x", velocity_x)?;
        validate_finite_field("velocity y", velocity_y)?;
        validate_finite_field("velocity z", velocity_z)?;

        let provider = self.context.provider();
        let velocity_x = provider.upload(velocity_x)?;
        let velocity_y = provider.upload(velocity_y)?;
        let velocity_z = provider.upload(velocity_z)?;
        let result = provider.alloc_zeroed(config.len)?;
        self.divergence_source.dispatch(
            provider,
            [
                WgslStorageBinding::new(1, &velocity_x),
                WgslStorageBinding::new(2, &velocity_y),
                WgslStorageBinding::new(3, &velocity_z),
                WgslStorageBinding::new(4, &result),
            ],
            &config.params,
            dispatch_grid(config.params.dimensions)?,
        )?;
        provider.download(&result, output)?;
        Ok(())
    }
}

fn dispatch_grid(dimensions: [u32; 4]) -> Result<DispatchGrid> {
    let [nx, ny, nz, _] = dimensions;
    Ok(DispatchGrid::covering_domain(
        [
            usize::try_from(nx).expect("invariant: u32 fits usize"),
            usize::try_from(ny).expect("invariant: u32 fits usize"),
            usize::try_from(nz).expect("invariant: u32 fits usize"),
        ],
        WORKGROUP,
    )?)
}

impl std::fmt::Debug for GpuVelocityKernel {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("GpuVelocityKernel")
            .finish_non_exhaustive()
    }
}
