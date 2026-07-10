//! Validated pressure operations and provider dispatch.

use crate::compute::gpu::kernels::{dispatch_grid_3d, validate_field_len, validate_finite_field};
use crate::compute::gpu::GpuContext;
use crate::error::{Error, Result};
use bytemuck::{Pod, Zeroable};
use hephaestus_wgpu::{
    ComputeDevice, MultiStorageKernel, WgslMultiStorageKernel, WgslStorageBinding,
    WgslStorageBindingLayout,
};
use std::sync::Arc;

const WORKGROUP: [usize; 3] = [8, 8, 1];

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct PressureParams {
    dimensions: [u32; 4],
    inverse_spacing_squared_relaxation: [f32; 4],
}

/// Validated grid and relaxation contract for pressure operations.
#[derive(Debug, Clone, Copy)]
pub struct PressureConfig {
    params: PressureParams,
    len: usize,
}

impl PressureConfig {
    /// Construct a validated pressure-operation configuration.
    ///
    /// Every axis must contain at least three points. Spacings must be finite
    /// and positive. Weighted-Jacobi relaxation must lie in `(0, 1]`.
    ///
    /// # Errors
    /// Returns [`Error::InvalidConfiguration`] when an invariant is violated,
    /// derived coefficients are not representable, the element count
    /// overflows `usize`, or an axis exceeds `u32::MAX`.
    pub fn new(dimensions: [usize; 3], spacing: [f32; 3], relaxation: f32) -> Result<Self> {
        if dimensions.iter().any(|extent| *extent < 3) {
            return Err(Error::InvalidConfiguration(format!(
                "pressure grid requires every axis to contain at least three points; got {dimensions:?}"
            )));
        }
        if spacing
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(Error::InvalidConfiguration(format!(
                "pressure spacing must be finite and positive; got {spacing:?}"
            )));
        }
        if !relaxation.is_finite() || relaxation <= 0.0 || relaxation > 1.0 {
            return Err(Error::InvalidConfiguration(format!(
                "weighted-Jacobi relaxation must be finite and in (0, 1]; got {relaxation}"
            )));
        }

        let inverse_spacing_squared = spacing.map(|value| value.recip().powi(2));
        if inverse_spacing_squared
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(Error::InvalidConfiguration(
                "pressure inverse spacing squared is not representable as finite f32 values"
                    .to_string(),
            ));
        }

        let [nx, ny, nz] = dimensions;
        let len = nx
            .checked_mul(ny)
            .and_then(|plane| plane.checked_mul(nz))
            .ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "pressure grid element count overflows usize: {dimensions:?}"
                ))
            })?;
        let axis = |value| {
            u32::try_from(value).map_err(|_| {
                Error::InvalidConfiguration(format!("pressure grid axis {value} exceeds u32::MAX"))
            })
        };
        let [nx, ny, nz] = [axis(nx)?, axis(ny)?, axis(nz)?];

        Ok(Self {
            params: PressureParams {
                dimensions: [nx, ny, nz, 0],
                inverse_spacing_squared_relaxation: [
                    inverse_spacing_squared[0],
                    inverse_spacing_squared[1],
                    inverse_spacing_squared[2],
                    relaxation,
                ],
            },
            len,
        })
    }

    /// Number of values required for every pressure/source field.
    #[must_use]
    pub const fn element_count(self) -> usize {
        self.len
    }
}

/// Hephaestus-backed weighted-Jacobi pressure and residual kernels.
pub struct GpuPressureKernel {
    context: Arc<GpuContext>,
    iteration: WgslMultiStorageKernel,
    residual: WgslMultiStorageKernel,
}

impl GpuPressureKernel {
    /// Compile the pressure shaders through Hephaestus.
    ///
    /// # Errors
    /// Returns [`Error::GpuProvider`] when either shader or binding contract is
    /// rejected by the provider.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        let layouts = [
            WgslStorageBindingLayout::read_only(1),
            WgslStorageBindingLayout::read_only(2),
            WgslStorageBindingLayout::read_write(3),
        ];
        let iteration = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Weighted Jacobi Pressure",
            include_str!("iteration.wgsl"),
            "pressure_iteration",
            &layouts,
            0,
        )?;
        let residual = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Pressure Residual",
            include_str!("residual.wgsl"),
            "pressure_residual",
            &layouts,
            0,
        )?;
        Ok(Self {
            context,
            iteration,
            residual,
        })
    }

    /// Perform one weighted-Jacobi pressure iteration.
    ///
    /// Boundary cells copy the nearest cell that is interior on every boundary
    /// axis, enforcing the homogeneous Neumann contract at edges and corners.
    ///
    /// # Errors
    /// Returns a typed dimension error for mismatched fields, a physical
    /// invariant error for non-finite input, or a provider error for
    /// transfer/dispatch failure.
    pub fn iterate(
        &self,
        pressure: &[f32],
        source: &[f32],
        config: PressureConfig,
        output: &mut [f32],
    ) -> Result<()> {
        self.execute(&self.iteration, pressure, source, config, output)
    }

    /// Evaluate `abs(Laplacian(pressure) - source)` at every interior cell.
    ///
    /// Boundary residuals are zero.
    ///
    /// # Errors
    /// Returns a typed dimension error for mismatched fields, a physical
    /// invariant error for non-finite input, or a provider error for
    /// transfer/dispatch failure.
    pub fn residual(
        &self,
        pressure: &[f32],
        source: &[f32],
        config: PressureConfig,
        output: &mut [f32],
    ) -> Result<()> {
        self.execute(&self.residual, pressure, source, config, output)
    }

    fn execute(
        &self,
        kernel: &WgslMultiStorageKernel,
        pressure: &[f32],
        source: &[f32],
        config: PressureConfig,
        output: &mut [f32],
    ) -> Result<()> {
        for actual in [pressure.len(), source.len(), output.len()] {
            validate_field_len(config.len, actual)?;
        }
        validate_finite_field("pressure", pressure)?;
        validate_finite_field("pressure source", source)?;

        let provider = self.context.provider();
        let pressure = provider.upload(pressure)?;
        let source = provider.upload(source)?;
        let result = provider.alloc_zeroed(config.len)?;
        kernel.dispatch(
            provider,
            [
                WgslStorageBinding::new(1, &pressure),
                WgslStorageBinding::new(2, &source),
                WgslStorageBinding::new(3, &result),
            ],
            &config.params,
            dispatch_grid_3d(config.params.dimensions, WORKGROUP)?,
        )?;
        provider.download(&result, output)?;
        Ok(())
    }
}

impl std::fmt::Debug for GpuPressureKernel {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("GpuPressureKernel")
            .finish_non_exhaustive()
    }
}
