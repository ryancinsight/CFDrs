//! Validated turbulence grids and Hephaestus dispatch.

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
struct TurbulenceParams {
    dimensions: [u32; 4],
    scales: [f32; 4],
}

/// Validated two-dimensional Cartesian grid for turbulence operations.
#[derive(Debug, Clone, Copy)]
pub struct TurbulenceGrid {
    dimensions: [u32; 4],
    scales: [f32; 3],
    len: usize,
}

impl TurbulenceGrid {
    /// Construct a turbulence grid with centered-difference support.
    ///
    /// Both dimensions must contain at least three points and both spacings
    /// must be finite and positive.
    ///
    /// # Errors
    /// Returns [`Error::InvalidConfiguration`] when an invariant is violated,
    /// derived scales are not representable, the element count overflows, or
    /// a dimension exceeds `u32::MAX`.
    pub fn new(dimensions: [usize; 2], spacing: [f32; 2]) -> Result<Self> {
        if dimensions.iter().any(|extent| *extent < 3) {
            return Err(Error::InvalidConfiguration(format!(
                "turbulence grid requires both axes to contain at least three points; got {dimensions:?}"
            )));
        }
        if spacing
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(Error::InvalidConfiguration(format!(
                "turbulence spacing must be finite and positive; got {spacing:?}"
            )));
        }
        let len = dimensions[0].checked_mul(dimensions[1]).ok_or_else(|| {
            Error::InvalidConfiguration(format!(
                "turbulence grid element count overflows usize: {dimensions:?}"
            ))
        })?;
        let axis = |value| {
            u32::try_from(value).map_err(|_| {
                Error::InvalidConfiguration(format!(
                    "turbulence grid axis {value} exceeds u32::MAX"
                ))
            })
        };
        let dimensions = [axis(dimensions[0])?, axis(dimensions[1])?, 1, 0];
        let filter_width = (spacing[0] * spacing[1]).sqrt();
        let scales = [0.5 / spacing[0], 0.5 / spacing[1], filter_width];
        if scales.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidConfiguration(
                "turbulence grid scales are not representable as finite f32 values".to_string(),
            ));
        }
        Ok(Self {
            dimensions,
            scales,
            len,
        })
    }

    /// Number of scalar values required by this grid.
    #[must_use]
    pub const fn element_count(self) -> usize {
        self.len
    }

    fn params(self, constant: f32) -> Result<TurbulenceParams> {
        if !constant.is_finite() || constant < 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "turbulence model constant must be finite and nonnegative; got {constant}"
            )));
        }
        Ok(TurbulenceParams {
            dimensions: self.dimensions,
            scales: [self.scales[0], self.scales[1], self.scales[2], constant],
        })
    }
}

/// Compiled turbulence operation family sharing one provider context.
pub(crate) struct TurbulenceKernels {
    context: Arc<GpuContext>,
    smagorinsky: WgslMultiStorageKernel,
    des_grid_scale: WgslMultiStorageKernel,
    wall_distance: WgslMultiStorageKernel,
}

impl TurbulenceKernels {
    pub(crate) fn new(context: Arc<GpuContext>) -> Result<Self> {
        let smagorinsky = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Smagorinsky SGS Viscosity",
            include_str!("smagorinsky.wgsl"),
            "smagorinsky_viscosity",
            &[
                WgslStorageBindingLayout::read_only(1),
                WgslStorageBindingLayout::read_only(2),
                WgslStorageBindingLayout::read_write(3),
            ],
            0,
        )?;
        let output_layout = [WgslStorageBindingLayout::read_write(1)];
        let des_grid_scale = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD DES Grid Length Scale",
            include_str!("des_grid_scale.wgsl"),
            "des_grid_scale",
            &output_layout,
            0,
        )?;
        let wall_distance = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Rectangular Wall Distance",
            include_str!("wall_distance.wgsl"),
            "wall_distance",
            &output_layout,
            0,
        )?;
        Ok(Self {
            context,
            smagorinsky,
            des_grid_scale,
            wall_distance,
        })
    }

    pub(crate) fn smagorinsky(
        &self,
        velocity_x: &[f32],
        velocity_y: &[f32],
        grid: TurbulenceGrid,
        constant: f32,
        output: &mut [f32],
    ) -> Result<()> {
        for actual in [velocity_x.len(), velocity_y.len(), output.len()] {
            validate_field_len(grid.len, actual)?;
        }
        validate_finite_field("turbulence velocity x", velocity_x)?;
        validate_finite_field("turbulence velocity y", velocity_y)?;
        let params = grid.params(constant)?;
        let provider = self.context.provider();
        let velocity_x = provider.upload(velocity_x)?;
        let velocity_y = provider.upload(velocity_y)?;
        let result = provider.alloc_zeroed(grid.len)?;
        self.smagorinsky.dispatch(
            provider,
            [
                WgslStorageBinding::new(1, &velocity_x),
                WgslStorageBinding::new(2, &velocity_y),
                WgslStorageBinding::new(3, &result),
            ],
            &params,
            dispatch_grid(grid)?,
        )?;
        provider.download(&result, output)?;
        Ok(())
    }

    pub(crate) fn des_grid_scale(
        &self,
        grid: TurbulenceGrid,
        constant: f32,
        output: &mut [f32],
    ) -> Result<()> {
        self.output_only(&self.des_grid_scale, grid, constant, output)
    }

    pub(crate) fn wall_distance(&self, grid: TurbulenceGrid, output: &mut [f32]) -> Result<()> {
        self.output_only(&self.wall_distance, grid, 0.0, output)
    }

    fn output_only(
        &self,
        kernel: &WgslMultiStorageKernel,
        grid: TurbulenceGrid,
        constant: f32,
        output: &mut [f32],
    ) -> Result<()> {
        validate_field_len(grid.len, output.len())?;
        let params = grid.params(constant)?;
        let provider = self.context.provider();
        let result = provider.alloc_zeroed(grid.len)?;
        kernel.dispatch(
            provider,
            [WgslStorageBinding::new(1, &result)],
            &params,
            dispatch_grid(grid)?,
        )?;
        provider.download(&result, output)?;
        Ok(())
    }
}

fn dispatch_grid(grid: TurbulenceGrid) -> Result<DispatchGrid> {
    Ok(DispatchGrid::covering_domain(
        [
            usize::try_from(grid.dimensions[0]).expect("invariant: u32 fits usize"),
            usize::try_from(grid.dimensions[1]).expect("invariant: u32 fits usize"),
            1,
        ],
        WORKGROUP,
    )?)
}
