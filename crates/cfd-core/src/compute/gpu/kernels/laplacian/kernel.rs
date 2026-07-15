//! Hephaestus-backed GPU implementation of the two-dimensional Laplacian.

use super::types::{BoundaryType, Laplacian2DUniforms};
use crate::compute::gpu::buffer::GpuBuffer;
use crate::compute::gpu::shaders::LAPLACIAN_2D_SHADER;
use crate::compute::gpu::GpuContext;
use crate::compute::traits::ComputeBuffer;
use crate::error::{Error, Result};
use hephaestus_wgpu::{
    ComputeDevice, DispatchGrid, MultiStorageKernel, WgslMultiStorageKernel, WgslStorageBinding,
    WgslStorageBindingLayout,
};
use std::sync::Arc;

const WORKGROUP: [usize; 3] = [8, 8, 1];

/// GPU kernel for the second-order two-dimensional Laplacian.
pub struct Laplacian2DKernel {
    context: Arc<GpuContext>,
    kernel: WgslMultiStorageKernel,
}

impl Laplacian2DKernel {
    /// Compile the Laplacian kernel through the Hephaestus provider.
    ///
    /// # Errors
    /// Returns [`Error::GpuProvider`] when the provider rejects the WGSL
    /// source or binding layout.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        let kernel = WgslMultiStorageKernel::new(
            context.provider(),
            "CFD Laplacian 2D",
            LAPLACIAN_2D_SHADER,
            "laplacian_2d",
            &[
                WgslStorageBindingLayout::read_only(0),
                WgslStorageBindingLayout::read_write(2),
            ],
            1,
        )?;
        Ok(Self { context, kernel })
    }

    /// Compute the Laplacian with homogeneous Dirichlet boundary conditions.
    ///
    /// # Errors
    /// Returns a typed configuration, dimension, or provider error when the
    /// grid contract, buffer lengths, transfer, or dispatch is invalid.
    pub fn execute(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        result: &mut [f32],
    ) -> Result<()> {
        self.execute_with_bc(field, nx, ny, dx, dy, BoundaryType::Dirichlet, result)
    }

    /// Compute the Laplacian with an explicit boundary condition.
    ///
    /// # Errors
    /// Returns a typed configuration, dimension, or provider error when the
    /// grid contract, buffer lengths, transfer, or dispatch is invalid.
    pub fn execute_with_bc(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        bc: BoundaryType,
        result: &mut [f32],
    ) -> Result<()> {
        let params = validate_contract(field.len(), result.len(), nx, ny, dx, dy, bc)?;
        let provider = self.context.provider();
        let input = provider.upload(field)?;
        let output = provider.alloc_zeroed(field.len())?;
        self.dispatch(&input, &output, &params, nx, ny)?;
        provider.download(&output, result)?;
        Ok(())
    }

    /// Dispatch directly over existing `f32` provider buffers.
    ///
    /// # Errors
    /// Returns a typed configuration, dimension, or provider error when the
    /// grid contract, buffer lengths, or dispatch is invalid.
    pub fn execute_on_gpu(
        &self,
        input: &GpuBuffer<f32>,
        output: &mut GpuBuffer<f32>,
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        bc: BoundaryType,
    ) -> Result<()> {
        let params = validate_contract(input.size(), output.size(), nx, ny, dx, dy, bc)?;
        self.dispatch(&input.buffer, &output.buffer, &params, nx, ny)
    }

    fn dispatch(
        &self,
        input: &hephaestus_wgpu::WgpuBuffer<f32>,
        output: &hephaestus_wgpu::WgpuBuffer<f32>,
        params: &Laplacian2DUniforms,
        nx: usize,
        ny: usize,
    ) -> Result<()> {
        let grid = DispatchGrid::covering_domain([nx, ny, 1], WORKGROUP)?;
        self.kernel.dispatch(
            self.context.provider(),
            [
                WgslStorageBinding::new(0, input),
                WgslStorageBinding::new(2, output),
            ],
            params,
            grid,
        )?;
        Ok(())
    }
}

fn validate_contract(
    input_len: usize,
    output_len: usize,
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    bc: BoundaryType,
) -> Result<Laplacian2DUniforms> {
    if nx < 2 || ny < 2 {
        return Err(Error::InvalidConfiguration(format!(
            "Laplacian grid axes must contain at least two points: nx={nx}, ny={ny}"
        )));
    }
    if !dx.is_finite() || dx <= 0.0 || !dy.is_finite() || dy <= 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "Laplacian spacing must be finite and positive: dx={dx}, dy={dy}"
        )));
    }
    let expected = nx.checked_mul(ny).ok_or_else(|| {
        Error::InvalidConfiguration(format!(
            "Laplacian grid size overflows usize: nx={nx}, ny={ny}"
        ))
    })?;
    if input_len != expected {
        return Err(Error::DimensionMismatch {
            expected,
            actual: input_len,
        });
    }
    if output_len != expected {
        return Err(Error::DimensionMismatch {
            expected,
            actual: output_len,
        });
    }
    let nx = u32::try_from(nx)
        .map_err(|_| Error::InvalidConfiguration(format!("Laplacian nx {nx} exceeds u32::MAX")))?;
    let ny = u32::try_from(ny)
        .map_err(|_| Error::InvalidConfiguration(format!("Laplacian ny {ny} exceeds u32::MAX")))?;

    Ok(Laplacian2DUniforms {
        dims_bc: [nx, ny, bc.as_u32(), 0],
        inv2: [dx.recip().powi(2), dy.recip().powi(2), 0.0, 0.0],
    })
}
