//! Hephaestus-backed GPU implementation of the two-dimensional Laplacian.
//!
//! The WGSL source and dispatch machinery now live in the Hephaestus provider;
//! this module is a thin typed consumer that validates the CFD grid contract
//! and forwards to `hephaestus_wgpu::Laplacian2DKernel`.

use super::{BoundaryCondition, LaplacianPolarity};
use crate::compute::gpu::buffer::GpuBuffer;
use crate::compute::gpu::GpuContext;
use crate::compute::traits::ComputeBuffer;
use crate::error::{Error, Result};
use aequitas::systems::si::quantities::Length;
use hephaestus_wgpu::{
    ComputeDevice, Laplacian2DKernel as HephaestusLaplacian2DKernel, Laplacian2DParams,
};
use std::sync::Arc;

/// GPU kernel for the second-order two-dimensional Laplacian.
pub struct Laplacian2DKernel {
    context: Arc<GpuContext>,
    kernel: HephaestusLaplacian2DKernel,
}

impl Laplacian2DKernel {
    /// Compile the Laplacian kernel through the Hephaestus provider.
    ///
    /// # Errors
    /// Returns [`Error::GpuProvider`] when the provider rejects the WGSL
    /// source or binding layout.
    pub fn new(context: Arc<GpuContext>) -> Result<Self> {
        let kernel = HephaestusLaplacian2DKernel::new(context.provider())?;
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
        dx: Length<f32>,
        dy: Length<f32>,
        result: &mut [f32],
    ) -> Result<()> {
        self.execute_with_bc(field, nx, ny, dx, dy, BoundaryCondition::Dirichlet, result)
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
        dx: Length<f32>,
        dy: Length<f32>,
        bc: BoundaryCondition,
        result: &mut [f32],
    ) -> Result<()> {
        let params = validate_contract(
            field.len(),
            result.len(),
            nx,
            ny,
            dx,
            dy,
            bc,
            LaplacianPolarity::Laplacian,
        )?;
        let provider = self.context.provider();
        let input = provider.upload(field)?;
        let output = provider.alloc_zeroed(field.len())?;
        self.kernel.dispatch(provider, &input, &output, &params)?;
        provider.download(&output, result)?;
        Ok(())
    }

    /// Dispatch the selected Laplacian polarity over existing `f32` provider
    /// buffers.
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
        dx: Length<f32>,
        dy: Length<f32>,
        bc: BoundaryCondition,
        polarity: LaplacianPolarity,
    ) -> Result<()> {
        let params = validate_contract(input.size(), output.size(), nx, ny, dx, dy, bc, polarity)?;
        self.kernel.dispatch(
            self.context.provider(),
            &input.buffer,
            &output.buffer,
            &params,
        )?;
        Ok(())
    }
}

fn validate_contract(
    input_len: usize,
    output_len: usize,
    nx: usize,
    ny: usize,
    dx: Length<f32>,
    dy: Length<f32>,
    bc: BoundaryCondition,
    polarity: LaplacianPolarity,
) -> Result<Laplacian2DParams> {
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

    Laplacian2DParams::new(nx, ny, dx, dy, bc, polarity)
        .map_err(|error| Error::InvalidConfiguration(error.to_string()))
}
