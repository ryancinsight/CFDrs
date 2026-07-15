//! Hephaestus-backed pointwise field arithmetic.

use super::GpuFieldOps;
use crate::error::{Error, Result};
use hephaestus_wgpu::{binary_elementwise, scalar_elementwise, AddOp, ComputeDevice, MulOp};

impl GpuFieldOps {
    /// Add two fields element-wise into `result`.
    ///
    /// # Errors
    /// Returns [`Error::DimensionMismatch`] when the host slices differ in
    /// length, or [`Error::GpuProvider`] when transfer or dispatch fails.
    pub fn add_fields(&self, left: &[f32], right: &[f32], result: &mut [f32]) -> Result<()> {
        validate_length(left.len(), right.len())?;
        validate_length(left.len(), result.len())?;
        if left.is_empty() {
            return Ok(());
        }

        let provider = self.context.provider();
        let left_buffer = provider.upload(left)?;
        let right_buffer = provider.upload(right)?;
        let output = binary_elementwise::<AddOp, f32>(provider, &left_buffer, &right_buffer)?;
        provider.download(&output, result)?;
        Ok(())
    }

    /// Multiply a field by `scalar` into `result`.
    ///
    /// # Errors
    /// Returns [`Error::DimensionMismatch`] when the host slices differ in
    /// length, or [`Error::GpuProvider`] when transfer or dispatch fails.
    pub fn multiply_field(&self, field: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        validate_length(field.len(), result.len())?;
        if field.is_empty() {
            return Ok(());
        }

        let provider = self.context.provider();
        let input = provider.upload(field)?;
        let output = scalar_elementwise::<MulOp, f32>(provider, &input, scalar)?;
        provider.download(&output, result)?;
        Ok(())
    }
}

fn validate_length(expected: usize, actual: usize) -> Result<()> {
    if expected == actual {
        Ok(())
    } else {
        Err(Error::DimensionMismatch { expected, actual })
    }
}
