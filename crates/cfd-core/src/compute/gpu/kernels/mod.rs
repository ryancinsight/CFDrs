//! GPU compute kernels for CFD operations

use hephaestus_wgpu::wgpu;
pub mod advection;
pub mod diffusion;
pub mod laplacian;
pub mod pressure;
pub mod turbulence;
pub mod velocity;

pub use advection::{AdvectionConfig, GpuAdvectionKernel};
pub use diffusion::{DiffusionConfig, GpuDiffusionKernel};
pub use laplacian::Laplacian2DKernel;
pub use velocity::{GpuVelocityKernel, VelocityConfig};

use crate::compute::traits::{ComputeKernel, KernelParams};
use crate::error::{Error, Result};
use eunomia::RealField;

fn validate_field_len(expected: usize, actual: usize) -> Result<()> {
    if actual == expected {
        Ok(())
    } else {
        Err(Error::DimensionMismatch { expected, actual })
    }
}

fn validate_finite_field(name: &str, values: &[f32]) -> Result<()> {
    if let Some(index) = values.iter().position(|value| !value.is_finite()) {
        Err(Error::PhysicsViolation(format!(
            "{name} field is non-finite at index {index}"
        )))
    } else {
        Ok(())
    }
}

/// Base trait for GPU kernels
pub trait GpuKernel<T: RealField + Copy>: ComputeKernel<T> {
    /// Get WGSL shader code
    fn shader_code(&self) -> &str;

    /// Create compute pipeline
    ///
    /// # Errors
    ///
    /// Returns an error if the compute pipeline creation fails due to shader
    /// compilation issues, device limitations, or invalid pipeline configuration.
    fn create_pipeline(&self, device: &wgpu::Device) -> Result<wgpu::ComputePipeline>;

    /// Dispatch compute operation
    fn dispatch(&self, encoder: &mut wgpu::CommandEncoder, params: KernelParams);
}
