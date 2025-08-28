//! GPU compute kernels for CFD operations

pub mod advection;
pub mod diffusion;
pub mod pressure;

use crate::compute::traits::{ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use std::sync::Arc;

/// Base trait for GPU kernels
pub trait GpuKernel<T: RealField + Copy>: ComputeKernel<T> {
    /// Get WGSL shader code
    fn shader_code(&self) -> &str;

    /// Create compute pipeline
    fn create_pipeline(&self, device: &wgpu::Device) -> Result<wgpu::ComputePipeline>;

    /// Dispatch compute operation
    fn dispatch(&self, encoder: &mut wgpu::CommandEncoder, params: KernelParams);
}
