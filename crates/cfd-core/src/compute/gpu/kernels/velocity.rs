//! GPU velocity correction kernel implementation

use super::GpuKernel;
use crate::compute::traits::{ComputeBackend, ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use std::marker::PhantomData;

/// GPU velocity correction kernel for SIMPLE algorithm
pub struct GpuVelocityKernel<T: RealField + Copy> {
    shader_module: Option<wgpu::ShaderModule>,
    _phantom: PhantomData<T>,
}

impl<T: RealField + Copy> GpuVelocityKernel<T> {
    /// Shader source code
    const SHADER_SOURCE: &'static str = include_str!("velocity.wgsl");

    pub fn new() -> Self {
        Self {
            shader_module: None,
            _phantom: PhantomData,
        }
    }
}

impl<T: RealField + Copy> ComputeKernel<T> for GpuVelocityKernel<T> {
    fn name(&self) -> &str {
        "GPU Velocity Correction (SIMPLE)"
    }

    fn execute(&self, _input: &[T], _output: &mut [T], _params: KernelParams) -> Result<()> {
        // GPU execution handled through GpuKernel trait
        Ok(())
    }

    fn complexity(&self, size: usize) -> usize {
        use crate::compute::gpu::constants::complexity;
        size * complexity::VELOCITY_UPDATE
    }

    fn supports_backend(&self, backend: &ComputeBackend) -> bool {
        matches!(backend, ComputeBackend::Gpu | ComputeBackend::Hybrid)
    }
}

impl<T: RealField + Copy> GpuKernel<T> for GpuVelocityKernel<T> {
    fn shader_code(&self) -> &str {
        Self::SHADER_SOURCE
    }

    fn create_pipeline(&self, device: &wgpu::Device) -> Result<wgpu::ComputePipeline> {
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Velocity Shader"),
            source: wgpu::ShaderSource::Wgsl(self.shader_code().into()),
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Velocity Pipeline Layout"),
            bind_group_layouts: &[],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Velocity Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader_module,
            entry_point: "velocity_correction",
        });

        Ok(pipeline)
    }

    fn dispatch(&self, encoder: &mut wgpu::CommandEncoder, params: KernelParams) {
        let (nx, ny, nz) = params.domain_params.grid_dims;
        let work_group_size = params.work_group_size;

        // Calculate dispatch dimensions
        let dispatch_x = (nx + work_group_size - 1) / work_group_size;
        let dispatch_y = (ny + work_group_size - 1) / work_group_size;
        let dispatch_z = (nz + work_group_size - 1) / work_group_size;

        let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("Velocity Pass"),
            timestamp_writes: None,
        });

        compute_pass.dispatch_workgroups(dispatch_x as u32, dispatch_y as u32, dispatch_z as u32);
    }
}

impl<T: RealField + Copy> std::fmt::Debug for GpuVelocityKernel<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuVelocityKernel").finish()
    }
}
