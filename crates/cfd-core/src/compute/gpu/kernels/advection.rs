//! GPU advection kernel implementation

use super::GpuKernel;
use crate::compute::traits::{ComputeBackend, ComputeKernel, KernelParams};
use crate::error::Result;
use nalgebra::RealField;
use std::marker::PhantomData;

/// GPU advection kernel using upwind scheme
pub struct GpuAdvectionKernel<T: RealField + Copy> {
    shader_module: Option<wgpu::ShaderModule>,
    _phantom: PhantomData<T>,
}

impl<T: RealField + Copy> Default for GpuAdvectionKernel<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> GpuAdvectionKernel<T> {
    /// Shader source code
    const SHADER_SOURCE: &'static str = include_str!("advection.wgsl");

    /// Creates a new GPU advection kernel
    #[must_use] pub fn new() -> Self {
        Self {
            shader_module: None,
            _phantom: PhantomData,
        }
    }

    /// Initialize the shader module with the given device
    pub fn initialize(&mut self, device: &wgpu::Device) {
        if self.shader_module.is_none() {
            let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Advection Shader"),
                source: wgpu::ShaderSource::Wgsl(Self::SHADER_SOURCE.into()),
            });
            self.shader_module = Some(module);
        }
    }

    /// Get the shader module, initializing if necessary
    ///
    /// # Panics
    /// Panics if shader compilation fails or if initialization has not been called
    pub fn get_shader_module(&mut self, device: &wgpu::Device) -> &wgpu::ShaderModule {
        self.initialize(device);
        self.shader_module
            .as_ref()
            .expect("Shader module should be initialized")
    }
}

impl<T: RealField + Copy> ComputeKernel<T> for GpuAdvectionKernel<T> {
    fn name(&self) -> &'static str {
        "GPU Advection (Upwind)"
    }

    fn execute(&self, _input: &[T], _output: &mut [T], _params: KernelParams) -> Result<()> {
        // GPU execution handled through GpuKernel trait
        Ok(())
    }

    fn complexity(&self, size: usize) -> usize {
        use crate::compute::gpu::constants::complexity;
        size * complexity::ADVECTION_UPWIND
    }

    fn supports_backend(&self, backend: &ComputeBackend) -> bool {
        matches!(backend, ComputeBackend::Gpu | ComputeBackend::Hybrid)
    }
}

impl<T: RealField + Copy> GpuKernel<T> for GpuAdvectionKernel<T> {
    fn shader_code(&self) -> &str {
        Self::SHADER_SOURCE
    }

    fn create_pipeline(&self, device: &wgpu::Device) -> Result<wgpu::ComputePipeline> {
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Advection Shader"),
            source: wgpu::ShaderSource::Wgsl(self.shader_code().into()),
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Advection Pipeline Layout"),
            bind_group_layouts: &[],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Advection Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader_module,
            entry_point: "advection_upwind",
        });

        Ok(pipeline)
    }

    fn dispatch(&self, encoder: &mut wgpu::CommandEncoder, params: KernelParams) {
        let (nx, ny, nz) = params.domain_params.grid_dims;
        let work_group_size = params.work_group_size;

        // Calculate dispatch dimensions using safe div_ceil
        let dispatch_x = nx.div_ceil(work_group_size);
        let dispatch_y = ny.div_ceil(work_group_size);
        let dispatch_z = nz.div_ceil(work_group_size);

        let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("Advection Pass"),
            timestamp_writes: None,
        });

        // Safe casting with bounds checking for GPU dispatch
        let dispatch_x = std::cmp::min(dispatch_x, u32::MAX as usize) as u32;
        let dispatch_y = std::cmp::min(dispatch_y, u32::MAX as usize) as u32;
        let dispatch_z = std::cmp::min(dispatch_z, u32::MAX as usize) as u32;

        compute_pass.dispatch_workgroups(dispatch_x, dispatch_y, dispatch_z);
    }
}

impl<T: RealField + Copy> std::fmt::Debug for GpuAdvectionKernel<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuAdvectionKernel").finish()
    }
}
