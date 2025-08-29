//! GPU compute pipeline management for CFD operations

use super::{GpuBuffer, GpuContext};
use crate::compute::traits::KernelParams;
use crate::error::{Error, Result};
use nalgebra::RealField;
use std::collections::HashMap;
use std::sync::Arc;
use wgpu::util::DeviceExt;

/// Manages GPU compute pipelines and resources
pub struct GpuPipelineManager {
    context: Arc<GpuContext>,
    pipelines: HashMap<String, wgpu::ComputePipeline>,
    bind_group_layouts: HashMap<String, wgpu::BindGroupLayout>,
}

impl GpuPipelineManager {
    /// Create a new pipeline manager
    pub fn new(context: Arc<GpuContext>) -> Self {
        Self {
            context,
            pipelines: HashMap::new(),
            bind_group_layouts: HashMap::new(),
        }
    }

    /// Register a compute pipeline
    pub fn register_pipeline(
        &mut self,
        name: &str,
        shader_source: &str,
        entry_point: &str,
    ) -> Result<()> {
        // Create shader module
        let shader_module =
            self.context
                .device
                .create_shader_module(wgpu::ShaderModuleDescriptor {
                    label: Some(&format!("{} Shader", name)),
                    source: wgpu::ShaderSource::Wgsl(shader_source.into()),
                });

        // Create bind group layout for standard CFD operations
        let bind_group_layout =
            self.context
                .device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some(&format!("{} Bind Group Layout", name)),
                    entries: &[
                        // Uniform buffer for grid parameters
                        wgpu::BindGroupLayoutEntry {
                            binding: 0,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // Input storage buffers
                        wgpu::BindGroupLayoutEntry {
                            binding: 1,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // Additional input buffers (velocity, divergence, etc.)
                        wgpu::BindGroupLayoutEntry {
                            binding: 2,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        // Output storage buffer
                        wgpu::BindGroupLayoutEntry {
                            binding: 3,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: false },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                    ],
                });

        // Create pipeline layout
        let pipeline_layout =
            self.context
                .device
                .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: Some(&format!("{} Pipeline Layout", name)),
                    bind_group_layouts: &[&bind_group_layout],
                    push_constant_ranges: &[],
                });

        // Create compute pipeline
        let pipeline =
            self.context
                .device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some(&format!("{} Pipeline", name)),
                    layout: Some(&pipeline_layout),
                    module: &shader_module,
                    entry_point,
                });

        self.bind_group_layouts
            .insert(name.to_string(), bind_group_layout);
        self.pipelines.insert(name.to_string(), pipeline);

        Ok(())
    }

    /// Execute a compute pipeline
    pub fn execute<T: RealField + Copy + bytemuck::Pod>(
        &self,
        pipeline_name: &str,
        input_buffers: &[&GpuBuffer<T>],
        output_buffer: &mut GpuBuffer<T>,
        params: &KernelParams,
    ) -> Result<()> {
        let pipeline = self
            .pipelines
            .get(pipeline_name)
            .ok_or_else(|| Error::InvalidInput(format!("Pipeline {} not found", pipeline_name)))?;

        let bind_group_layout = self.bind_group_layouts.get(pipeline_name).ok_or_else(|| {
            Error::InvalidInput(format!("Bind group layout {} not found", pipeline_name))
        })?;

        // Create uniform buffer for parameters
        let uniform_data = create_uniform_data(params);
        let uniform_buffer =
            self.context
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Uniform Buffer"),
                    contents: bytemuck::cast_slice(&uniform_data),
                    usage: wgpu::BufferUsages::UNIFORM,
                });

        // Create bind group
        let mut entries = vec![wgpu::BindGroupEntry {
            binding: 0,
            resource: uniform_buffer.as_entire_binding(),
        }];

        // Add input buffers
        for (i, buffer) in input_buffers.iter().enumerate() {
            entries.push(wgpu::BindGroupEntry {
                binding: (i + 1) as u32,
                resource: buffer.buffer.as_entire_binding(),
            });
        }

        // Add output buffer
        entries.push(wgpu::BindGroupEntry {
            binding: (input_buffers.len() + 1) as u32,
            resource: output_buffer.buffer.as_entire_binding(),
        });

        let bind_group = self
            .context
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Compute Bind Group"),
                layout: bind_group_layout,
                entries: &entries,
            });

        // Create command encoder and dispatch
        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Compute Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Compute Pass"),
                timestamp_writes: None,
            });

            compute_pass.set_pipeline(pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);

            // Calculate dispatch dimensions
            let (nx, ny, nz) = params.domain_params.grid_dims;
            let workgroup_size = params.work_group_size as u32;

            let dispatch_x = (nx as u32 + workgroup_size - 1) / workgroup_size;
            let dispatch_y = (ny as u32 + workgroup_size - 1) / workgroup_size;
            let dispatch_z = (nz as u32 + workgroup_size - 1) / workgroup_size;

            compute_pass.dispatch_workgroups(dispatch_x, dispatch_y, dispatch_z);
        }

        // Submit command buffer
        self.context.queue.submit(std::iter::once(encoder.finish()));

        Ok(())
    }

    /// Get GPU context
    pub fn context(&self) -> &Arc<GpuContext> {
        &self.context
    }
}

/// Create uniform data from kernel parameters
fn create_uniform_data(params: &KernelParams) -> Vec<f32> {
    let (nx, ny, nz) = params.domain_params.grid_dims;
    let (dx, dy, dz) = params.domain_params.grid_spacing;

    vec![
        nx as f32,
        ny as f32,
        nz as f32,
        0.0, // padding for alignment
        dx as f32,
        dy as f32,
        dz as f32,
        params.domain_params.dt as f32,
    ]
}
