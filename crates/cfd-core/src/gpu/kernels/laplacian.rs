//! Laplacian operator GPU implementation

use crate::compute::gpu::GpuContext;
use crate::gpu::shaders::LAPLACIAN_2D_SHADER;
use bytemuck::{Pod, Zeroable};
use std::sync::Arc;
use wgpu::util::DeviceExt;

/// Uniform parameters for 2D Laplacian
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct Laplacian2DUniforms {
    nx: u32,
    ny: u32,
    dx_inv2: f32,
    dy_inv2: f32,
}

/// GPU kernel for 2D Laplacian computation
pub struct Laplacian2DKernel {
    context: Arc<GpuContext>,
    pipeline: wgpu::ComputePipeline,
}

impl Laplacian2DKernel {
    /// Create new 2D Laplacian kernel
    #[must_use]
    pub fn new(context: Arc<GpuContext>) -> Self {
        let shader_module = context
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Laplacian 2D Shader"),
                source: wgpu::ShaderSource::Wgsl(LAPLACIAN_2D_SHADER.into()),
            });

        let bind_group_layout =
            context
                .device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Laplacian Bind Group Layout"),
                    entries: &[
                        wgpu::BindGroupLayoutEntry {
                            binding: 0,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        wgpu::BindGroupLayoutEntry {
                            binding: 1,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        wgpu::BindGroupLayoutEntry {
                            binding: 2,
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

        let pipeline_layout =
            context
                .device
                .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: Some("Laplacian Pipeline Layout"),
                    bind_group_layouts: &[&bind_group_layout],
                    push_constant_ranges: &[],
                });

        let pipeline = context
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Laplacian Pipeline"),
                layout: Some(&pipeline_layout),
                module: &shader_module,
                entry_point: "laplacian_2d",
            });

        Self { context, pipeline }
    }

    /// Compute 2D Laplacian of a field
    pub fn execute(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        result: &mut [f32],
    ) {
        assert_eq!(field.len(), nx * ny);
        assert_eq!(field.len(), result.len());

        let uniforms = Laplacian2DUniforms {
            nx: nx as u32,
            ny: ny as u32,
            dx_inv2: 1.0 / (dx * dx),
            dy_inv2: 1.0 / (dy * dy),
        };

        // Create GPU buffers
        let uniforms_buffer =
            self.context
                .device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Uniforms Buffer"),
                    contents: bytemuck::cast_slice(&[uniforms]),
                    usage: wgpu::BufferUsages::UNIFORM,
                });

        let field_buffer =
            self.context
                .device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Field Buffer"),
                    contents: bytemuck::cast_slice(field),
                    usage: wgpu::BufferUsages::STORAGE,
                });

        let result_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Result Buffer"),
            size: (field.len() * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Create bind group
        let bind_group_layout =
            self.context
                .device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Laplacian Bind Group Layout"),
                    entries: &[
                        wgpu::BindGroupLayoutEntry {
                            binding: 0,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        wgpu::BindGroupLayoutEntry {
                            binding: 1,
                            visibility: wgpu::ShaderStages::COMPUTE,
                            ty: wgpu::BindingType::Buffer {
                                ty: wgpu::BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        },
                        wgpu::BindGroupLayoutEntry {
                            binding: 2,
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

        let bind_group = self
            .context
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Laplacian Bind Group"),
                layout: &bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: field_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: uniforms_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: result_buffer.as_entire_binding(),
                    },
                ],
            });

        // Execute compute pass
        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Laplacian Command Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Laplacian Compute Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups((nx as u32).div_ceil(8), (ny as u32).div_ceil(8), 1);
        }

        // Copy results to staging buffer
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: (field.len() * 4) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        encoder.copy_buffer_to_buffer(
            &result_buffer,
            0,
            &staging_buffer,
            0,
            (field.len() * 4) as u64,
        );
        self.context.queue.submit(std::iter::once(encoder.finish()));

        // Read back results
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r); // Channel may be dropped if GPU context is destroyed
        });
        self.context.device.poll(wgpu::Maintain::Wait);

        // Handle potential GPU readback failure
        match rx.recv() {
            Ok(Ok(())) => {
                let data = buffer_slice.get_mapped_range();
                let float_data: &[f32] = bytemuck::cast_slice(&data);
                result.copy_from_slice(float_data);
                staging_buffer.unmap();
            }
            Ok(Err(_)) => {
                // GPU mapping failed - fill with NaN to indicate error
                result.fill(f32::NAN);
            }
            Err(_) => {
                // Channel communication failed - fill with NaN
                result.fill(f32::NAN);
            }
        }
    }
}
