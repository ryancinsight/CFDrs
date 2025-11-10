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
                    label: Some("Laplacian 2D Bind Group Layout"),
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
                    label: Some("Laplacian 2D Pipeline Layout"),
                    bind_group_layouts: &[&bind_group_layout],
                    push_constant_ranges: &[],
                });

        let pipeline = context
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Laplacian 2D Pipeline"),
                layout: Some(&pipeline_layout),
                module: &shader_module,
                entry_point: "laplacian_2d",
            });

        Self { context, pipeline }
    }

    /// Execute 2D Laplacian computation
    pub fn execute(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        result: &mut [f32],
    ) {
        // CPU fallback for small arrays or GPU initialization failures
        if field.len() < 64 || self.context.device.limits().max_storage_buffer_binding_size == 0 {
            self.execute_cpu(field, nx, ny, dx, dy, result);
            return;
        }

        match self.execute_gpu(field, nx, ny, dx, dy) {
            Ok(gpu_result) => result.copy_from_slice(&gpu_result),
            Err(_) => {
                // Fallback to CPU implementation
                self.execute_cpu(field, nx, ny, dx, dy, result);
            }
        }
    }

    fn execute_cpu(&self, field: &[f32], nx: usize, ny: usize, dx: f32, dy: f32, result: &mut [f32]) {
        let dx_inv2 = 1.0 / (dx * dx);
        let dy_inv2 = 1.0 / (dy * dy);
        
        for y in 1..ny-1 {
            for x in 1..nx-1 {
                let idx = y * nx + x;
                let laplacian = dx_inv2 * (field[(y-1) * nx + x] + field[(y+1) * nx + x] - 2.0 * field[idx]) +
                               dy_inv2 * (field[y * nx + (x-1)] + field[y * nx + (x+1)] - 2.0 * field[idx]);
                result[idx] = laplacian;
            }
        }
    }

    fn execute_gpu(&self, field: &[f32], nx: usize, ny: usize, dx: f32, dy: f32) -> Result<Vec<f32>, Box<dyn std::error::Error>> {
        let mut result = vec![0.0; field.len()];

        // Create uniform buffer for dimensions and grid spacing
        let uniforms = Laplacian2DUniforms {
            nx: nx as u32,
            ny: ny as u32,
            dx_inv2: 1.0 / (dx * dx),
            dy_inv2: 1.0 / (dy * dy),
        };
        let uniform_buffer = self.context.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Uniforms Buffer"),
            contents: bytemuck::cast_slice(&[uniforms]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        // Create field buffer
        let field_buffer = self.context.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Field Buffer"),
            contents: bytemuck::cast_slice(field),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });

        // Create result buffer
        let result_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Result Buffer"),
            size: (field.len() * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Create bind group using the pipeline's layout
        let bind_group = self.context.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Laplacian 2D Bind Group"),
            layout: &self.pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: field_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: uniform_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: result_buffer.as_entire_binding(),
                },
            ],
        });

        // Create staging buffer
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: (field.len() * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        // Execute compute pass
        let mut encoder = self
            .context
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Laplacian 2D Encoder"),
            });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Laplacian 2D Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups(((nx * ny) as f32 / 256.0).ceil() as u32, 1, 1);
        }

        encoder.copy_buffer_to_buffer(&result_buffer, 0, &staging_buffer, 0, (field.len() * std::mem::size_of::<f32>()) as u64);
        self.context.queue.submit(std::iter::once(encoder.finish()));

        // Read back results with timeout
        let (tx, rx) = std::sync::mpsc::channel();
        staging_buffer
            .slice(..)
            .map_async(wgpu::MapMode::Read, move |result| {
                tx.send(result).unwrap();
            });
        
        self.context.device.poll(wgpu::Maintain::Wait);
        
        // Wait for result with timeout
        match rx.recv_timeout(std::time::Duration::from_millis(100)) {
            Ok(Ok(())) => {
                let data = staging_buffer.slice(..).get_mapped_range();
                result.copy_from_slice(bytemuck::cast_slice(&data));
                staging_buffer.unmap();
                Ok(result)
            }
            _ => {
                staging_buffer.unmap();
                Err("GPU operation timed out or failed".into())
            }
        }
    }
}
