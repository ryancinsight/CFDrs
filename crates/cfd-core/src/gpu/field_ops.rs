//! GPU-accelerated field operations

use super::GpuContext;
use bytemuck::{Pod, Zeroable};
use std::sync::Arc;

/// GPU-accelerated 2D field operations
pub struct GpuFieldOps {
    context: Arc<GpuContext>,
    add_pipeline: wgpu::ComputePipeline,
    #[allow(dead_code)] // Will be used when multiply_field is implemented
    mul_pipeline: wgpu::ComputePipeline,
    #[allow(dead_code)] // Will be used when laplacian_2d is implemented
    laplacian_pipeline: wgpu::ComputePipeline,
}

impl GpuFieldOps {
    /// Create new GPU field operations handler
    pub fn new(context: Arc<GpuContext>) -> Self {
        // Create bind group layout for add operation (3 storage buffers)
        let add_bind_group_layout =
            context
                .device()
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Add Bind Group Layout"),
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
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
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

        // Create bind group layout for multiply operation (1 storage, 1 uniform, 1 storage)
        let mul_bind_group_layout =
            context
                .device()
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Multiply Bind Group Layout"),
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

        // Create bind group layout for laplacian (1 uniform, 1 storage, 1 storage)
        let laplacian_bind_group_layout =
            context
                .device()
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Laplacian Bind Group Layout"),
                    entries: &[
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

        let add_pipeline = context.create_compute_pipeline_with_layout(
            FIELD_ADD_SHADER,
            "add_fields",
            &add_bind_group_layout,
        );
        let mul_pipeline = context.create_compute_pipeline_with_layout(
            FIELD_MUL_SHADER,
            "multiply_field",
            &mul_bind_group_layout,
        );
        let laplacian_pipeline = context.create_compute_pipeline_with_layout(
            LAPLACIAN_SHADER,
            "compute_laplacian",
            &laplacian_bind_group_layout,
        );

        Self {
            context,
            add_pipeline,
            mul_pipeline,
            laplacian_pipeline,
        }
    }

    /// Add two fields on GPU
    pub fn add_fields(&self, a: &[f32], b: &[f32], result: &mut [f32]) {
        assert_eq!(a.len(), b.len());
        assert_eq!(a.len(), result.len());

        let size = a.len() as u32;

        // Create GPU buffers
        let buffer_a = self.context.create_buffer_init(a);
        let buffer_b = self.context.create_buffer_init(b);
        let buffer_result = self
            .context
            .device()
            .create_buffer(&wgpu::BufferDescriptor {
                label: Some("Result Buffer"),
                size: (size * 4) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });

        // Create bind group
        let bind_group_layout =
            self.context
                .device()
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Add Bind Group Layout"),
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
                                ty: wgpu::BufferBindingType::Storage { read_only: true },
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
            .device()
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Add Bind Group"),
                layout: &bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: buffer_a.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: buffer_b.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: buffer_result.as_entire_binding(),
                    },
                ],
            });

        // Dispatch compute
        let mut encoder =
            self.context
                .device()
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Add Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Add Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.add_pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups((size + 63) / 64, 1, 1);
        }

        // Copy result back
        let staging_buffer = self
            .context
            .device()
            .create_buffer(&wgpu::BufferDescriptor {
                label: Some("Staging Buffer"),
                size: (size * 4) as u64,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });

        encoder.copy_buffer_to_buffer(&buffer_result, 0, &staging_buffer, 0, (size * 4) as u64);

        self.context.queue().submit(Some(encoder.finish()));

        // Read back result
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |r| tx.send(r).unwrap());

        self.context.device().poll(wgpu::Maintain::Wait);
        rx.recv().unwrap().unwrap();

        {
            let data = buffer_slice.get_mapped_range();
            let float_data: &[f32] = bytemuck::cast_slice(&data);
            result.copy_from_slice(float_data);
        }
        staging_buffer.unmap();
    }

    /// Compute 2D Laplacian on GPU
    pub fn laplacian_2d(
        &self,
        field: &[f32],
        nx: u32,
        ny: u32,
        dx: f32,
        dy: f32,
        result: &mut [f32],
    ) {
        assert_eq!(field.len(), (nx * ny) as usize);
        assert_eq!(field.len(), result.len());

        // Create uniforms buffer
        #[repr(C)]
        #[derive(Copy, Clone, Pod, Zeroable)]
        struct Uniforms {
            nx: u32,
            ny: u32,
            dx_inv2: f32,
            dy_inv2: f32,
        }

        let uniforms = Uniforms {
            nx,
            ny,
            dx_inv2: 1.0 / (dx * dx),
            dy_inv2: 1.0 / (dy * dy),
        };

        let uniforms_buffer = self.context.create_buffer_init(&[uniforms]);
        let field_buffer = self.context.create_buffer_init(field);
        let result_buffer = self
            .context
            .device()
            .create_buffer(&wgpu::BufferDescriptor {
                label: Some("Laplacian Result"),
                size: (field.len() * 4) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });

        // Create bind group matching the laplacian pipeline layout
        let bind_group = self
            .context
            .device()
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Laplacian Bind Group"),
                layout: &self.laplacian_pipeline.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: uniforms_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: field_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: result_buffer.as_entire_binding(),
                    },
                ],
            });

        // Dispatch compute shader
        let mut encoder =
            self.context
                .device()
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Laplacian Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Laplacian Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.laplacian_pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups((nx + 7) / 8, (ny + 7) / 8, 1);
        }

        // Copy result to staging buffer and read back
        let staging_buffer = self
            .context
            .device()
            .create_buffer(&wgpu::BufferDescriptor {
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

        self.context
            .queue()
            .submit(std::iter::once(encoder.finish()));

        // Read back results
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |r| tx.send(r).unwrap());
        self.context.device().poll(wgpu::Maintain::Wait);
        rx.recv().unwrap().unwrap();

        {
            let data = buffer_slice.get_mapped_range();
            let float_data: &[f32] = bytemuck::cast_slice(&data);
            result.copy_from_slice(float_data);
        }
        staging_buffer.unmap();
    }
}

// WGSL shader for field addition
const FIELD_ADD_SHADER: &str = r#"
@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(64)
fn add_fields(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= arrayLength(&a)) {
        return;
    }
    result[index] = a[index] + b[index];
}
"#;

// WGSL shader for scalar multiplication
const FIELD_MUL_SHADER: &str = r#"
@group(0) @binding(0) var<storage, read> field: array<f32>;
@group(0) @binding(1) var<uniform> scalar: f32;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(64)
fn multiply_field(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= arrayLength(&field)) {
        return;
    }
    result[index] = field[index] * scalar;
}
"#;

// WGSL shader for 2D Laplacian
const LAPLACIAN_SHADER: &str = r#"
struct Uniforms {
    nx: u32,
    ny: u32,
    dx_inv2: f32,
    dy_inv2: f32,
}

@group(0) @binding(0) var<uniform> uniforms: Uniforms;
@group(0) @binding(1) var<storage, read> field: array<f32>;
@group(0) @binding(2) var<storage, read_write> result: array<f32>;

@compute @workgroup_size(8, 8)
fn compute_laplacian(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    
    if (i >= uniforms.nx || j >= uniforms.ny) {
        return;
    }
    
    let idx = j * uniforms.nx + i;
    
    // Handle boundaries with zero gradient
    var laplacian = 0.0;
    
    // X direction
    if (i > 0u && i < uniforms.nx - 1u) {
        let left = field[j * uniforms.nx + (i - 1u)];
        let center = field[idx];
        let right = field[j * uniforms.nx + (i + 1u)];
        laplacian += (left - 2.0 * center + right) * uniforms.dx_inv2;
    }
    
    // Y direction
    if (j > 0u && j < uniforms.ny - 1u) {
        let bottom = field[(j - 1u) * uniforms.nx + i];
        let center = field[idx];
        let top = field[(j + 1u) * uniforms.nx + i];
        laplacian += (bottom - 2.0 * center + top) * uniforms.dy_inv2;
    }
    
    result[idx] = laplacian;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_field_add() {
        if let Some(context) = GpuContext::new() {
            let ops = GpuFieldOps::new(Arc::new(context));

            let a = vec![1.0_f32; 1024];
            let b = vec![2.0_f32; 1024];
            let mut result = vec![0.0_f32; 1024];

            ops.add_fields(&a, &b, &mut result);

            for val in result {
                assert!((val - 3.0).abs() < 1e-6);
            }
        }
    }
}
