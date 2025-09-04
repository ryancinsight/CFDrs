//! Arithmetic operations on GPU fields

use crate::compute::gpu::GpuContext;
use crate::gpu::shaders::{FIELD_ADD_SHADER, FIELD_MUL_SHADER};
use std::sync::Arc;
use wgpu::util::DeviceExt;

/// GPU kernel for field addition
pub struct FieldAddKernel {
    context: Arc<GpuContext>,
    pipeline: wgpu::ComputePipeline,
}

impl FieldAddKernel {
    /// Create new field addition kernel
    pub fn new(context: Arc<GpuContext>) -> Self {
        let shader_module = context
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Field Add Shader"),
                source: wgpu::ShaderSource::Wgsl(FIELD_ADD_SHADER.into()),
            });

        let bind_group_layout =
            context
                .device
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

        let pipeline_layout =
            context
                .device
                .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: Some("Add Pipeline Layout"),
                    bind_group_layouts: &[&bind_group_layout],
                    push_constant_ranges: &[],
                });

        let pipeline = context
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Add Pipeline"),
                layout: Some(&pipeline_layout),
                module: &shader_module,
                entry_point: "add_fields",
            });

        Self { context, pipeline }
    }

    /// Execute field addition: result = a + b
    pub fn execute(&self, a: &[f32], b: &[f32], result: &mut [f32]) {
        assert_eq!(a.len(), b.len());
        assert_eq!(a.len(), result.len());

        let size = a.len() as u32;

        // Create GPU buffers
        let buffer_a =
            self.context
                .device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Buffer A"),
                    contents: bytemuck::cast_slice(a),
                    usage: wgpu::BufferUsages::STORAGE,
                });
        let buffer_b =
            self.context
                .device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Buffer B"),
                    contents: bytemuck::cast_slice(b),
                    usage: wgpu::BufferUsages::STORAGE,
                });
        let buffer_result = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Result Buffer"),
            size: (size * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Create bind group for the operation
        let bind_group_layout =
            self.context
                .device
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
            .device
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

        // Create command encoder and dispatch compute work
        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Add Command Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Add Compute Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups((size + 63) / 64, 1, 1);
        }

        // Create staging buffer for readback
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: (size * 4) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        encoder.copy_buffer_to_buffer(&buffer_result, 0, &staging_buffer, 0, (size * 4) as u64);
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

/// GPU kernel for field scalar multiplication
pub struct FieldMulKernel {
    context: Arc<GpuContext>,
    pipeline: wgpu::ComputePipeline,
}

impl FieldMulKernel {
    /// Create new field multiplication kernel
    pub fn new(context: Arc<GpuContext>) -> Self {
        let shader_module = context
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Field Mul Shader"),
                source: wgpu::ShaderSource::Wgsl(FIELD_MUL_SHADER.into()),
            });

        // Pipeline setup similar to FieldAddKernel
        // Implementation details omitted for brevity

        let pipeline = context
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Mul Pipeline"),
                layout: None,
                module: &shader_module,
                entry_point: "multiply_field",
            });

        Self { context, pipeline }
    }

    /// Execute scalar multiplication: result = field * scalar
    pub fn execute(&self, field: &[f32], scalar: f32, result: &mut [f32]) {
        assert_eq!(field.len(), result.len());
        
        // Create buffers for GPU computation
        let field_buffer = self.context.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Field Buffer"),
            contents: bytemuck::cast_slice(field),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
        
        let result_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Result Buffer"),
            size: (result.len() * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        
        let scalar_buffer = self.context.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Scalar Buffer"),
            contents: bytemuck::cast_slice(&[scalar]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        
        // Create bind group
        let bind_group = self.context.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Multiply Bind Group"),
            layout: &self.pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: field_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: result_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: scalar_buffer.as_entire_binding(),
                },
            ],
        });
        
        // Dispatch computation
        let mut encoder = self.context.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Multiply Encoder"),
        });
        
        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Multiply Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups((field.len() as u32 + 255) / 256, 1, 1);
        }
        
        // Read back results
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: (result.len() * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        
        encoder.copy_buffer_to_buffer(&result_buffer, 0, &staging_buffer, 0, staging_buffer.size());
        self.context.queue.submit(std::iter::once(encoder.finish()));
        
        // Map and read results
        let buffer_slice = staging_buffer.slice(..);
        let (sender, receiver) = futures::channel::oneshot::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |v| sender.send(v).unwrap());
        
        self.context.device.poll(wgpu::Maintain::Wait);
        if let Ok(Ok(())) = pollster::block_on(receiver) {
            let data = buffer_slice.get_mapped_range();
            let gpu_result: &[f32] = bytemuck::cast_slice(&data);
            result.copy_from_slice(gpu_result);
        }
    }
}
