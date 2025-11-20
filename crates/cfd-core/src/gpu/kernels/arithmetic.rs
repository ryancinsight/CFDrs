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
    #[must_use]
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

    /// Add two fields element-wise
    pub fn execute(&self, a: &[f32], b: &[f32], result: &mut [f32]) {
        assert_eq!(a.len(), b.len());
        assert_eq!(a.len(), result.len());

        // For small arrays, use CPU computation
        if a.len() < 64 {
            for i in 0..a.len() {
                result[i] = a[i] + b[i];
            }
            return;
        }

        // Try GPU computation with error handling
        match self.execute_gpu(a, b, result) {
            Ok(()) => {}
            Err(_) => {
                // Fall back to CPU computation on GPU failure
                for i in 0..a.len() {
                    result[i] = a[i] + b[i];
                }
            }
        }
    }

    /// Execute GPU computation with proper error handling
    fn execute_gpu(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<(), &'static str> {
        // Create GPU buffers with proper alignment
        let data_size = (a.len() * std::mem::size_of::<f32>()) as u64;
        let aligned_size = ((data_size + 255) / 256) * 256; // Align to 256 bytes

        let buffer_a = self
            .context
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Field A Buffer"),
                contents: bytemuck::cast_slice(a),
                usage: wgpu::BufferUsages::STORAGE,
            });

        let buffer_b = self
            .context
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Field B Buffer"),
                contents: bytemuck::cast_slice(b),
                usage: wgpu::BufferUsages::STORAGE,
            });

        let buffer_result = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Result Buffer"),
            size: aligned_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Create bind group using the pipeline's layout
        let bind_group = self
            .context
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Field Add Bind Group"),
                layout: &self.pipeline.get_bind_group_layout(0),
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

        // Execute compute pass with proper workgroup dispatch
        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Field Add Command Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Field Add Compute Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            // Ensure we don't dispatch more workgroups than needed
            let workgroups = ((a.len() + 63) / 64) as u32;
            compute_pass.dispatch_workgroups(workgroups, 1, 1);
        }

        // Copy results to staging buffer with proper alignment
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: aligned_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        encoder.copy_buffer_to_buffer(&buffer_result, 0, &staging_buffer, 0, data_size);
        self.context.queue.submit(std::iter::once(encoder.finish()));

        // Read back results with timeout
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });

        // Poll with timeout to prevent hanging
        let start = std::time::Instant::now();
        while rx.try_recv().is_err() {
            self.context.device.poll(wgpu::Maintain::Wait);
            if start.elapsed().as_secs() > 5 {
                return Err("GPU readback timeout".into());
            }
        }

        match rx.recv() {
            Ok(Ok(())) => {
                let data = buffer_slice.get_mapped_range();
                let float_data: &[f32] = bytemuck::cast_slice(&data);
                result.copy_from_slice(float_data);
                staging_buffer.unmap();
                Ok(())
            }
            Ok(Err(_)) | Err(_) => Err("GPU readback failed"),
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
    #[must_use]
    pub fn new(context: Arc<GpuContext>) -> Self {
        let shader_module = context
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Field Mul Shader"),
                source: wgpu::ShaderSource::Wgsl(FIELD_MUL_SHADER.into()),
            });

        let bind_group_layout =
            context
                .device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some("Field Mul Bind Group Layout"),
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
                    label: Some("Field Mul Pipeline Layout"),
                    bind_group_layouts: &[&bind_group_layout],
                    push_constant_ranges: &[],
                });

        let pipeline = context
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Field Mul Pipeline"),
                layout: Some(&pipeline_layout),
                module: &shader_module,
                entry_point: "multiply_field",
            });

        Self { context, pipeline }
    }

    /// Execute scalar multiplication: result = field * scalar
    pub fn execute(&self, field: &[f32], scalar: f32, result: &mut [f32]) {
        assert_eq!(field.len(), result.len());

        // CPU fallback for small arrays or GPU initialization failures
        if field.len() < 64 || self.context.device.limits().max_storage_buffer_binding_size == 0 {
            for i in 0..field.len() {
                result[i] = field[i] * scalar;
            }
            return;
        }

        match self.execute_gpu(field, scalar, result) {
            Ok(()) => {}
            Err(_) => {
                // Fallback to CPU implementation
                for i in 0..field.len() {
                    result[i] = field[i] * scalar;
                }
            }
        }
    }

    fn execute_gpu(
        &self,
        field: &[f32],
        scalar: f32,
        result: &mut [f32],
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create buffers with proper alignment
        let data_size = (field.len() * std::mem::size_of::<f32>()) as u64;
        let aligned_size = ((data_size + 255) / 256) * 256; // Align to 256 bytes

        let field_buffer =
            self.context
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Field Buffer"),
                    contents: bytemuck::cast_slice(field),
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                });

        let result_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Result Buffer"),
            size: aligned_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Create uniform buffer with proper alignment (16 bytes minimum for uniforms)
        let scalar_buffer =
            self.context
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Scalar Buffer"),
                    contents: bytemuck::cast_slice(&[scalar, 0.0, 0.0, 0.0]), // Pad to 16 bytes
                    usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
                });

        // Create bind group using the pipeline's layout
        let bind_group = self
            .context
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Multiply Bind Group"),
                layout: &self.pipeline.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: field_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: scalar_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: result_buffer.as_entire_binding(),
                    },
                ],
            });

        // Dispatch computation
        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Multiply Encoder"),
                });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Multiply Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            // Use consistent workgroup dispatch with FieldAddKernel
            let workgroups = ((field.len() + 63) / 64) as u32;
            compute_pass.dispatch_workgroups(workgroups, 1, 1);
        }

        // Read back results with proper alignment
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: aligned_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        encoder.copy_buffer_to_buffer(&result_buffer, 0, &staging_buffer, 0, data_size);
        self.context.queue.submit(std::iter::once(encoder.finish()));

        // Read back results with timeout
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });

        // Poll with timeout to prevent hanging
        let start = std::time::Instant::now();
        while rx.try_recv().is_err() {
            self.context.device.poll(wgpu::Maintain::Wait);
            if start.elapsed().as_secs() > 5 {
                return Err("GPU readback timeout".into());
            }
        }

        match rx.recv() {
            Ok(Ok(())) => {
                let data = buffer_slice.get_mapped_range();
                let gpu_result: &[f32] = bytemuck::cast_slice(&data);
                result.copy_from_slice(gpu_result);
                staging_buffer.unmap();
                Ok(())
            }
            Ok(Err(_)) | Err(_) => Err("GPU readback failed".into()),
        }
    }
}
