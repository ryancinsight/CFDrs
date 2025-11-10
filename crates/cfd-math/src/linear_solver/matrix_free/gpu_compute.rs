//! GPU compute acceleration for matrix-free linear solvers.
//!
//! This module provides GPU-accelerated implementations using wgpu compute shaders
//! for CFD operators. It manages GPU context, buffer allocation, and shader execution.

#[cfg(feature = "gpu")]
use wgpu::{self, util::DeviceExt};
#[cfg(feature = "gpu")]
use std::sync::Arc;

/// GPU compute context for managing wgpu resources
#[cfg(feature = "gpu")]
#[derive(Debug)]
pub struct GpuComputeContext {
    /// WGPU instance
    instance: wgpu::Instance,
    /// GPU adapter
    adapter: wgpu::Adapter,
    /// GPU device
    device: wgpu::Device,
    /// Command queue
    queue: wgpu::Queue,
}

#[cfg(feature = "gpu")]
impl GpuComputeContext {
    /// Create a new GPU compute context
    pub async fn new() -> Result<Self, wgpu::Error> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::PRIMARY,
            ..Default::default()
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .ok_or(wgpu::Error::RequestAdapterFailed)?;

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    label: Some("CFD GPU Compute Device"),
                },
                None,
            )
            .await?;

        Ok(Self {
            instance,
            adapter,
            device,
            queue,
        })
    }

    /// Get device reference
    pub fn device(&self) -> &wgpu::Device {
        &self.device
    }

    /// Get queue reference
    pub fn queue(&self) -> &wgpu::Queue {
        &self.queue
    }

    /// Create a GPU buffer from data
    pub fn create_buffer_init<T: bytemuck::Pod>(
        &self,
        usage: wgpu::BufferUsages,
        data: &[T],
        label: Option<&str>,
    ) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label,
            contents: bytemuck::cast_slice(data),
            usage,
        })
    }

    /// Create an empty GPU buffer
    pub fn create_buffer(&self, size: u64, usage: wgpu::BufferUsages, label: Option<&str>) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label,
            size,
            usage,
            mapped_at_creation: false,
        })
    }
}

/// GPU buffer for CFD data
#[cfg(feature = "gpu")]
#[derive(Debug)]
pub struct GpuBuffer<T: bytemuck::Pod> {
    buffer: wgpu::Buffer,
    size: usize,
    _phantom: std::marker::PhantomData<T>,
}

#[cfg(feature = "gpu")]
impl<T: bytemuck::Pod> GpuBuffer<T> {
    /// Create a new GPU buffer from data
    pub fn new(ctx: &GpuComputeContext, data: &[T], usage: wgpu::BufferUsages, label: Option<&str>) -> Self {
        let buffer = ctx.create_buffer_init(usage, data, label);
        Self {
            buffer,
            size: data.len(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create an empty GPU buffer
    pub fn new_empty(ctx: &GpuComputeContext, size: usize, usage: wgpu::BufferUsages, label: Option<&str>) -> Self {
        let buffer = ctx.create_buffer((size * std::mem::size_of::<T>()) as u64, usage, label);
        Self {
            buffer,
            size,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Copy data from CPU to GPU buffer
    pub fn write(&self, ctx: &GpuComputeContext, data: &[T]) {
        assert_eq!(data.len(), self.size, "Data size mismatch");
        ctx.queue.write_buffer(&self.buffer, 0, bytemuck::cast_slice(data));
    }

    /// Copy data from GPU buffer to CPU (async)
    pub async fn read(&self, ctx: &GpuComputeContext) -> Vec<T> {
        let size = (self.size * std::mem::size_of::<T>()) as u64;
        let staging_buffer = ctx.create_buffer(size, wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, None);

        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Read Buffer Encoder"),
        });

        encoder.copy_buffer_to_buffer(&self.buffer, 0, &staging_buffer, 0, size);
        ctx.queue.submit(std::iter::once(encoder.finish()));

        let buffer_slice = staging_buffer.slice(..);
        let (sender, receiver) = futures::channel::oneshot::channel();

        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).unwrap();
        });

        ctx.device.poll(wgpu::Maintain::Wait);
        receiver.await.unwrap().unwrap();

        let data = buffer_slice.get_mapped_range();
        let result: Vec<T> = bytemuck::cast_slice(&data).to_vec();
        drop(data);
        staging_buffer.unmap();

        result
    }

    /// Get buffer reference
    pub fn buffer(&self) -> &wgpu::Buffer {
        &self.buffer
    }

    /// Get buffer size in elements
    pub fn size(&self) -> usize {
        self.size
    }
}

/// Compute shader for CFD operations
#[cfg(feature = "gpu")]
#[derive(Debug)]
pub struct ComputeShader {
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

#[cfg(feature = "gpu")]
impl ComputeShader {
    /// Create a compute shader from WGSL source
    pub fn new(ctx: &GpuComputeContext, shader_source: &str, entry_point: &str) -> Self {
        let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Compute Shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        let bind_group_layout = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Compute Bind Group Layout"),
            entries: &[
                // Uniform buffer (params)
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
                // Input buffer
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
                // Additional input buffer (optional)
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
                // Output buffer
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

        let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Compute Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Compute Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point,
        });

        Self {
            pipeline,
            bind_group_layout,
        }
    }

    /// Execute the compute shader
    pub fn execute(
        &self,
        ctx: &GpuComputeContext,
        params_buffer: &wgpu::Buffer,
        input_buffers: &[&wgpu::Buffer],
        output_buffer: &wgpu::Buffer,
        workgroups: (u32, u32, u32),
    ) {
        let mut bind_group_entries = vec![
            wgpu::BindGroupEntry {
                binding: 0,
                resource: params_buffer.as_entire_binding(),
            },
        ];

        // Add input buffers
        for (i, buffer) in input_buffers.iter().enumerate() {
            bind_group_entries.push(wgpu::BindGroupEntry {
                binding: 1 + i as u32,
                resource: buffer.as_entire_binding(),
            });
        }

        // Add output buffer (binding 3)
        if input_buffers.len() < 3 {
            for _ in input_buffers.len()..2 {
                bind_group_entries.push(wgpu::BindGroupEntry {
                    binding: 1 + bind_group_entries.len() as u32 - 1,
                    resource: input_buffers[0].as_entire_binding(), // Dummy
                });
            }
        }
        bind_group_entries.push(wgpu::BindGroupEntry {
            binding: 3,
            resource: output_buffer.as_entire_binding(),
        });

        let bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Compute Bind Group"),
            layout: &self.bind_group_layout,
            entries: &bind_group_entries,
        });

        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Compute Encoder"),
        });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Compute Pass"),
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups(workgroups.0, workgroups.1, workgroups.2);
        }

        ctx.queue.submit(std::iter::once(encoder.finish()));
    }
}

/// Stub implementations for when GPU feature is disabled
#[cfg(not(feature = "gpu"))]
#[derive(Debug)]
pub struct GpuComputeContext;

#[cfg(not(feature = "gpu"))]
impl GpuComputeContext {
    pub async fn new() -> Result<Self, std::io::Error> {
        Err(std::io::Error::new(std::io::ErrorKind::Unsupported, "GPU feature not enabled"))
    }
}

#[cfg(not(feature = "gpu"))]
#[derive(Debug)]
pub struct GpuBuffer<T> {
    _phantom: std::marker::PhantomData<T>,
}

#[cfg(not(feature = "gpu"))]
impl<T> GpuBuffer<T> {
    pub fn new(_ctx: &GpuComputeContext, _data: &[T], _usage: (), _label: Option<&str>) -> Self {
        Self { _phantom: std::marker::PhantomData }
    }
}

#[cfg(not(feature = "gpu"))]
#[derive(Debug)]
pub struct ComputeShader;





