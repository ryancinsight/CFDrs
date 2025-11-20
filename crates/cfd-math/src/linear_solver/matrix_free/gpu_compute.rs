//! GPU compute acceleration for matrix-free linear solvers.
//!
//! This module provides GPU-accelerated implementations using wgpu compute shaders
//! for CFD operators. It manages GPU context, buffer allocation, and shader execution.
//!
//! # Status: EXPERIMENTAL (MAJOR-011)
//!
//! **⚠️ Warning**: This module is currently EXPERIMENTAL and NOT integrated with
//! the main solver pipeline. GPU acceleration features are under active development.
//!
//! To enable GPU features (experimental):
//! ```bash
//! cargo build --features gpu
//! ```
//!
//! ## Integration Status
//! - ✅ GPU context management implemented
//! - ✅ Buffer allocation and transfer complete
//! - ✅ Compute shader infrastructure ready
//! - ❌ NOT YET integrated with iterative solvers
//! - ❌ Performance benchmarks pending
//!
//! ## Future Work
//! - Integration with GMRES/BiCGSTAB solvers
//! - Sparse matrix-vector product (SpMV) kernels
//! - Preconditioner application on GPU
//! - CPU/GPU performance crossover analysis

#[cfg(feature = "gpu")]
use std::sync::Arc;
#[cfg(feature = "gpu")]
use wgpu::{self, util::DeviceExt};

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
/// GPU configuration options
#[cfg(feature = "gpu")]
#[derive(Debug, Clone)]
pub struct Config {
    pub backends: wgpu::Backends,
    pub power_preference: wgpu::PowerPreference,
    pub enable_timestamps: bool,
}

/// Default GPU configuration: auto-detect backends, high performance adapter, timestamps enabled
#[cfg(feature = "gpu")]
impl Default for Config {
    fn default() -> Self {
        Self {
            backends: wgpu::Backends::PRIMARY,
            power_preference: wgpu::PowerPreference::HighPerformance,
            enable_timestamps: true,
        }
    }
}

#[cfg(feature = "gpu")]
impl GpuComputeContext {
    /// Create a new GPU compute context with default configuration
    pub async fn new() -> Result<Self, std::io::Error> {
        Self::new_with_config(Config::default()).await
    }

    /// Create a new GPU compute context with explicit configuration
    pub async fn new_with_config(config: Config) -> Result<Self, std::io::Error> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: config.backends,
            ..Default::default()
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: config.power_preference,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "No suitable GPU adapter found",
                )
            })?;

        let supported = adapter.features();
        let want_ts =
            config.enable_timestamps && supported.contains(wgpu::Features::TIMESTAMP_QUERY);
        let required_features = if want_ts {
            wgpu::Features::TIMESTAMP_QUERY
        } else {
            wgpu::Features::empty()
        };

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    required_features,
                    required_limits: wgpu::Limits::default(),
                    label: Some("CFD GPU Compute Device"),
                },
                None,
            )
            .await
            .map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("request_device failed: {:?}", e),
                )
            })?;

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
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label,
                contents: bytemuck::cast_slice(data),
                usage,
            })
    }

    /// Create an empty GPU buffer
    pub fn create_buffer(
        &self,
        size: u64,
        usage: wgpu::BufferUsages,
        label: Option<&str>,
    ) -> wgpu::Buffer {
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
    pub fn new(
        ctx: &GpuComputeContext,
        data: &[T],
        usage: wgpu::BufferUsages,
        label: Option<&str>,
    ) -> Self {
        let buffer = ctx.create_buffer_init(usage, data, label);
        Self {
            buffer,
            size: data.len(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create an empty GPU buffer
    pub fn new_empty(
        ctx: &GpuComputeContext,
        size: usize,
        usage: wgpu::BufferUsages,
        label: Option<&str>,
    ) -> Self {
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
        ctx.queue
            .write_buffer(&self.buffer, 0, bytemuck::cast_slice(data));
    }

    /// Copy data from GPU buffer to CPU (async)
    pub async fn read(&self, ctx: &GpuComputeContext) -> Vec<T> {
        let size = (self.size * std::mem::size_of::<T>()) as u64;
        let staging_buffer = ctx.create_buffer(
            size,
            wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            None,
        );

        let mut encoder = ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
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
        let shader = ctx
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Compute Shader"),
                source: wgpu::ShaderSource::Wgsl(shader_source.into()),
            });

        let bind_group_layout =
            ctx.device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
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

        let pipeline_layout = ctx
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("Compute Pipeline Layout"),
                bind_group_layouts: &[&bind_group_layout],
                push_constant_ranges: &[],
            });

        let pipeline = ctx
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
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
        let primary_input = input_buffers.get(0).expect("input buffer required");
        let secondary_input = input_buffers.get(1).unwrap_or(primary_input);

        let bind_group_entries = [
            wgpu::BindGroupEntry {
                binding: 0,
                resource: params_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: primary_input.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: secondary_input.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: output_buffer.as_entire_binding(),
            },
        ];

        let bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Compute Bind Group"),
            layout: &self.bind_group_layout,
            entries: &bind_group_entries,
        });

        let mut encoder = ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Compute Encoder"),
            });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Compute Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.pipeline);
            compute_pass.set_bind_group(0, &bind_group, &[]);
            compute_pass.dispatch_workgroups(workgroups.0, workgroups.1, workgroups.2);
        }

        ctx.queue.submit(std::iter::once(encoder.finish()));
    }

    /// Execute the compute shader and return simple wall-clock metrics
    pub fn execute_with_metrics(
        &self,
        ctx: &GpuComputeContext,
        params_buffer: &wgpu::Buffer,
        input_buffers: &[&wgpu::Buffer],
        output_buffer: &wgpu::Buffer,
        workgroups: (u32, u32, u32),
    ) -> DispatchMetrics {
        let start = std::time::Instant::now();
        self.execute(ctx, params_buffer, input_buffers, output_buffer, workgroups);
        ctx.device.poll(wgpu::Maintain::Wait);
        let elapsed = start.elapsed();
        DispatchMetrics {
            duration_ms: elapsed.as_secs_f64() * 1e3,
            timestamp_supported: ctx
                .adapter
                .features()
                .contains(wgpu::Features::TIMESTAMP_QUERY),
        }
    }
}

/// Metrics collected for a compute dispatch
#[cfg(feature = "gpu")]
#[derive(Debug, Clone, Copy)]
pub struct DispatchMetrics {
    pub duration_ms: f64,
    pub timestamp_supported: bool,
}

/// Stub implementations for when GPU feature is disabled
#[cfg(not(feature = "gpu"))]
#[derive(Debug)]
pub struct GpuComputeContext;

#[cfg(not(feature = "gpu"))]
impl GpuComputeContext {
    pub async fn new() -> Result<Self, std::io::Error> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "GPU feature not enabled",
        ))
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
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

#[cfg(not(feature = "gpu"))]
#[derive(Debug)]
pub struct ComputeShader;

