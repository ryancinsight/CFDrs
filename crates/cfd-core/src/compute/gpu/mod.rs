//! GPU compute backend using wgpu

use crate::error::{Error, Result};
use std::sync::Arc;

pub mod buffer;
pub mod constants;
pub mod field_ops;
pub mod kernels;
pub mod pipeline;
pub mod poisson_solver;
pub mod shaders;
pub mod turbulence_compute;
pub mod validation_tests;

pub use buffer::GpuBuffer;
pub use field_ops::GpuFieldOps;
pub use poisson_solver::{GpuPoissonSolver, PoissonParams};
pub use turbulence_compute::{GpuTurbulenceCompute, TurbulencePerformanceInfo};

/// GPU context for managing device and queue
pub struct GpuContext {
    /// wgpu instance
    instance: wgpu::Instance,
    /// GPU adapter
    adapter: wgpu::Adapter,
    /// GPU device
    pub device: Arc<wgpu::Device>,
    /// Command queue
    pub queue: Arc<wgpu::Queue>,
    /// Device limits
    pub limits: wgpu::Limits,
}

impl GpuContext {
    /// Create a new GPU context
    ///
    /// # Errors
    /// Returns error if GPU device initialization fails or no suitable adapter found
    pub fn create() -> Result<Self> {
        pollster::block_on(Self::create_async())
    }

    /// Async initialization with support for integrated graphics
    async fn create_async() -> Result<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        // Try high performance first (discrete GPU if available)
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await
            // Fall back to low power (integrated graphics)
            .or_else(|| {
                pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::LowPower,
                    force_fallback_adapter: false,
                    compatible_surface: None,
                }))
            })
            // Last resort: software fallback
            .or_else(|| {
                pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::LowPower,
                    force_fallback_adapter: true,
                    compatible_surface: None,
                }))
            })
            .ok_or_else(|| {
                Error::InvalidConfiguration(
                    "No GPU adapter found (tried discrete, integrated, and software)".to_string(),
                )
            })?;

        // Log adapter info
        let info = adapter.get_info();
        tracing::info!(
            "GPU adapter selected: {} ({:?}) - Backend: {:?}",
            info.name,
            info.device_type,
            info.backend
        );

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("CFD GPU Device"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                },
                None,
            )
            .await
            .map_err(|e| {
                Error::InvalidConfiguration(format!("Failed to create GPU device: {e}"))
            })?;

        let limits = device.limits();

        Ok(Self {
            instance,
            adapter,
            device: Arc::new(device),
            queue: Arc::new(queue),
            limits,
        })
    }

    /// Get the wgpu instance for surface creation
    #[must_use]
    pub fn instance(&self) -> &wgpu::Instance {
        &self.instance
    }

    /// Get the GPU adapter information
    #[must_use]
    pub fn adapter(&self) -> &wgpu::Adapter {
        &self.adapter
    }

    /// Get adapter info
    #[must_use]
    pub fn adapter_info(&self) -> wgpu::AdapterInfo {
        self.adapter.get_info()
    }

    /// Check if GPU supports required features
    #[must_use]
    pub fn supports_features(&self, features: wgpu::Features) -> bool {
        self.adapter.features().contains(features)
    }

    /// Get maximum work group size
    #[must_use]
    pub fn max_work_group_size(&self) -> usize {
        self.limits.max_compute_invocations_per_workgroup as usize
    }

    /// Get maximum buffer size
    #[must_use]
    pub fn max_buffer_size(&self) -> usize {
        self.limits.max_buffer_size as usize
    }

    /// Create compute pipeline with explicit bind group layout
    #[must_use]
    pub fn create_compute_pipeline_with_layout(
        &self,
        shader_source: &str,
        entry_point: &str,
        bind_group_layout: &wgpu::BindGroupLayout,
    ) -> wgpu::ComputePipeline {
        let shader = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some(entry_point),
                source: wgpu::ShaderSource::Wgsl(shader_source.into()),
            });

        let pipeline_layout = self
            .device
            .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(&format!("{entry_point} Pipeline Layout")),
                bind_group_layouts: &[bind_group_layout],
                push_constant_ranges: &[],
            });

        self.device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(&format!("{entry_point} Pipeline")),
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point,
            })
    }
}
