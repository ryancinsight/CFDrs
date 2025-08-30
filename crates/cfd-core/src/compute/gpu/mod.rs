//! GPU compute backend using wgpu

use crate::error::{Error, Result};
use std::sync::Arc;

pub mod buffer;
pub mod constants;
pub mod kernels;
pub mod pipeline;

pub use buffer::GpuBuffer;

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
                Error::InvalidConfiguration(format!("Failed to create GPU device: {}", e))
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

    /// Get adapter info
    pub fn adapter_info(&self) -> wgpu::AdapterInfo {
        self.adapter.get_info()
    }

    /// Check if GPU supports required features
    pub fn supports_features(&self, features: wgpu::Features) -> bool {
        self.adapter.features().contains(features)
    }

    /// Get maximum work group size
    pub fn max_work_group_size(&self) -> usize {
        self.limits.max_compute_invocations_per_workgroup as usize
    }

    /// Get maximum buffer size
    pub fn max_buffer_size(&self) -> usize {
        self.limits.max_buffer_size as usize
    }
}
