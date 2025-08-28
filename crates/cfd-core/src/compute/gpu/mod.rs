//! GPU compute backend using wgpu

use crate::error::{Error, Result};
use std::sync::Arc;
use wgpu::util::DeviceExt;

pub mod buffer;
pub mod kernels;

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
    pub fn new() -> Result<Self> {
        pollster::block_on(Self::new_async())
    }

    /// Async initialization
    async fn new_async() -> Result<Self> {
        let instance = wgpu::Instance::default();

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await
            .ok_or_else(|| {
                Error::InvalidConfiguration("No suitable GPU adapter found".to_string())
            })?;

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
