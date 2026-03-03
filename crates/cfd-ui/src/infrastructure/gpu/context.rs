//! GPU viewport context — independent wgpu instance for 3D mesh rendering.
//!
//! This context is completely separate from any UI framework's GPU context.
//! It owns its own `wgpu::Instance`, `Adapter`, `Device`, and `Queue`, ensuring
//! that 3D rendering does not interfere with the UI compositor.

use anyhow::Context;

/// Independent GPU context for the 3D viewport renderer.
pub struct GpuViewportContext {
    /// The wgpu device for creating buffers, textures, and pipelines.
    pub device: wgpu::Device,
    /// The command queue for submitting GPU work.
    pub queue: wgpu::Queue,
}

impl GpuViewportContext {
    /// Create a new GPU context by requesting the best available adapter.
    ///
    /// Falls back through: discrete GPU -> integrated GPU -> software renderer.
    pub fn new() -> anyhow::Result<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .context("no suitable GPU adapter found")?;

        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("cfd-ui viewport"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
            },
            None,
        ))
        .context("failed to create GPU device")?;

        Ok(Self { device, queue })
    }
}
