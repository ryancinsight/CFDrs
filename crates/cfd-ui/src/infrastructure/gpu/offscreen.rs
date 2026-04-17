//! Offscreen render target — texture + depth buffer + staging buffer for CPU readback.
//!
//! # Theorem — Pixel-Perfect Readback
//!
//! The staging buffer is sized to `padded_row_bytes * height` where
//! `padded_row_bytes = align_up(width * 4, COPY_BYTES_PER_ROW_ALIGNMENT)`.
//! After `copy_texture_to_buffer` and `map_async`, every pixel is guaranteed
//! to be available in the mapped slice because we call `device.poll(Maintain::Wait)`
//! which blocks until the GPU queue drains completely.  ∎

/// Offscreen render target for wgpu-based 3D rendering.
pub struct OffscreenTarget {
    /// Color attachment texture (`Bgra8UnormSrgb`).
    pub color_texture: wgpu::Texture,
    /// View into the color texture.
    pub color_view: wgpu::TextureView,
    /// Depth/stencil attachment texture (`Depth32Float`).
    pub depth_texture: wgpu::Texture,
    /// View into the depth texture.
    pub depth_view: wgpu::TextureView,
    /// Staging buffer for CPU readback of rendered pixels.
    staging_buffer: wgpu::Buffer,
    /// Render target width in pixels.
    pub width: u32,
    /// Render target height in pixels.
    pub height: u32,
    /// Row stride in the staging buffer (padded to alignment).
    padded_row_bytes: u32,
}

impl OffscreenTarget {
    /// Create a new offscreen target with the given dimensions.
    #[must_use]
    pub fn new(device: &wgpu::Device, width: u32, height: u32) -> Self {
        let (color_texture, color_view) = create_color_texture(device, width, height);
        let (depth_texture, depth_view) = create_depth_texture(device, width, height);
        let padded_row_bytes = padded_bytes_per_row(width);
        let staging_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("offscreen staging"),
            size: u64::from(padded_row_bytes * height),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        Self {
            color_texture,
            color_view,
            depth_texture,
            depth_view,
            staging_buffer,
            width,
            height,
            padded_row_bytes,
        }
    }

    /// Resize the offscreen target. Recreates all GPU resources.
    pub fn resize(&mut self, device: &wgpu::Device, width: u32, height: u32) {
        if width == self.width && height == self.height {
            return;
        }
        *self = Self::new(device, width, height);
    }

    /// Copy the color texture to the staging buffer and read back BGRA pixels.
    ///
    /// Returns a `Vec<u8>` of length `width * height * 4` in BGRA8 format.
    pub fn read_pixels(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> Vec<u8> {
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("offscreen readback"),
        });
        encoder.copy_texture_to_buffer(
            wgpu::ImageCopyTexture {
                texture: &self.color_texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::ImageCopyBuffer {
                buffer: &self.staging_buffer,
                layout: wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(self.padded_row_bytes),
                    rows_per_image: Some(self.height),
                },
            },
            wgpu::Extent3d {
                width: self.width,
                height: self.height,
                depth_or_array_layers: 1,
            },
        );
        queue.submit(std::iter::once(encoder.finish()));

        let buffer_slice = self.staging_buffer.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = sender.send(result);
        });
        device.poll(wgpu::Maintain::Wait);
        receiver.recv().expect("map_async channel closed").expect("map_async failed");

        let data = buffer_slice.get_mapped_range();
        let unpadded_row = self.width * 4;
        let mut pixels = Vec::with_capacity((unpadded_row * self.height) as usize);
        for row in 0..self.height {
            let start = (row * self.padded_row_bytes) as usize;
            let end = start + unpadded_row as usize;
            pixels.extend_from_slice(&data[start..end]);
        }
        drop(data);
        self.staging_buffer.unmap();
        pixels
    }
}

fn create_color_texture(
    device: &wgpu::Device,
    width: u32,
    height: u32,
) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("offscreen color"),
        size: wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Bgra8UnormSrgb,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    (texture, view)
}

fn create_depth_texture(
    device: &wgpu::Device,
    width: u32,
    height: u32,
) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("offscreen depth"),
        size: wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    (texture, view)
}

/// Compute padded bytes per row, aligned to `COPY_BYTES_PER_ROW_ALIGNMENT`.
fn padded_bytes_per_row(width: u32) -> u32 {
    let unpadded = width * 4;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    unpadded.div_ceil(align) * align
}
