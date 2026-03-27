//! Pick target — offscreen render target for GPU color-ID picking.
//!
//! A separate render target that stores encoded entity IDs as pixel colors.
//! After rendering the pick pass, a single pixel readback identifies the
//! entity under the cursor.

/// Offscreen target for pick-pass rendering with CPU pixel readback.
pub struct PickTarget {
    pub color_texture: wgpu::Texture,
    pub color_view: wgpu::TextureView,
    pub depth_texture: wgpu::Texture,
    pub depth_view: wgpu::TextureView,
    staging_buffer: wgpu::Buffer,
    pub width: u32,
    pub height: u32,
    padded_row_bytes: u32,
}

impl PickTarget {
    /// Create a new pick target with the given dimensions.
    pub fn new(device: &wgpu::Device, width: u32, height: u32) -> Self {
        let (color_texture, color_view) = create_pick_color_texture(device, width, height);
        let (depth_texture, depth_view) = create_pick_depth_texture(device, width, height);
        let padded_row_bytes = padded_bytes_per_row(width);
        let staging_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pick staging"),
            size: (padded_row_bytes * height) as u64,
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

    /// Resize the pick target.
    pub fn resize(&mut self, device: &wgpu::Device, width: u32, height: u32) {
        if width == self.width && height == self.height {
            return;
        }
        *self = Self::new(device, width, height);
    }

    /// Read the encoded ID at a single pixel. Returns `(node_idx, face_idx)`.
    pub fn read_pixel(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        x: u32,
        y: u32,
    ) -> Option<(u32, u32)> {
        if x >= self.width || y >= self.height {
            return None;
        }

        let pixels = self.read_all_pixels(device, queue);
        let offset = ((y * self.width + x) * 4) as usize;
        if offset + 3 >= pixels.len() {
            return None;
        }

        decode_pick_pixel([pixels[offset], pixels[offset + 1], pixels[offset + 2], pixels[offset + 3]])
    }

    /// Read the entire pick buffer back to CPU.
    fn read_all_pixels(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> Vec<u8> {
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("pick readback"),
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
        receiver
            .recv()
            .expect("pick map_async channel closed")
            .expect("pick map_async failed");

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

/// Decode a BGRA pick pixel into `(node_idx, face_idx)`.
///
/// Encoding (from pick.wgsl):
/// - R = node & 0xFF, G = (node >> 8) & 0xFF
/// - B = face & 0xFF, A = (face >> 8) & 0xFF
///
/// But the texture is Bgra8UnormSrgb, so BGRA byte order:
/// byte[0]=B, byte[1]=G, byte[2]=R, byte[3]=A
fn decode_pick_pixel(bgra: [u8; 4]) -> Option<(u32, u32)> {
    let r = bgra[2] as u32;
    let g = bgra[1] as u32;
    let b = bgra[0] as u32;
    let a = bgra[3] as u32;

    let node_idx = r | (g << 8);
    let face_idx = b | (a << 8);

    // Background clear color (0,0,0,0 or all-zero) means no hit.
    if node_idx == 0 && face_idx == 0 {
        return None;
    }

    Some((node_idx, face_idx))
}

fn create_pick_color_texture(
    device: &wgpu::Device,
    width: u32,
    height: u32,
) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("pick color"),
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

fn create_pick_depth_texture(
    device: &wgpu::Device,
    width: u32,
    height: u32,
) -> (wgpu::Texture, wgpu::TextureView) {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("pick depth"),
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

fn padded_bytes_per_row(width: u32) -> u32 {
    let unpadded = width * 4;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    unpadded.div_ceil(align) * align
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_zero_is_none() {
        assert_eq!(decode_pick_pixel([0, 0, 0, 0]), None);
    }

    #[test]
    fn decode_nonzero_round_trips() {
        // Encode node=5, face=10 in BGRA order
        let r = 5u8; // node & 0xFF
        let g = 0u8; // (node >> 8) & 0xFF
        let b = 10u8; // face & 0xFF
        let a = 0u8; // (face >> 8) & 0xFF
        // BGRA: [B, G, R, A]
        let bgra = [b, g, r, a];
        assert_eq!(decode_pick_pixel(bgra), Some((5, 10)));
    }
}
