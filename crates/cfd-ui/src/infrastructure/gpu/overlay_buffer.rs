//! Overlay GPU buffers — vertex-colored geometry for overlays.

use wgpu::util::DeviceExt;

/// A vertex for overlay rendering (position + RGBA color).
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct OverlayVertex {
    pub position: [f32; 3],
    pub color: [f32; 4],
}

/// A colored line segment for overlay rendering.
#[derive(Clone, Debug)]
pub struct OverlayLine {
    pub start: [f32; 3],
    pub end: [f32; 3],
    pub color: [f32; 4],
}

/// A colored triangle for overlay rendering (e.g. view cube faces).
#[derive(Clone, Debug)]
pub struct OverlayTriangle {
    pub v0: [f32; 3],
    pub v1: [f32; 3],
    pub v2: [f32; 3],
    pub color: [f32; 4],
}

/// GPU-resident overlay geometry.
pub struct OverlayBuffer {
    pub vertex_buffer: wgpu::Buffer,
    pub vertex_count: u32,
}

impl OverlayBuffer {
    /// Build a line-list overlay buffer from line segments.
    #[must_use]
    pub fn from_lines(lines: &[OverlayLine], device: &wgpu::Device) -> Self {
        let mut vertices = Vec::with_capacity(lines.len() * 2);
        for line in lines {
            vertices.push(OverlayVertex {
                position: line.start,
                color: line.color,
            });
            vertices.push(OverlayVertex {
                position: line.end,
                color: line.color,
            });
        }
        Self::from_vertices(&vertices, device)
    }

    /// Build a triangle-list overlay buffer from triangles.
    #[must_use]
    pub fn from_triangles(tris: &[OverlayTriangle], device: &wgpu::Device) -> Self {
        let mut vertices = Vec::with_capacity(tris.len() * 3);
        for tri in tris {
            vertices.push(OverlayVertex {
                position: tri.v0,
                color: tri.color,
            });
            vertices.push(OverlayVertex {
                position: tri.v1,
                color: tri.color,
            });
            vertices.push(OverlayVertex {
                position: tri.v2,
                color: tri.color,
            });
        }
        Self::from_vertices(&vertices, device)
    }

    /// Build from raw vertices.
    #[must_use]
    pub fn from_vertices(vertices: &[OverlayVertex], device: &wgpu::Device) -> Self {
        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("overlay vertex buffer"),
            contents: bytemuck::cast_slice(vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });
        Self {
            vertex_buffer,
            vertex_count: vertices.len() as u32,
        }
    }

    /// Vertex buffer layout for `OverlayVertex`.
    #[must_use]
    pub fn vertex_buffer_layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<OverlayVertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                // position: vec3<f32>
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 0,
                    format: wgpu::VertexFormat::Float32x3,
                },
                // color: vec4<f32>
                wgpu::VertexAttribute {
                    offset: 12,
                    shader_location: 1,
                    format: wgpu::VertexFormat::Float32x4,
                },
            ],
        }
    }
}
