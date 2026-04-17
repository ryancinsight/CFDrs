//! GPU mesh buffer — uploads `IndexedMesh` vertex/index data to the GPU.

use crate::domain::mesh_adapter::{self, MeshVertex};
use cfd_mesh::IndexedMesh;
use wgpu::util::DeviceExt;

/// GPU-resident vertex and index buffers for a single mesh.
pub struct GpuMeshBuffer {
    /// Vertex buffer containing interleaved `MeshVertex` data.
    pub vertex_buffer: wgpu::Buffer,
    /// Index buffer containing `u32` triangle indices.
    pub index_buffer: wgpu::Buffer,
    /// Number of indices (3 per triangle).
    pub index_count: u32,
}

impl GpuMeshBuffer {
    /// Upload an `IndexedMesh` to the GPU, converting to GPU-ready format.
    #[must_use]
    pub fn from_indexed_mesh(mesh: &IndexedMesh<f64>, device: &wgpu::Device) -> Self {
        let gpu_data = mesh_adapter::convert_mesh(mesh);

        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("mesh vertex buffer"),
            contents: bytemuck::cast_slice(&gpu_data.vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });

        let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("mesh index buffer"),
            contents: bytemuck::cast_slice(&gpu_data.indices),
            usage: wgpu::BufferUsages::INDEX,
        });

        Self {
            vertex_buffer,
            index_buffer,
            index_count: gpu_data.indices.len() as u32,
        }
    }

    /// The vertex buffer layout describing `MeshVertex` for the render pipeline.
    #[must_use]
    pub fn vertex_buffer_layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<MeshVertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                // position: vec3<f32>
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 0,
                    format: wgpu::VertexFormat::Float32x3,
                },
                // normal: vec3<f32>
                wgpu::VertexAttribute {
                    offset: 12,
                    shader_location: 1,
                    format: wgpu::VertexFormat::Float32x3,
                },
                // region_id: u32
                wgpu::VertexAttribute {
                    offset: 24,
                    shader_location: 2,
                    format: wgpu::VertexFormat::Uint32,
                },
                // field_value: f32
                wgpu::VertexAttribute {
                    offset: 28,
                    shader_location: 3,
                    format: wgpu::VertexFormat::Float32,
                },
            ],
        }
    }
}
