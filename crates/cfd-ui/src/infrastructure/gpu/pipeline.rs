//! Mesh render pipeline — wgpu render pipeline for Phong-shaded mesh rendering
//! with up to 6 clip planes for section views.

use bytemuck::Zeroable;
use crate::infrastructure::gpu::clip_uniform::ClipUniforms;
use crate::infrastructure::gpu::mesh_buffer::GpuMeshBuffer;
use wgpu::util::DeviceExt;

/// Shading mode for the viewport renderer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u32)]
pub enum ShadingMode {
    /// Flat shading (per-face normals).
    Flat = 0,
    /// Smooth shading (interpolated vertex normals).
    Smooth = 1,
    /// Wireframe overlay.
    Wireframe = 2,
}

/// Uniform data sent to the vertex/fragment shaders each frame.
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct ViewUniforms {
    /// Combined view-projection matrix (column-major).
    pub view_proj: [[f32; 4]; 4],
    /// Camera eye position in world space.
    pub camera_pos: [f32; 3],
    pub _pad0: f32,
    /// Directional light direction (normalized, world space).
    pub light_dir: [f32; 3],
    pub _pad1: f32,
    /// Shading mode (0=flat, 1=smooth, 2=wireframe).
    pub shading_mode: u32,
    /// Field visualization min value.
    pub field_min: f32,
    /// Field visualization max value.
    pub field_max: f32,
    /// Whether field visualization is active (0=off, 1=on).
    pub field_active: u32,
}

/// The mesh render pipeline with two bind groups: view uniforms and clip planes.
pub struct MeshRenderPipeline {
    /// The wgpu render pipeline.
    pub pipeline: wgpu::RenderPipeline,
    /// Bind group layout for view uniforms (group 0).
    pub view_bind_group_layout: wgpu::BindGroupLayout,
    /// Bind group layout for clip uniforms (group 1).
    pub clip_bind_group_layout: wgpu::BindGroupLayout,
    /// View uniform buffer.
    pub uniform_buffer: wgpu::Buffer,
    /// Clip uniform buffer.
    pub clip_uniform_buffer: wgpu::Buffer,
    /// Bind group for view uniforms.
    pub bind_group: wgpu::BindGroup,
    /// Bind group for clip uniforms.
    pub clip_bind_group: wgpu::BindGroup,
}

impl MeshRenderPipeline {
    /// Create the mesh render pipeline with clip plane support.
    pub fn new(device: &wgpu::Device) -> Self {
        let shader_source = include_str!("../shaders/mesh.wgsl");
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("mesh shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        // Group 0: view uniforms.
        let view_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("mesh view uniforms layout"),
                entries: &[uniform_binding_entry(0)],
            });

        // Group 1: clip uniforms.
        let clip_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("mesh clip uniforms layout"),
                entries: &[uniform_binding_entry(0)],
            });

        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("mesh view uniforms"),
            contents: bytemuck::bytes_of(&ViewUniforms::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let clip_uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("mesh clip uniforms"),
            contents: bytemuck::bytes_of(&ClipUniforms::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("mesh view bind group"),
            layout: &view_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let clip_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("mesh clip bind group"),
            layout: &clip_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: clip_uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("mesh pipeline layout"),
            bind_group_layouts: &[&view_bind_group_layout, &clip_bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("mesh pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader_module,
                entry_point: "vs_main",
                buffers: &[GpuMeshBuffer::vertex_buffer_layout()],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader_module,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: wgpu::TextureFormat::Bgra8UnormSrgb,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        Self {
            pipeline,
            view_bind_group_layout,
            clip_bind_group_layout,
            uniform_buffer,
            clip_uniform_buffer,
            bind_group,
            clip_bind_group,
        }
    }

    /// Update the view uniform buffer.
    pub fn update_uniforms(&self, queue: &wgpu::Queue, uniforms: &ViewUniforms) {
        queue.write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(uniforms));
    }

    /// Update the clip plane uniform buffer.
    pub fn update_clip_uniforms(&self, queue: &wgpu::Queue, clip: &ClipUniforms) {
        queue.write_buffer(&self.clip_uniform_buffer, 0, bytemuck::bytes_of(clip));
    }
}

/// Helper to create a uniform buffer binding entry.
fn uniform_binding_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}
