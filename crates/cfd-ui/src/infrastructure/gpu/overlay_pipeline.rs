//! Overlay render pipeline — unlit vertex-colored geometry rendered on top.
//!
//! Supports two topologies: `LineList` for axis indicators, measurement lines,
//! and wireframes; `TriangleList` for view cube faces and selection highlights.
//! Both share the same shader and uniform layout.

use bytemuck::Zeroable;
use wgpu::util::DeviceExt;

use super::overlay_buffer::OverlayBuffer;

/// Uniform data for the overlay pipeline (just view-projection).
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct OverlayUniforms {
    pub view_proj: [[f32; 4]; 4],
}

/// Which primitive topology the overlay pass should use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OverlayTopology {
    Lines,
    Triangles,
}

/// Overlay render pipeline with two sub-pipelines (lines and triangles).
pub struct OverlayPipeline {
    pub line_pipeline: wgpu::RenderPipeline,
    pub triangle_pipeline: wgpu::RenderPipeline,
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub uniform_buffer: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
}

impl OverlayPipeline {
    /// Create the overlay pipeline with line and triangle sub-pipelines.
    #[must_use]
    pub fn new(device: &wgpu::Device) -> Self {
        let shader_source = include_str!("../shaders/overlay.wgsl");
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("overlay shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("overlay uniforms layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("overlay uniforms"),
            contents: bytemuck::bytes_of(&OverlayUniforms::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("overlay bind group"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("overlay pipeline layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let line_pipeline =
            Self::create_sub_pipeline(device, &shader_module, &pipeline_layout, OverlayTopology::Lines);
        let triangle_pipeline =
            Self::create_sub_pipeline(device, &shader_module, &pipeline_layout, OverlayTopology::Triangles);

        Self {
            line_pipeline,
            triangle_pipeline,
            bind_group_layout,
            uniform_buffer,
            bind_group,
        }
    }

    /// Update the uniform buffer with a new view-projection matrix.
    pub fn update_uniforms(&self, queue: &wgpu::Queue, uniforms: &OverlayUniforms) {
        queue.write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(uniforms));
    }

    /// Draw overlay geometry into an existing render pass.
    pub fn draw<'a>(
        &'a self,
        pass: &mut wgpu::RenderPass<'a>,
        buffer: &'a OverlayBuffer,
        topology: OverlayTopology,
    ) {
        let pipeline = match topology {
            OverlayTopology::Lines => &self.line_pipeline,
            OverlayTopology::Triangles => &self.triangle_pipeline,
        };
        pass.set_pipeline(pipeline);
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.set_vertex_buffer(0, buffer.vertex_buffer.slice(..));
        pass.draw(0..buffer.vertex_count, 0..1);
    }

    fn create_sub_pipeline(
        device: &wgpu::Device,
        shader_module: &wgpu::ShaderModule,
        pipeline_layout: &wgpu::PipelineLayout,
        topology: OverlayTopology,
    ) -> wgpu::RenderPipeline {
        let prim_topology = match topology {
            OverlayTopology::Lines => wgpu::PrimitiveTopology::LineList,
            OverlayTopology::Triangles => wgpu::PrimitiveTopology::TriangleList,
        };

        device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some(match topology {
                OverlayTopology::Lines => "overlay line pipeline",
                OverlayTopology::Triangles => "overlay triangle pipeline",
            }),
            layout: Some(pipeline_layout),
            vertex: wgpu::VertexState {
                module: shader_module,
                entry_point: "vs_main",
                buffers: &[OverlayBuffer::vertex_buffer_layout()],
            },
            fragment: Some(wgpu::FragmentState {
                module: shader_module,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: wgpu::TextureFormat::Bgra8UnormSrgb,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: prim_topology,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None, // Overlays are double-sided
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: false,
                depth_compare: wgpu::CompareFunction::Always,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        })
    }
}
