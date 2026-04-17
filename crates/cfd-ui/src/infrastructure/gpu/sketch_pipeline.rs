//! Sketch render pipeline — 2D anti-aliased lines on a work plane.

use bytemuck::Zeroable;
use crate::infrastructure::gpu::overlay_buffer::OverlayBuffer;
use wgpu::util::DeviceExt;

/// Uniform data for sketch rendering.
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct SketchUniforms {
    /// Combined view-projection matrix.
    pub view_proj: [[f32; 4]; 4],
    /// Work plane model matrix (sketch 2D → world 3D).
    pub model: [[f32; 4]; 4],
    /// Viewport dimensions in pixels (for line width computation).
    pub viewport_size: [f32; 2],
    /// Line width in pixels.
    pub line_width_px: f32,
    pub pad: f32,
}

/// Sketch render pipeline with solid and dashed sub-pipelines.
pub struct SketchPipeline {
    pub solid_pipeline: wgpu::RenderPipeline,
    pub dashed_pipeline: wgpu::RenderPipeline,
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub uniform_buffer: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
}

impl SketchPipeline {
    /// Create the sketch pipeline.
    #[must_use]
    pub fn new(device: &wgpu::Device) -> Self {
        let shader_source = include_str!("../shaders/sketch.wgsl");
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("sketch shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("sketch uniforms layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("sketch uniforms"),
            contents: bytemuck::bytes_of(&SketchUniforms::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("sketch bind group"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("sketch pipeline layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let solid_pipeline =
            Self::create_sub_pipeline(device, &shader_module, &pipeline_layout, "fs_main");
        let dashed_pipeline =
            Self::create_sub_pipeline(device, &shader_module, &pipeline_layout, "fs_dashed");

        Self {
            solid_pipeline,
            dashed_pipeline,
            bind_group_layout,
            uniform_buffer,
            bind_group,
        }
    }

    /// Update the uniform buffer.
    pub fn update_uniforms(&self, queue: &wgpu::Queue, uniforms: &SketchUniforms) {
        queue.write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(uniforms));
    }

    /// Draw solid sketch geometry.
    pub fn draw_solid<'a>(
        &'a self,
        pass: &mut wgpu::RenderPass<'a>,
        buffer: &'a OverlayBuffer,
    ) {
        pass.set_pipeline(&self.solid_pipeline);
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.set_vertex_buffer(0, buffer.vertex_buffer.slice(..));
        pass.draw(0..buffer.vertex_count, 0..1);
    }

    /// Draw dashed construction geometry.
    pub fn draw_dashed<'a>(
        &'a self,
        pass: &mut wgpu::RenderPass<'a>,
        buffer: &'a OverlayBuffer,
    ) {
        pass.set_pipeline(&self.dashed_pipeline);
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.set_vertex_buffer(0, buffer.vertex_buffer.slice(..));
        pass.draw(0..buffer.vertex_count, 0..1);
    }

    fn create_sub_pipeline(
        device: &wgpu::Device,
        shader_module: &wgpu::ShaderModule,
        pipeline_layout: &wgpu::PipelineLayout,
        fragment_entry: &str,
    ) -> wgpu::RenderPipeline {
        device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("sketch pipeline"),
            layout: Some(pipeline_layout),
            vertex: wgpu::VertexState {
                module: shader_module,
                entry_point: "vs_main",
                buffers: &[OverlayBuffer::vertex_buffer_layout()],
            },
            fragment: Some(wgpu::FragmentState {
                module: shader_module,
                entry_point: fragment_entry,
                targets: &[Some(wgpu::ColorTargetState {
                    format: wgpu::TextureFormat::Bgra8UnormSrgb,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: false,
                depth_compare: wgpu::CompareFunction::LessEqual,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        })
    }
}
