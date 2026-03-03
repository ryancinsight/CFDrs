//! Viewport renderer — renders the 3D scene to a pixel buffer for display.

use crate::domain::scene::camera::OrbitalCamera;
use crate::infrastructure::gpu::context::GpuViewportContext;
use crate::infrastructure::gpu::mesh_buffer::GpuMeshBuffer;
use crate::infrastructure::gpu::offscreen::OffscreenTarget;
use crate::infrastructure::gpu::pipeline::{MeshRenderPipeline, ShadingMode, ViewUniforms};
use cfd_mesh::IndexedMesh;

/// The 3D viewport renderer. Manages GPU resources and renders meshes offscreen.
pub struct ViewportRenderer {
    gpu: GpuViewportContext,
    pipeline: MeshRenderPipeline,
    target: OffscreenTarget,
    mesh_buffers: Vec<GpuMeshBuffer>,
    shading_mode: ShadingMode,
}

impl ViewportRenderer {
    /// Create a new viewport renderer with the given initial dimensions.
    pub fn new(width: u32, height: u32) -> anyhow::Result<Self> {
        let gpu = GpuViewportContext::new()?;
        let pipeline = MeshRenderPipeline::new(&gpu.device);
        let target = OffscreenTarget::new(&gpu.device, width, height);
        Ok(Self {
            gpu,
            pipeline,
            target,
            mesh_buffers: Vec::new(),
            shading_mode: ShadingMode::Smooth,
        })
    }

    /// Upload meshes to the GPU.
    pub fn upload_meshes(&mut self, meshes: &[&IndexedMesh<f64>]) {
        self.mesh_buffers.clear();
        for mesh in meshes {
            self.mesh_buffers.push(GpuMeshBuffer::from_indexed_mesh(mesh, &self.gpu.device));
        }
    }

    /// Set the shading mode.
    pub fn set_shading_mode(&mut self, mode: ShadingMode) {
        self.shading_mode = mode;
    }

    /// Current shading mode.
    #[must_use]
    pub fn shading_mode(&self) -> ShadingMode {
        self.shading_mode
    }

    /// Resize the offscreen render target.
    pub fn resize(&mut self, width: u32, height: u32) {
        self.target.resize(&self.gpu.device, width, height);
    }

    /// Render the scene and return BGRA pixel data.
    pub fn render(&self, camera: &OrbitalCamera) -> Vec<u8> {
        let aspect = self.target.width as f64 / self.target.height.max(1) as f64;
        let view = camera.view_matrix();
        let proj = camera.projection_matrix(aspect);
        let view_proj = proj * view;

        // Convert to f32 4x4 array.
        let mut vp: [[f32; 4]; 4] = [[0.0; 4]; 4];
        for r in 0..4 {
            for c in 0..4 {
                vp[r][c] = view_proj[(r, c)] as f32;
            }
        }

        let eye = camera.eye_position();
        let light_dir = nalgebra::Vector3::new(0.3_f64, 1.0, 0.5).normalize();

        let uniforms = ViewUniforms {
            view_proj: vp,
            camera_pos: [eye.x as f32, eye.y as f32, eye.z as f32],
            _pad0: 0.0,
            light_dir: [light_dir.x as f32, light_dir.y as f32, light_dir.z as f32],
            _pad1: 0.0,
            shading_mode: self.shading_mode as u32,
            field_min: 0.0,
            field_max: 1.0,
            field_active: 0,
        };

        self.pipeline.update_uniforms(&self.gpu.queue, &uniforms);

        // Encode render pass.
        let mut encoder = self.gpu.device.create_command_encoder(
            &wgpu::CommandEncoderDescriptor { label: Some("viewport render") },
        );

        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("mesh pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.target.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.18, g: 0.20, b: 0.25, a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.target.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            pass.set_pipeline(&self.pipeline.pipeline);
            pass.set_bind_group(0, &self.pipeline.bind_group, &[]);

            for buf in &self.mesh_buffers {
                pass.set_vertex_buffer(0, buf.vertex_buffer.slice(..));
                pass.set_index_buffer(buf.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..buf.index_count, 0, 0..1);
            }
        }

        self.gpu.queue.submit(std::iter::once(encoder.finish()));

        self.target.read_pixels(&self.gpu.device, &self.gpu.queue)
    }

    /// Width of the render target.
    #[must_use]
    pub fn width(&self) -> u32 {
        self.target.width
    }

    /// Height of the render target.
    #[must_use]
    pub fn height(&self) -> u32 {
        self.target.height
    }
}
