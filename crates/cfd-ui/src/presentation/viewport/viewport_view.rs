//! Viewport renderer — multi-pass 3D rendering with overlays and picking.
//!
//! Render order:
//! 1. Main mesh pass (Phong-shaded geometry)
//! 2. Pick pass (flat color-ID for selection, only on click)
//! 3. Selection highlight overlay (translucent face/edge highlights)
//! 4. Axis indicator overlay (bottom-left scissor rect)
//! 5. View cube overlay (top-right scissor rect)

use crate::domain::clipping::ClipPlaneSet;
use crate::domain::scene::camera::orbital::OrbitalCamera;
use crate::infrastructure::gpu::clip_uniform::ClipUniforms;
use crate::infrastructure::gpu::context::GpuViewportContext;
use crate::infrastructure::gpu::mesh_buffer::GpuMeshBuffer;
use crate::infrastructure::gpu::offscreen::OffscreenTarget;
use crate::infrastructure::gpu::overlay_buffer::OverlayBuffer;
use crate::infrastructure::gpu::overlay_pipeline::{
    OverlayPipeline, OverlayTopology, OverlayUniforms,
};
use crate::infrastructure::gpu::pick_pipeline::PickRenderPipeline;
use crate::infrastructure::gpu::pick_target::PickTarget;
use crate::infrastructure::gpu::pipeline::{MeshRenderPipeline, ShadingMode, ViewUniforms};
use cfd_mesh::IndexedMesh;

/// The 3D viewport renderer. Manages GPU resources and renders meshes offscreen.
pub struct ViewportRenderer {
    gpu: GpuViewportContext,
    pipeline: MeshRenderPipeline,
    overlay_pipeline: OverlayPipeline,
    pick_pipeline: PickRenderPipeline,
    target: OffscreenTarget,
    pick_target: PickTarget,
    mesh_buffers: Vec<GpuMeshBuffer>,
    shading_mode: ShadingMode,
}

impl ViewportRenderer {
    /// Create a new viewport renderer with the given initial dimensions.
    pub fn new(width: u32, height: u32) -> anyhow::Result<Self> {
        let gpu = GpuViewportContext::new()?;
        let pipeline = MeshRenderPipeline::new(&gpu.device);
        let overlay_pipeline = OverlayPipeline::new(&gpu.device);
        let pick_pipeline = PickRenderPipeline::new(&gpu.device);
        let target = OffscreenTarget::new(&gpu.device, width, height);
        let pick_target = PickTarget::new(&gpu.device, width, height);
        Ok(Self {
            gpu,
            pipeline,
            overlay_pipeline,
            pick_pipeline,
            target,
            pick_target,
            mesh_buffers: Vec::new(),
            shading_mode: ShadingMode::Smooth,
        })
    }

    /// Upload meshes to the GPU.
    pub fn upload_meshes(&mut self, meshes: &[&IndexedMesh<f64>]) {
        self.mesh_buffers.clear();
        for mesh in meshes {
            self.mesh_buffers
                .push(GpuMeshBuffer::from_indexed_mesh(mesh, &self.gpu.device));
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

    /// Resize the offscreen render targets.
    pub fn resize(&mut self, width: u32, height: u32) {
        self.target.resize(&self.gpu.device, width, height);
        self.pick_target.resize(&self.gpu.device, width, height);
    }

    /// Access the pick target for readback.
    #[must_use]
    pub fn pick_target(&self) -> &PickTarget {
        &self.pick_target
    }

    /// Access the GPU device.
    #[must_use]
    pub fn device(&self) -> &wgpu::Device {
        &self.gpu.device
    }

    /// Access the GPU queue.
    #[must_use]
    pub fn queue(&self) -> &wgpu::Queue {
        &self.gpu.queue
    }

    /// Render the full scene and return BGRA pixel data.
    ///
    /// `clip_planes` provides active section planes for mesh clipping.
    /// `overlay_buffers` contains optional overlay geometry for each pass.
    pub fn render(
        &self,
        camera: &OrbitalCamera,
        clip_planes: &ClipPlaneSet,
        overlays: &RenderOverlays,
    ) -> Vec<u8> {
        let aspect = f64::from(self.target.width) / f64::from(self.target.height.max(1));
        let view = camera.view_matrix();
        let proj = camera.projection_matrix(aspect);
        let view_proj = proj * view;

        let vp = mat4_to_f32(&view_proj);
        let eye = camera.eye_position();
        let light_dir = nalgebra::Vector3::new(0.3_f64, 1.0, 0.5).normalize();

        let uniforms = ViewUniforms {
            view_proj: vp,
            camera_pos: [eye.x as f32, eye.y as f32, eye.z as f32],
            pad0: 0.0,
            light_dir: [light_dir.x as f32, light_dir.y as f32, light_dir.z as f32],
            pad1: 0.0,
            shading_mode: self.shading_mode as u32,
            field_min: 0.0,
            field_max: 1.0,
            field_active: 0,
        };

        self.pipeline.update_uniforms(&self.gpu.queue, &uniforms);

        let clip_uniforms = ClipUniforms::from_clip_plane_set(clip_planes);
        self.pipeline
            .update_clip_uniforms(&self.gpu.queue, &clip_uniforms);

        let overlay_uniforms = OverlayUniforms { view_proj: vp };
        self.overlay_pipeline
            .update_uniforms(&self.gpu.queue, &overlay_uniforms);

        let mut encoder = self
            .gpu
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("viewport render"),
            });

        // Pass 1: Main mesh pass.
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("mesh pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.target.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.18,
                            g: 0.20,
                            b: 0.25,
                            a: 1.0,
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
            pass.set_bind_group(1, &self.pipeline.clip_bind_group, &[]);

            for buf in &self.mesh_buffers {
                pass.set_vertex_buffer(0, buf.vertex_buffer.slice(..));
                pass.set_index_buffer(buf.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..buf.index_count, 0, 0..1);
            }
        }

        // Pass 2: Selection highlight overlay.
        if let Some(sel_buf) = &overlays.selection_tris {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("selection overlay pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.target.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.target.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            self.overlay_pipeline
                .draw(&mut pass, sel_buf, OverlayTopology::Triangles);
        }

        // Pass 3: Axis indicator overlay.
        if let Some(axis_buf) = &overlays.axis_lines {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("axis indicator pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.target.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.target.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            self.overlay_pipeline
                .draw(&mut pass, axis_buf, OverlayTopology::Lines);
        }

        // Pass 4: View cube overlay.
        if let Some(cube_buf) = &overlays.view_cube_tris {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("view cube pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.target.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.target.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            self.overlay_pipeline
                .draw(&mut pass, cube_buf, OverlayTopology::Triangles);
        }

        self.gpu.queue.submit(std::iter::once(encoder.finish()));

        self.target.read_pixels(&self.gpu.device, &self.gpu.queue)
    }

    /// Execute a pick pass and read back the ID buffer.
    ///
    /// This is called separately from `render()` and only when a click is detected.
    pub fn render_pick_pass(&self, camera: &OrbitalCamera) {
        let aspect = f64::from(self.target.width) / f64::from(self.target.height.max(1));
        let view = camera.view_matrix();
        let proj = camera.projection_matrix(aspect);
        let view_proj = proj * view;

        let vp = mat4_to_f32(&view_proj);
        let eye = camera.eye_position();

        let uniforms = ViewUniforms {
            view_proj: vp,
            camera_pos: [eye.x as f32, eye.y as f32, eye.z as f32],
            pad0: 0.0,
            light_dir: [0.0, 1.0, 0.0],
            pad1: 0.0,
            shading_mode: 0,
            field_min: 0.0,
            field_max: 1.0,
            field_active: 0,
        };

        self.pick_pipeline
            .update_uniforms(&self.gpu.queue, &uniforms);

        let mut encoder = self
            .gpu
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("pick render"),
            });

        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("pick pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.pick_target.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.0,
                            g: 0.0,
                            b: 0.0,
                            a: 0.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.pick_target.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            pass.set_pipeline(&self.pick_pipeline.pipeline);
            pass.set_bind_group(0, &self.pick_pipeline.bind_group, &[]);

            for buf in &self.mesh_buffers {
                pass.set_vertex_buffer(0, buf.vertex_buffer.slice(..));
                pass.set_index_buffer(buf.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..buf.index_count, 0, 0..1);
            }
        }

        self.gpu.queue.submit(std::iter::once(encoder.finish()));
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

/// Optional overlay buffers to render in each frame.
#[derive(Default)]
pub struct RenderOverlays<'a> {
    /// Selection highlight triangles.
    pub selection_tris: Option<&'a OverlayBuffer>,
    /// Axis indicator lines.
    pub axis_lines: Option<&'a OverlayBuffer>,
    /// View cube face triangles.
    pub view_cube_tris: Option<&'a OverlayBuffer>,
}

/// Convert a nalgebra Matrix4<f64> to a [[f32;4];4] for GPU uniforms.
fn mat4_to_f32(m: &nalgebra::Matrix4<f64>) -> [[f32; 4]; 4] {
    let mut out = [[0.0f32; 4]; 4];
    for r in 0..4 {
        for c in 0..4 {
            out[r][c] = m[(r, c)] as f32;
        }
    }
    out
}
