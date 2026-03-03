//! Viewport element — displays the 3D scene and handles mouse interaction.

use std::sync::Arc;
use gpui::{
    div, img, Context, ImageSource, InteractiveElement, IntoElement, MouseButton,
    MouseDownEvent, MouseMoveEvent, MouseUpEvent, ParentElement, ScrollWheelEvent, Styled,
};

use super::image_bridge::bgra_to_render_image;
use super::workspace::{gpui_rgba, Workspace};
use crate::domain::scene::graph::MeshHandle;
use crate::presentation::viewport::camera_controller;

impl Workspace {
    /// Render the viewport area with the cached image and mouse handlers.
    pub(crate) fn build_viewport_element(
        &self,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let theme = self.theme;
        let mut vp = div()
            .id("viewport")
            .flex_1()
            .min_w_0()
            .bg(gpui_rgba(theme.viewport_bg))
            .overflow_hidden()
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(Self::on_viewport_left_down),
            )
            .on_mouse_down(
                MouseButton::Middle,
                cx.listener(Self::on_viewport_middle_down),
            )
            .on_mouse_down(
                MouseButton::Right,
                cx.listener(Self::on_viewport_right_down),
            )
            .on_mouse_up(
                MouseButton::Left,
                cx.listener(Self::on_viewport_mouse_up),
            )
            .on_mouse_up(
                MouseButton::Middle,
                cx.listener(Self::on_viewport_mouse_up),
            )
            .on_mouse_up(
                MouseButton::Right,
                cx.listener(Self::on_viewport_mouse_up),
            )
            .on_mouse_move(cx.listener(Self::on_viewport_mouse_move))
            .on_scroll_wheel(cx.listener(Self::on_viewport_scroll));

        if let Some(image) = &self.viewport_image {
            vp = vp.child(
                img(ImageSource::Render(Arc::clone(image)))
                    .size_full(),
            );
        }

        vp
    }

    /// Re-render the scene via wgpu and cache the result as a `RenderImage`.
    pub(crate) fn refresh_viewport(&mut self) {
        let Some(renderer) = &mut self.renderer else {
            return;
        };

        // Collect visible mesh handles first (Copy types), then resolve to
        // references so we can hand them to the GPU uploader.
        let handles: Vec<MeshHandle> = self
            .state
            .document
            .scene
            .visible_meshes()
            .map(|(_idx, handle)| handle)
            .collect();

        let mesh_refs: Vec<&cfd_mesh::IndexedMesh<f64>> = handles
            .iter()
            .filter_map(|h| self.state.document.mesh(*h))
            .collect();

        renderer.upload_meshes(&mesh_refs);

        let camera = self.state.document.scene.camera();
        let bgra = renderer.render(camera);
        let w = renderer.width();
        let h = renderer.height();
        self.viewport_image = Some(bgra_to_render_image(bgra, w, h));
    }

    // -- Mouse event handlers -------------------------------------------------

    fn on_viewport_left_down(
        &mut self,
        event: &MouseDownEvent,
        _window: &mut gpui::Window,
        _cx: &mut Context<Self>,
    ) {
        let pos = event.position;
        self.state.camera_controller.mouse_down(
            camera_controller::MouseButton::Left,
            pos.x.to_f64(),
            pos.y.to_f64(),
        );
    }

    fn on_viewport_middle_down(
        &mut self,
        event: &MouseDownEvent,
        _window: &mut gpui::Window,
        _cx: &mut Context<Self>,
    ) {
        let pos = event.position;
        self.state.camera_controller.mouse_down(
            camera_controller::MouseButton::Middle,
            pos.x.to_f64(),
            pos.y.to_f64(),
        );
    }

    fn on_viewport_right_down(
        &mut self,
        event: &MouseDownEvent,
        _window: &mut gpui::Window,
        _cx: &mut Context<Self>,
    ) {
        let pos = event.position;
        self.state.camera_controller.mouse_down(
            camera_controller::MouseButton::Right,
            pos.x.to_f64(),
            pos.y.to_f64(),
        );
    }

    fn on_viewport_mouse_up(
        &mut self,
        _event: &MouseUpEvent,
        _window: &mut gpui::Window,
        _cx: &mut Context<Self>,
    ) {
        self.state.camera_controller.mouse_up();
    }

    fn on_viewport_mouse_move(
        &mut self,
        event: &MouseMoveEvent,
        _window: &mut gpui::Window,
        cx: &mut Context<Self>,
    ) {
        let pos = event.position;
        let camera = self.state.document.scene.camera_mut();
        let changed = self.state.camera_controller.mouse_move(
            pos.x.to_f64(),
            pos.y.to_f64(),
            camera,
        );
        if changed {
            self.refresh_viewport();
            cx.notify();
        }
    }

    fn on_viewport_scroll(
        &mut self,
        event: &ScrollWheelEvent,
        _window: &mut gpui::Window,
        cx: &mut Context<Self>,
    ) {
        let delta_pixels = event.delta.pixel_delta(gpui::px(1.0));
        let delta = delta_pixels.y.to_f64() * 0.01;
        let camera = self.state.document.scene.camera_mut();
        self.state.camera_controller.scroll(delta, camera);
        self.refresh_viewport();
        cx.notify();
    }
}
