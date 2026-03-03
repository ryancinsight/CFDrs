//! Workspace root view — owns application state and composes the CAD layout.

use std::sync::Arc;
use gpui::{
    div, rgba, Context, FocusHandle, Focusable, InteractiveElement, IntoElement,
    ParentElement, Render, Styled, Window,
};

use crate::presentation::theme::{ThemeColors, DARK_THEME};
use crate::presentation::workspace::WorkspaceState;
use crate::presentation::viewport::viewport_view::ViewportRenderer;

/// The root gpui view. Owns the workspace state, renderer, and cached viewport
/// image for display.
pub struct Workspace {
    pub(crate) state: WorkspaceState,
    pub(crate) renderer: Option<ViewportRenderer>,
    pub(crate) viewport_image: Option<Arc<gpui::RenderImage>>,
    pub(crate) theme: ThemeColors,
    pub(crate) focus_handle: FocusHandle,
}

impl Workspace {
    /// Create a new workspace, initializing the GPU renderer if available.
    pub fn new(cx: &mut Context<Self>) -> Self {
        let state = WorkspaceState::new();
        let renderer = match ViewportRenderer::new(800, 600) {
            Ok(r) => Some(r),
            Err(e) => {
                tracing::warn!("failed to create viewport renderer: {e}");
                None
            }
        };

        let mut ws = Self {
            state,
            renderer,
            viewport_image: None,
            theme: DARK_THEME,
            focus_handle: cx.focus_handle(),
        };
        ws.refresh_viewport();
        ws
    }
}

/// Convert a `[u8; 4]` RGBA color to a gpui `Rgba` via a packed u32.
pub(crate) fn gpui_rgba(color: [u8; 4]) -> gpui::Rgba {
    rgba(
        ((color[0] as u32) << 24)
            | ((color[1] as u32) << 16)
            | ((color[2] as u32) << 8)
            | (color[3] as u32),
    )
}

impl Focusable for Workspace {
    fn focus_handle(&self, _cx: &gpui::App) -> FocusHandle {
        self.focus_handle.clone()
    }
}

impl Render for Workspace {
    fn render(&mut self, _window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let theme = self.theme;

        // Build the toolbar.
        let toolbar = super::toolbar::render_toolbar(&theme);

        // Build the viewport area with mouse handlers.
        let viewport_area = self.build_viewport_element(cx);

        // Build side panels.
        let model_tree = super::panels::render_model_tree(&self.state, &theme);
        let properties = super::panels::render_properties(&self.state, &theme);

        // Build the console.
        let console = super::panels::render_console(&self.state.console, &theme);

        // Compose the CAD layout with action handlers on the root element.
        div()
            .id("workspace-root")
            .track_focus(&self.focus_handle)
            .flex()
            .flex_col()
            .size_full()
            .bg(gpui_rgba(theme.viewport_bg))
            .on_action(cx.listener(Self::handle_undo))
            .on_action(cx.listener(Self::handle_redo))
            .on_action(cx.listener(Self::handle_fit_view))
            .on_action(cx.listener(Self::handle_toggle_wireframe))
            .on_action(cx.listener(Self::handle_new_project))
            .on_action(cx.listener(Self::handle_open_project))
            .on_action(cx.listener(Self::handle_save_project))
            .on_action(cx.listener(Self::handle_import_stl))
            .on_action(cx.listener(Self::handle_import_mesh))
            .on_action(cx.listener(Self::handle_export_stl))
            .on_action(cx.listener(Self::handle_export_obj))
            .on_action(cx.listener(Self::handle_export_ply))
            .on_action(cx.listener(Self::handle_export_glb))
            .on_action(cx.listener(Self::handle_export_dxf))
            .on_action(cx.listener(Self::handle_export_drawing))
            .on_action(cx.listener(Self::handle_export_openfoam))
            .on_action(cx.listener(Self::handle_create_cube))
            .on_action(cx.listener(Self::handle_create_cylinder))
            .on_action(cx.listener(Self::handle_create_sphere))
            .on_action(cx.listener(Self::handle_create_cone))
            .on_action(cx.listener(Self::handle_create_torus))
            .on_action(cx.listener(Self::handle_create_pipe))
            .on_action(cx.listener(Self::handle_create_ellipsoid))
            .on_action(cx.listener(Self::handle_create_capsule))
            .on_action(cx.listener(Self::handle_create_frustum))
            .on_action(cx.listener(Self::handle_create_geodesic_sphere))
            .on_action(cx.listener(Self::handle_create_rounded_cube))
            .on_action(cx.listener(Self::handle_create_elbow))
            .on_action(cx.listener(Self::handle_csg_union))
            .on_action(cx.listener(Self::handle_csg_intersection))
            .on_action(cx.listener(Self::handle_csg_difference))
            .on_action(cx.listener(Self::handle_delete))
            .on_action(cx.listener(Self::handle_select_all))
            .on_action(cx.listener(Self::handle_run_simulation))
            .on_action(cx.listener(Self::handle_stop_simulation))
            // Toolbar row.
            .child(toolbar)
            // Middle row: model tree | viewport | properties.
            .child(
                div()
                    .flex()
                    .flex_row()
                    .flex_1()
                    .min_h_0()
                    .child(model_tree)
                    .child(viewport_area)
                    .child(properties),
            )
            // Bottom console.
            .child(console)
    }
}
