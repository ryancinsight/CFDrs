//! Application entry — initializes the workspace and renders the window.
//!
//! This module is where gpui integration will connect. Currently it provides
//! the `CfdApp` struct that manages the `WorkspaceState` and viewport renderer.

use crate::presentation::workspace::WorkspaceState;
use crate::presentation::viewport::viewport_view::{RenderOverlays, ViewportRenderer};

/// The main application container.
pub struct CfdApp {
    /// Central workspace state.
    pub state: WorkspaceState,
    /// The 3D viewport renderer (optional, requires GPU).
    pub renderer: Option<ViewportRenderer>,
}

impl CfdApp {
    /// Create a new application instance.
    pub fn new() -> anyhow::Result<Self> {
        let state = WorkspaceState::new();
        let renderer = match ViewportRenderer::new(800, 600) {
            Ok(r) => Some(r),
            Err(e) => {
                tracing::warn!("failed to create viewport renderer: {e}");
                None
            }
        };
        Ok(Self { state, renderer })
    }

    /// Render the viewport to BGRA pixels.
    pub fn render_viewport(&self) -> Option<Vec<u8>> {
        let renderer = self.renderer.as_ref()?;
        let camera = self.state.document.scene.camera();
        Some(renderer.render(camera, &self.state.clip_planes, &RenderOverlays::default()))
    }
}
