//! Camera controller — interprets mouse events for orbit/pan/zoom and pick.
//!
//! Supports two camera backends: orbital (spherical coordinates) and trackball
//! (quaternion rotation). Left-click without drag (< 3px movement) emits a
//! pick event instead of orbiting.

use crate::domain::scene::camera::orbital::OrbitalCamera;
use crate::domain::scene::camera::trackball::TrackballCamera;

/// Mouse button state for camera interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MouseButton {
    Left,
    Middle,
    Right,
}

/// Camera interaction mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InteractionMode {
    /// Not currently interacting.
    Idle,
    /// Orbiting (left mouse drag).
    Orbiting,
    /// Panning (middle mouse drag).
    Panning,
}

/// Which camera backend is active.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CameraMode {
    #[default]
    Orbital,
    Trackball,
}

/// Translates mouse input into camera manipulations.
pub struct CameraController {
    mode: InteractionMode,
    camera_mode: CameraMode,
    last_x: f64,
    last_y: f64,
    down_x: f64,
    down_y: f64,
    orbit_sensitivity: f64,
    drag_threshold_px: f64,
    pending_pick: Option<(u32, u32)>,
}

impl CameraController {
    /// Create a new camera controller with default sensitivity.
    #[must_use]
    pub fn new() -> Self {
        Self {
            mode: InteractionMode::Idle,
            camera_mode: CameraMode::default(),
            last_x: 0.0,
            last_y: 0.0,
            down_x: 0.0,
            down_y: 0.0,
            orbit_sensitivity: 0.005,
            drag_threshold_px: 3.0,
            pending_pick: None,
        }
    }

    /// Set the active camera mode.
    pub fn set_camera_mode(&mut self, mode: CameraMode) {
        self.camera_mode = mode;
    }

    /// Current camera mode.
    #[must_use]
    pub fn camera_mode(&self) -> CameraMode {
        self.camera_mode
    }

    /// Handle mouse button press.
    pub fn mouse_down(&mut self, button: MouseButton, x: f64, y: f64) {
        self.last_x = x;
        self.last_y = y;
        self.down_x = x;
        self.down_y = y;
        self.pending_pick = None;
        self.mode = match button {
            MouseButton::Left => InteractionMode::Orbiting,
            MouseButton::Middle => InteractionMode::Panning,
            MouseButton::Right => InteractionMode::Panning,
        };
    }

    /// Handle mouse button release.
    pub fn mouse_up(&mut self) {
        // If the mouse barely moved since mouse_down, this is a pick (click).
        if self.mode == InteractionMode::Orbiting {
            let dx = self.last_x - self.down_x;
            let dy = self.last_y - self.down_y;
            let dist = (dx * dx + dy * dy).sqrt();
            if dist < self.drag_threshold_px {
                self.pending_pick = Some((self.down_x as u32, self.down_y as u32));
            }
        }
        self.mode = InteractionMode::Idle;
    }

    /// Handle mouse movement with an orbital camera. Returns true if modified.
    pub fn mouse_move_orbital(&mut self, x: f64, y: f64, camera: &mut OrbitalCamera) -> bool {
        let dx = x - self.last_x;
        let dy = y - self.last_y;
        self.last_x = x;
        self.last_y = y;

        match self.mode {
            InteractionMode::Idle => false,
            InteractionMode::Orbiting => {
                camera.orbit(-dx * self.orbit_sensitivity, -dy * self.orbit_sensitivity);
                true
            }
            InteractionMode::Panning => {
                camera.pan(dx, dy);
                true
            }
        }
    }

    /// Handle mouse movement with a trackball camera. Returns true if modified.
    pub fn mouse_move_trackball(&mut self, x: f64, y: f64, camera: &mut TrackballCamera) -> bool {
        let dx = x - self.last_x;
        let dy = y - self.last_y;
        self.last_x = x;
        self.last_y = y;

        match self.mode {
            InteractionMode::Idle => false,
            InteractionMode::Orbiting => {
                camera.rotate(dx * self.orbit_sensitivity, dy * self.orbit_sensitivity);
                true
            }
            InteractionMode::Panning => {
                camera.pan(dx, dy);
                true
            }
        }
    }

    /// Handle scroll wheel for orbital camera. Returns true if modified.
    pub fn scroll_orbital(&self, delta: f64, camera: &mut OrbitalCamera) -> bool {
        camera.zoom(delta);
        true
    }

    /// Handle scroll wheel for trackball camera. Returns true if modified.
    pub fn scroll_trackball(&self, delta: f64, camera: &mut TrackballCamera) -> bool {
        camera.zoom(delta);
        true
    }

    /// Consume the pending pick event, if any.
    pub fn take_pending_pick(&mut self) -> Option<(u32, u32)> {
        self.pending_pick.take()
    }

    /// Current interaction mode.
    #[must_use]
    pub fn mode(&self) -> InteractionMode {
        self.mode
    }
}

impl Default for CameraController {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn click_without_drag_generates_pick() {
        let mut ctrl = CameraController::new();
        ctrl.mouse_down(MouseButton::Left, 100.0, 200.0);
        ctrl.mouse_up();
        assert_eq!(ctrl.take_pending_pick(), Some((100, 200)));
    }

    #[test]
    fn drag_does_not_generate_pick() {
        let mut ctrl = CameraController::new();
        ctrl.mouse_down(MouseButton::Left, 100.0, 200.0);
        let mut cam = OrbitalCamera::default();
        ctrl.mouse_move_orbital(200.0, 300.0, &mut cam);
        ctrl.mouse_up();
        assert_eq!(ctrl.take_pending_pick(), None);
    }

    #[test]
    fn middle_click_does_not_orbit() {
        let mut ctrl = CameraController::new();
        ctrl.mouse_down(MouseButton::Middle, 100.0, 200.0);
        assert_eq!(ctrl.mode(), InteractionMode::Panning);
    }
}
