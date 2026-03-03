//! Camera controller — interprets mouse events for orbit/pan/zoom.

use crate::domain::scene::camera::OrbitalCamera;

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

/// Translates mouse input into camera manipulations.
pub struct CameraController {
    mode: InteractionMode,
    last_x: f64,
    last_y: f64,
    orbit_sensitivity: f64,
}

impl CameraController {
    /// Create a new camera controller with default sensitivity.
    #[must_use]
    pub fn new() -> Self {
        Self {
            mode: InteractionMode::Idle,
            last_x: 0.0,
            last_y: 0.0,
            orbit_sensitivity: 0.005,
        }
    }

    /// Handle mouse button press.
    pub fn mouse_down(&mut self, button: MouseButton, x: f64, y: f64) {
        self.last_x = x;
        self.last_y = y;
        self.mode = match button {
            MouseButton::Left => InteractionMode::Orbiting,
            MouseButton::Middle => InteractionMode::Panning,
            MouseButton::Right => InteractionMode::Panning,
        };
    }

    /// Handle mouse button release.
    pub fn mouse_up(&mut self) {
        self.mode = InteractionMode::Idle;
    }

    /// Handle mouse movement. Returns true if the camera was modified.
    pub fn mouse_move(&mut self, x: f64, y: f64, camera: &mut OrbitalCamera) -> bool {
        let dx = x - self.last_x;
        let dy = y - self.last_y;
        self.last_x = x;
        self.last_y = y;

        match self.mode {
            InteractionMode::Idle => false,
            InteractionMode::Orbiting => {
                camera.orbit(
                    -dx * self.orbit_sensitivity,
                    -dy * self.orbit_sensitivity,
                );
                true
            }
            InteractionMode::Panning => {
                camera.pan(dx, dy);
                true
            }
        }
    }

    /// Handle scroll wheel. Returns true if the camera was modified.
    pub fn scroll(&self, delta: f64, camera: &mut OrbitalCamera) -> bool {
        camera.zoom(delta);
        true
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
