//! Camera module — projection trait and concrete camera implementations.

use nalgebra::{Matrix4, Point3};

pub mod animation;
pub mod orbital;
pub mod trackball;

pub use orbital::OrbitalCamera;
pub use trackball::TrackballCamera;

/// Abstraction over camera types for view/projection computation (DIP).
pub trait CameraProjection {
    /// World-to-camera 4x4 view matrix.
    fn view_matrix(&self) -> Matrix4<f64>;
    /// Perspective or orthographic projection matrix.
    fn projection_matrix(&self, aspect: f64) -> Matrix4<f64>;
    /// Camera eye position in world space.
    fn eye_position(&self) -> Point3<f64>;
    /// Reposition the camera so that the given AABB fills the viewport.
    fn fit_to_bounds(&mut self, min: &Point3<f64>, max: &Point3<f64>);
}
