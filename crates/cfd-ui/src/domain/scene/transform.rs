//! Transform helpers for scene operations.

use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

/// Create a translation-only isometry.
#[must_use]
pub fn translation(x: f64, y: f64, z: f64) -> Isometry3<f64> {
    Isometry3::from_parts(Translation3::new(x, y, z), UnitQuaternion::identity())
}

/// Create a rotation-only isometry around the Y axis.
#[must_use]
pub fn rotation_y(angle_rad: f64) -> Isometry3<f64> {
    Isometry3::from_parts(
        Translation3::identity(),
        UnitQuaternion::from_axis_angle(&Vector3::y_axis(), angle_rad),
    )
}

/// Create a uniform scale vector (isometry does not support scale directly).
#[must_use]
pub fn uniform_scale(s: f64) -> Vector3<f64> {
    Vector3::new(s, s, s)
}
