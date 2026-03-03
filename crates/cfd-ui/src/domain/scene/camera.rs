//! OrbitalCamera — interactive camera with orbit, pan, and zoom control.
//!
//! # Theorem — View Matrix Orthonormality
//!
//! The view matrix produced by `view_matrix()` is always orthogonal (its
//! 3x3 upper-left block has orthonormal columns). **Proof sketch**: the
//! matrix is constructed from `look_at_rh` using unit-length right/up/forward
//! vectors derived from spherical coordinates, which form an orthonormal
//! basis in R^3 by construction.  ∎

use nalgebra::{Matrix4, Point3, Vector3};

/// An orbital camera that rotates around a target point.
///
/// Coordinates follow right-handed Y-up convention:
/// - Azimuth rotates around the Y axis (0 = +Z direction).
/// - Elevation rotates from the XZ plane toward +Y.
#[derive(Clone, Debug)]
pub struct OrbitalCamera {
    /// The point the camera orbits around.
    pub target: Point3<f64>,
    /// Distance from the target to the camera eye.
    pub distance: f64,
    /// Horizontal rotation angle in radians (0 = looking along +Z).
    pub azimuth: f64,
    /// Vertical rotation angle in radians (clamped to ±89°).
    pub elevation: f64,
    /// Vertical field of view in radians.
    pub fov_y: f64,
    /// Near clipping plane distance.
    pub near: f64,
    /// Far clipping plane distance.
    pub far: f64,
}

impl Default for OrbitalCamera {
    fn default() -> Self {
        Self {
            target: Point3::origin(),
            distance: 5.0,
            azimuth: std::f64::consts::FRAC_PI_4,
            elevation: std::f64::consts::FRAC_PI_6,
            fov_y: std::f64::consts::FRAC_PI_4,
            near: 0.01,
            far: 1000.0,
        }
    }
}

impl OrbitalCamera {
    /// Compute the eye position in world space from spherical coordinates.
    #[must_use]
    pub fn eye_position(&self) -> Point3<f64> {
        let cos_elev = self.elevation.cos();
        let x = self.distance * cos_elev * self.azimuth.sin();
        let y = self.distance * self.elevation.sin();
        let z = self.distance * cos_elev * self.azimuth.cos();
        self.target + Vector3::new(x, y, z)
    }

    /// Produce the 4x4 view matrix (world-to-camera transform).
    #[must_use]
    pub fn view_matrix(&self) -> Matrix4<f64> {
        let eye = self.eye_position();
        Matrix4::look_at_rh(&eye, &self.target, &Vector3::y())
    }

    /// Produce the 4x4 perspective projection matrix.
    #[must_use]
    pub fn projection_matrix(&self, aspect: f64) -> Matrix4<f64> {
        Matrix4::new_perspective(aspect, self.fov_y, self.near, self.far)
    }

    /// Rotate the camera around the target.
    pub fn orbit(&mut self, delta_azimuth: f64, delta_elevation: f64) {
        self.azimuth += delta_azimuth;
        self.elevation = (self.elevation + delta_elevation)
            .clamp(-MAX_ELEVATION, MAX_ELEVATION);
    }

    /// Pan the camera (translate target) in the screen plane.
    pub fn pan(&mut self, delta_x: f64, delta_y: f64) {
        let view = self.view_matrix();
        let right = Vector3::new(view[(0, 0)], view[(1, 0)], view[(2, 0)]);
        let up = Vector3::new(view[(0, 1)], view[(1, 1)], view[(2, 1)]);
        let scale = self.distance * 0.002;
        self.target += right * (-delta_x * scale) + up * (delta_y * scale);
    }

    /// Zoom by adjusting the camera distance.
    pub fn zoom(&mut self, delta: f64) {
        self.distance = (self.distance * (1.0 - delta * ZOOM_SENSITIVITY))
            .clamp(MIN_DISTANCE, MAX_DISTANCE);
    }

    /// Fit the camera so that an axis-aligned bounding box fills the view.
    pub fn fit_to_bounds(&mut self, min: &Point3<f64>, max: &Point3<f64>) {
        let center = nalgebra::center(min, max);
        let half_diag = (max - min).norm() * 0.5;
        self.target = center;
        self.distance = half_diag / (self.fov_y * 0.5).tan();
        self.near = self.distance * 0.01;
        self.far = self.distance * 100.0;
    }
}

/// Maximum elevation angle (89° in radians) to prevent gimbal lock.
const MAX_ELEVATION: f64 = 89.0_f64 * std::f64::consts::PI / 180.0;

/// Minimum eye distance to prevent clipping through the target.
const MIN_DISTANCE: f64 = 0.001;

/// Maximum eye distance.
const MAX_DISTANCE: f64 = 1e6;

/// Zoom wheel sensitivity factor.
const ZOOM_SENSITIVITY: f64 = 0.1;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn default_camera_has_positive_distance() {
        let cam = OrbitalCamera::default();
        assert!(cam.distance > 0.0);
    }

    #[test]
    fn view_matrix_upper_left_is_orthonormal() {
        let cam = OrbitalCamera::default();
        let v = cam.view_matrix();
        let r = Vector3::new(v[(0, 0)], v[(0, 1)], v[(0, 2)]);
        let u = Vector3::new(v[(1, 0)], v[(1, 1)], v[(1, 2)]);
        let f = Vector3::new(v[(2, 0)], v[(2, 1)], v[(2, 2)]);
        assert_relative_eq!(r.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(u.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(f.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(r.dot(&u), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn elevation_is_clamped() {
        let mut cam = OrbitalCamera::default();
        cam.orbit(0.0, 100.0);
        assert!(cam.elevation <= MAX_ELEVATION);
        cam.orbit(0.0, -200.0);
        assert!(cam.elevation >= -MAX_ELEVATION);
    }

    #[test]
    fn zoom_clamps_distance() {
        let mut cam = OrbitalCamera::default();
        cam.zoom(1e10);
        assert!(cam.distance >= MIN_DISTANCE);
        cam.zoom(-1e10);
        assert!(cam.distance <= MAX_DISTANCE);
    }
}
