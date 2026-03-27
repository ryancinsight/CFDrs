//! TrackballCamera — quaternion-based camera with no gimbal lock.
//!
//! # Theorem — Quaternion Rotation Preserves Unit Norm
//!
//! Given unit quaternions `q_current` and `q_delta`, their product
//! `q_delta * q_current` remains unit because `|q*p| = |q|*|p| = 1`.
//! Periodic renormalization guards against floating-point drift.  ∎

use nalgebra::{Matrix4, Point3, UnitQuaternion, Vector3};

use super::CameraProjection;

/// A trackball camera using quaternion orientation (Shoemake arcball).
///
/// Unlike `OrbitalCamera`, this camera stores orientation as a unit
/// quaternion, eliminating gimbal lock and enabling arbitrary rotation.
#[derive(Clone, Debug)]
pub struct TrackballCamera {
    /// The point the camera orbits around.
    pub target: Point3<f64>,
    /// Distance from the target to the camera eye.
    pub distance: f64,
    /// Camera orientation as a unit quaternion.
    pub orientation: UnitQuaternion<f64>,
    /// Vertical field of view in radians.
    pub fov_y: f64,
    /// Near clipping plane distance.
    pub near: f64,
    /// Far clipping plane distance.
    pub far: f64,
}

impl Default for TrackballCamera {
    fn default() -> Self {
        // Default orientation: looking from +Z toward origin, Y up.
        let eye = Point3::new(0.0, 0.0, 5.0);
        let target = Point3::origin();
        let orientation =
            UnitQuaternion::face_towards(&(target - eye), &Vector3::y());
        Self {
            target,
            distance: 5.0,
            orientation,
            fov_y: std::f64::consts::FRAC_PI_4,
            near: 0.01,
            far: 1000.0,
        }
    }
}

impl TrackballCamera {
    /// Rotate the camera by a screen-space delta (Shoemake arcball).
    ///
    /// `dx` and `dy` are in radians-equivalent screen motion.
    pub fn rotate(&mut self, dx: f64, dy: f64) {
        let q_yaw = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), -dx);
        let local_right = self.orientation * Vector3::x();
        let right_axis = nalgebra::Unit::new_normalize(local_right);
        let q_pitch = UnitQuaternion::from_axis_angle(&right_axis, -dy);
        self.orientation = q_pitch * q_yaw * self.orientation;
        self.orientation.renormalize();
    }

    /// Pan the camera in the screen plane.
    pub fn pan(&mut self, dx: f64, dy: f64) {
        let right = self.orientation * Vector3::x();
        let up = self.orientation * Vector3::y();
        let scale = self.distance * 0.002;
        self.target += right * (-dx * scale) + up * (dy * scale);
    }

    /// Zoom by adjusting the distance.
    pub fn zoom(&mut self, delta: f64) {
        self.distance = (self.distance * (1.0 - delta * ZOOM_SENSITIVITY))
            .clamp(MIN_DISTANCE, MAX_DISTANCE);
    }
}

impl CameraProjection for TrackballCamera {
    fn view_matrix(&self) -> Matrix4<f64> {
        let eye = self.eye_position();
        let up = self.orientation * Vector3::y();
        Matrix4::look_at_rh(&eye, &self.target, &up)
    }

    fn projection_matrix(&self, aspect: f64) -> Matrix4<f64> {
        Matrix4::new_perspective(aspect, self.fov_y, self.near, self.far)
    }

    fn eye_position(&self) -> Point3<f64> {
        let forward = self.orientation * (-Vector3::z());
        self.target - forward * self.distance
    }

    fn fit_to_bounds(&mut self, min: &Point3<f64>, max: &Point3<f64>) {
        let center = nalgebra::center(min, max);
        let half_diag = (max - min).norm() * 0.5;
        self.target = center;
        self.distance = half_diag / (self.fov_y * 0.5).tan();
        self.near = self.distance * 0.01;
        self.far = self.distance * 100.0;
    }
}

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
    fn default_trackball_has_positive_distance() {
        let cam = TrackballCamera::default();
        assert!(cam.distance > 0.0);
    }

    #[test]
    fn view_matrix_is_orthonormal() {
        let cam = TrackballCamera::default();
        let v = cam.view_matrix();
        let r = Vector3::new(v[(0, 0)], v[(0, 1)], v[(0, 2)]);
        let u = Vector3::new(v[(1, 0)], v[(1, 1)], v[(1, 2)]);
        assert_relative_eq!(r.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(u.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(r.dot(&u), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn rotate_preserves_unit_quaternion() {
        let mut cam = TrackballCamera::default();
        cam.rotate(1.0, 0.5);
        cam.rotate(-0.3, 2.0);
        assert_relative_eq!(cam.orientation.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn zoom_clamps_distance() {
        let mut cam = TrackballCamera::default();
        cam.zoom(1e10);
        assert!(cam.distance >= MIN_DISTANCE);
        cam.zoom(-1e10);
        assert!(cam.distance <= MAX_DISTANCE);
    }
}
