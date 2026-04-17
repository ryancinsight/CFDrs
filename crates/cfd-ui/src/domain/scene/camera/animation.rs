//! Camera animation — smooth interpolation between camera states.
//!
//! # Theorem — SLERP Preserves Unit Norm
//!
//! Spherical linear interpolation between two unit quaternions `q0` and `q1`
//! produces a unit quaternion for all `t in [0, 1]`. The interpolated path
//! lies on the great arc of the 3-sphere.  ∎

use nalgebra::{Point3, UnitQuaternion, Vector3};

use super::OrbitalCamera;

/// Snapshot of camera state for interpolation.
#[derive(Clone, Debug)]
pub struct CameraSnapshot {
    /// Target point the camera orbits around.
    pub target: Point3<f64>,
    /// Distance from target to eye.
    pub distance: f64,
    /// Camera orientation as a unit quaternion.
    pub orientation: UnitQuaternion<f64>,
    /// Vertical field of view in radians.
    pub fov_y: f64,
}

impl CameraSnapshot {
    /// Capture the current state of an `OrbitalCamera`.
    #[must_use]
    pub fn from_orbital(cam: &OrbitalCamera) -> Self {
        let eye = cam.eye_position();
        let forward = (cam.target - eye).normalize();
        let orientation =
            UnitQuaternion::face_towards(&forward, &Vector3::y());
        Self {
            target: cam.target,
            distance: cam.distance,
            orientation,
            fov_y: cam.fov_y,
        }
    }

    /// Create a snapshot from azimuth/elevation angles (for named views).
    #[must_use]
    pub fn from_spherical(
        target: Point3<f64>,
        distance: f64,
        azimuth: f64,
        elevation: f64,
        fov_y: f64,
    ) -> Self {
        let cos_e = elevation.cos();
        let dir = Vector3::new(
            -(cos_e * azimuth.sin()),
            -elevation.sin(),
            -(cos_e * azimuth.cos()),
        );
        let orientation = UnitQuaternion::face_towards(&dir, &Vector3::y());
        Self {
            target,
            distance,
            orientation,
            fov_y,
        }
    }
}

/// Animates a camera between two snapshots over a fixed number of frames.
#[derive(Clone, Debug)]
pub struct CameraAnimator {
    start: Option<CameraSnapshot>,
    end: Option<CameraSnapshot>,
    current_frame: u32,
    total_frames: u32,
}

impl CameraAnimator {
    /// Create an idle animator.
    #[must_use]
    pub fn new() -> Self {
        Self {
            start: None,
            end: None,
            current_frame: 0,
            total_frames: 0,
        }
    }

    /// Begin animating from `start` to `end` over `duration_frames` frames.
    pub fn animate_to(&mut self, start: CameraSnapshot, end: CameraSnapshot, duration_frames: u32) {
        self.start = Some(start);
        self.end = Some(end);
        self.current_frame = 0;
        self.total_frames = duration_frames.max(1);
    }

    /// Whether the animator is currently running.
    #[must_use]
    pub fn is_active(&self) -> bool {
        self.start.is_some() && self.current_frame < self.total_frames
    }

    /// Advance one frame and apply the interpolated state to the camera.
    ///
    /// Returns `true` if the animation is still running after this tick.
    pub fn tick(&mut self, camera: &mut OrbitalCamera) -> bool {
        let (Some(start), Some(end)) = (&self.start, &self.end) else {
            return false;
        };

        self.current_frame += 1;
        let t = ease_in_out(f64::from(self.current_frame) / f64::from(self.total_frames));

        camera.target = lerp_point(&start.target, &end.target, t);
        camera.distance = start.distance + (end.distance - start.distance) * t;
        camera.fov_y = start.fov_y + (end.fov_y - start.fov_y) * t;

        // Recover azimuth/elevation from the interpolated orientation.
        let q = start.orientation.slerp(&end.orientation, t);
        let forward = q * (-Vector3::z());
        camera.azimuth = forward.x.atan2(forward.z);
        camera.elevation = (-forward.y).asin().clamp(
            -89.0_f64.to_radians(),
            89.0_f64.to_radians(),
        );

        let still_running = self.current_frame < self.total_frames;
        if !still_running {
            self.start = None;
            self.end = None;
        }
        still_running
    }
}

impl Default for CameraAnimator {
    fn default() -> Self {
        Self::new()
    }
}

/// Smooth ease-in/ease-out (cubic Hermite).
fn ease_in_out(t: f64) -> f64 {
    if t < 0.5 {
        4.0 * t * t * t
    } else {
        1.0 - (-2.0 * t + 2.0).powi(3) / 2.0
    }
}

/// Linear interpolation between two points.
fn lerp_point(a: &Point3<f64>, b: &Point3<f64>, t: f64) -> Point3<f64> {
    Point3::from(a.coords * (1.0 - t) + b.coords * t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn ease_in_out_endpoints() {
        assert_relative_eq!(ease_in_out(0.0), 0.0, epsilon = 1e-12);
        assert_relative_eq!(ease_in_out(1.0), 1.0, epsilon = 1e-12);
        assert_relative_eq!(ease_in_out(0.5), 0.5, epsilon = 1e-12);
    }

    #[test]
    fn animator_completes_in_given_frames() {
        let mut anim = CameraAnimator::new();
        let start = CameraSnapshot::from_orbital(&OrbitalCamera::default());
        let end = CameraSnapshot::from_spherical(
            Point3::origin(),
            5.0,
            0.0,
            0.0,
            std::f64::consts::FRAC_PI_4,
        );
        anim.animate_to(start, end, 10);
        let mut cam = OrbitalCamera::default();
        for _ in 0..9 {
            assert!(anim.tick(&mut cam));
        }
        assert!(!anim.tick(&mut cam));
        assert!(!anim.is_active());
    }
}
