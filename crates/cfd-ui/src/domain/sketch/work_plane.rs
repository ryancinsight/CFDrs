//! Work plane — 2D coordinate frame embedded in 3D space.
//!
//! # Theorem — Projection/Unprojection Inverses
//!
//! Given an orthonormal frame `(u, v, n)` with origin `O`, the projection
//! `P(w) = ((w - O) · u, (w - O) · v)` and unprojection
//! `U(s, t) = O + s·u + t·v` satisfy `P(U(s, t)) = (s, t)` for all
//! `(s, t) ∈ R²`, because `u · u = 1`, `u · v = 0`, etc.  ∎

use nalgebra::{Matrix4, Point2, Point3, Vector3};

/// A 2D coordinate frame embedded in 3D space for sketch drawing.
#[derive(Clone, Debug)]
pub struct WorkPlane {
    /// Origin point in world space.
    pub origin: Point3<f64>,
    /// Normal vector (out-of-plane direction).
    pub normal: Vector3<f64>,
    /// U-axis (horizontal in sketch space).
    pub u_axis: Vector3<f64>,
    /// V-axis (vertical in sketch space).
    pub v_axis: Vector3<f64>,
}

impl WorkPlane {
    /// Create a work plane from origin, normal, and U-axis.
    ///
    /// The V-axis is computed as `normal × u_axis` to ensure a right-handed frame.
    #[must_use]
    pub fn new(origin: Point3<f64>, normal: Vector3<f64>, u_axis: Vector3<f64>) -> Self {
        let n = normal.normalize();
        let u = (u_axis - n * u_axis.dot(&n)).normalize();
        let v = n.cross(&u);
        Self { origin, normal: n, u_axis: u, v_axis: v }
    }

    /// XY origin plane (normal = +Z, u = +X).
    #[must_use]
    pub fn xy() -> Self {
        Self {
            origin: Point3::origin(),
            normal: Vector3::z(),
            u_axis: Vector3::x(),
            v_axis: Vector3::y(),
        }
    }

    /// XZ origin plane (normal = +Y, u = +X).
    #[must_use]
    pub fn xz() -> Self {
        Self {
            origin: Point3::origin(),
            normal: Vector3::y(),
            u_axis: Vector3::x(),
            v_axis: -Vector3::z(),
        }
    }

    /// YZ origin plane (normal = +X, u = +Y).
    #[must_use]
    pub fn yz() -> Self {
        Self {
            origin: Point3::origin(),
            normal: Vector3::x(),
            u_axis: Vector3::y(),
            v_axis: Vector3::z(),
        }
    }

    /// Create a work plane from a triangular mesh face.
    ///
    /// The normal is the face normal, U-axis is along the first edge.
    #[must_use]
    pub fn from_face(v0: &Point3<f64>, v1: &Point3<f64>, v2: &Point3<f64>) -> Option<Self> {
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let n = e1.cross(&e2);
        let len = n.norm();
        if len < 1e-12 {
            return None; // Degenerate face.
        }
        Some(Self::new(*v0, n / len, e1.normalize()))
    }

    /// Project a 3D world point onto the sketch plane as 2D `(u, v)`.
    #[must_use]
    pub fn project(&self, world_point: &Point3<f64>) -> Point2<f64> {
        let d = world_point - self.origin;
        Point2::new(d.dot(&self.u_axis), d.dot(&self.v_axis))
    }

    /// Unproject a 2D sketch point back to 3D world coordinates.
    #[must_use]
    pub fn unproject(&self, sketch_point: &Point2<f64>) -> Point3<f64> {
        self.origin + self.u_axis * sketch_point.x + self.v_axis * sketch_point.y
    }

    /// The 4x4 model matrix that transforms sketch 2D coordinates to world 3D.
    ///
    /// Columns: [u_axis, v_axis, normal, origin] (homogeneous).
    #[must_use]
    pub fn model_matrix(&self) -> Matrix4<f64> {
        Matrix4::new(
            self.u_axis.x, self.v_axis.x, self.normal.x, self.origin.x,
            self.u_axis.y, self.v_axis.y, self.normal.y, self.origin.y,
            self.u_axis.z, self.v_axis.z, self.normal.z, self.origin.z,
            0.0,           0.0,           0.0,           1.0,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn xy_plane_project_unproject_roundtrip() {
        let wp = WorkPlane::xy();
        let world = Point3::new(3.0, 4.0, 0.0);
        let sketch = wp.project(&world);
        assert_relative_eq!(sketch.x, 3.0, epsilon = 1e-12);
        assert_relative_eq!(sketch.y, 4.0, epsilon = 1e-12);
        let back = wp.unproject(&sketch);
        assert_relative_eq!(back.x, world.x, epsilon = 1e-12);
        assert_relative_eq!(back.y, world.y, epsilon = 1e-12);
    }

    #[test]
    fn from_face_produces_valid_frame() {
        let v0 = Point3::new(0.0, 0.0, 0.0);
        let v1 = Point3::new(1.0, 0.0, 0.0);
        let v2 = Point3::new(0.0, 1.0, 0.0);
        let wp = WorkPlane::from_face(&v0, &v1, &v2).expect("valid face");
        assert_relative_eq!(wp.normal.dot(&wp.u_axis), 0.0, epsilon = 1e-12);
        assert_relative_eq!(wp.normal.dot(&wp.v_axis), 0.0, epsilon = 1e-12);
        assert_relative_eq!(wp.u_axis.dot(&wp.v_axis), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn model_matrix_fourth_column_is_origin() {
        let wp = WorkPlane::new(
            Point3::new(1.0, 2.0, 3.0),
            Vector3::z(),
            Vector3::x(),
        );
        let m = wp.model_matrix();
        assert_relative_eq!(m[(0, 3)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(m[(1, 3)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(m[(2, 3)], 3.0, epsilon = 1e-12);
    }
}
