//! Plane representation and polygon splitting.
//!
//! Adapted from csgrs's BSP plane, but operating on indexed vertices rather
//! than owned `Polygon<S>` structs. The plane is defined by the Hessian
//! normal form: `n · x + d = 0`.

use crate::core::scalar::{Real, Point3r, Vector3r, TOLERANCE};

/// Classification of a point relative to a plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointClassification {
    /// On the positive (front) side of the plane.
    Front,
    /// On the negative (back) side of the plane.
    Back,
    /// Within tolerance of the plane.
    Coplanar,
}

/// An oriented plane in Hessian normal form: `normal · x + w = 0`.
#[derive(Clone, Copy, Debug)]
pub struct Plane {
    /// Unit normal vector.
    pub normal: Vector3r,
    /// Signed distance from origin (w = -normal · point_on_plane).
    pub w: Real,
}

impl Plane {
    /// Create a plane from a normal and signed distance.
    pub fn new(normal: Vector3r, w: Real) -> Self {
        Self { normal, w }
    }

    /// Create a plane from a normal and a point on the plane.
    pub fn from_normal_and_point(normal: Vector3r, point: &Point3r) -> Self {
        let n = normal.normalize();
        let w = -n.dot(&point.coords);
        Self { normal: n, w }
    }

    /// Create a plane from three non-collinear points (CCW winding → outward normal).
    pub fn from_three_points(a: &Point3r, b: &Point3r, c: &Point3r) -> Option<Self> {
        let ab = b - a;
        let ac = c - a;
        let cross = ab.cross(&ac);
        let len = cross.norm();
        if len < TOLERANCE {
            return None; // Degenerate (collinear points)
        }
        let normal = cross / len;
        let w = -normal.dot(&a.coords);
        Some(Self { normal, w })
    }

    /// Signed distance from a point to this plane.
    ///
    /// Positive = front side, negative = back side, ~0 = coplanar.
    #[inline]
    pub fn signed_distance(&self, point: &Point3r) -> Real {
        self.normal.dot(&point.coords) + self.w
    }

    /// Classify a point relative to this plane (using the default TOLERANCE).
    #[inline]
    pub fn classify_point(&self, point: &Point3r) -> PointClassification {
        self.classify_point_with_eps(point, TOLERANCE)
    }

    /// Classify a point with a custom epsilon (used by BSP operations).
    #[inline]
    pub fn classify_point_with_eps(&self, point: &Point3r, eps: Real) -> PointClassification {
        let dist = self.signed_distance(point);
        if dist > eps {
            PointClassification::Front
        } else if dist < -eps {
            PointClassification::Back
        } else {
            PointClassification::Coplanar
        }
    }

    /// Flip the plane (reverse normal and w).
    pub fn flip(&self) -> Self {
        Self {
            normal: -self.normal,
            w: -self.w,
        }
    }

    /// Compute the intersection parameter `t` along the line segment `a → b`.
    ///
    /// Returns `None` if the segment is parallel to the plane.
    pub fn intersect_segment(&self, a: &Point3r, b: &Point3r) -> Option<Real> {
        self.intersect_segment_with_eps(a, b, TOLERANCE)
    }

    /// Compute the intersection parameter with a custom epsilon.
    pub fn intersect_segment_with_eps(&self, a: &Point3r, b: &Point3r, eps: Real) -> Option<Real> {
        let da = self.signed_distance(a);
        let db = self.signed_distance(b);
        let denom = da - db;
        if denom.abs() < eps {
            return None;
        }
        let t = da / denom;
        // Clamp to [0,1] to avoid extrapolation from accumulated FP error
        Some(t.max(0.0).min(1.0))
    }
}

impl PartialEq for Plane {
    fn eq(&self, other: &Self) -> bool {
        (self.normal - other.normal).norm() < TOLERANCE
            && (self.w - other.w).abs() < TOLERANCE
    }
}
