//! Distance queries between geometric primitives.

use crate::core::scalar::{Real, Point3r, Vector3r};
use crate::geometry::plane::Plane;

/// Closest point on a line segment `a–b` to point `p`.
pub fn closest_point_on_segment(p: &Point3r, a: &Point3r, b: &Point3r) -> Point3r {
    let ab = b - a;
    let ap = p - a;
    let len_sq = ab.norm_squared();
    if len_sq < Real::EPSILON {
        return *a;
    }
    let t = (ap.dot(&ab) / len_sq).clamp(0.0, 1.0);
    Point3r::from(a.coords + ab * t)
}

/// Distance from a point to a line segment.
pub fn point_to_segment_distance(p: &Point3r, a: &Point3r, b: &Point3r) -> Real {
    (p - closest_point_on_segment(p, a, b)).norm()
}

/// Closest point on a triangle `(a, b, c)` to point `p`.
///
/// Uses the Barycentric projection method.
pub fn closest_point_on_triangle(
    p: &Point3r,
    a: &Point3r,
    b: &Point3r,
    c: &Point3r,
) -> Point3r {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;

    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return *a;
    }

    let bp = p - b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return *b;
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return Point3r::from(a.coords + ab * v);
    }

    let cp = p - c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return *c;
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return Point3r::from(a.coords + ac * w);
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return Point3r::from(b.coords + (c - b) * w);
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    Point3r::from(a.coords + ab * v + ac * w)
}

/// Distance from a point to a plane.
pub fn point_to_plane_distance(p: &Point3r, plane: &Plane) -> Real {
    plane.signed_distance(p).abs()
}

/// Project a point onto a plane.
pub fn project_point_to_plane(p: &Point3r, plane: &Plane) -> Point3r {
    let dist = plane.signed_distance(p);
    Point3r::from(p.coords - plane.normal * dist)
}

/// Ray-triangle intersection (Möller–Trumbore algorithm).
///
/// Returns `Some(t)` where `t` is the parameter along the ray `origin + t * dir`.
/// Returns `None` if no intersection.
pub fn ray_triangle_intersection(
    origin: &Point3r,
    dir: &Vector3r,
    a: &Point3r,
    b: &Point3r,
    c: &Point3r,
) -> Option<Real> {
    let edge1 = b - a;
    let edge2 = c - a;
    let h = dir.cross(&edge2);
    let det = edge1.dot(&h);

    if det.abs() < Real::EPSILON {
        return None; // Parallel
    }

    let inv_det = 1.0 / det;
    let s = origin - a;
    let u = inv_det * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = s.cross(&edge1);
    let v = inv_det * dir.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = inv_det * edge2.dot(&q);
    if t > Real::EPSILON {
        Some(t)
    } else {
        None
    }
}
