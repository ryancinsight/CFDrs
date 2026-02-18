//! Normal computation for triangles and polygon fans.

use crate::core::scalar::{Real, Point3r, Vector3r, TOLERANCE};

/// Compute the face normal of a triangle (CCW winding → outward).
///
/// Returns `None` if the triangle is degenerate.
pub fn triangle_normal(a: &Point3r, b: &Point3r, c: &Point3r) -> Option<Vector3r> {
    let ab = b - a;
    let ac = c - a;
    let cross = ab.cross(&ac);
    let len = cross.norm();
    if len < TOLERANCE {
        return None;
    }
    Some(cross / len)
}

/// Compute the area-weighted normal of a triangle (not normalized).
///
/// The magnitude is twice the triangle area.
pub fn triangle_area_normal(a: &Point3r, b: &Point3r, c: &Point3r) -> Vector3r {
    let ab = b - a;
    let ac = c - a;
    ab.cross(&ac)
}

/// Newell's method for computing the normal of a polygon with ≥3 vertices.
///
/// This is robust for non-planar polygons and handles concave cases.
/// Returns `None` if the polygon is degenerate.
pub fn newell_normal(vertices: &[Point3r]) -> Option<Vector3r> {
    if vertices.len() < 3 {
        return None;
    }
    let mut normal = Vector3r::zeros();
    let n = vertices.len();
    for i in 0..n {
        let curr = &vertices[i];
        let next = &vertices[(i + 1) % n];
        normal.x += (curr.y - next.y) * (curr.z + next.z);
        normal.y += (curr.z - next.z) * (curr.x + next.x);
        normal.z += (curr.x - next.x) * (curr.y + next.y);
    }
    let len = normal.norm();
    if len < TOLERANCE {
        return None;
    }
    Some(normal / len)
}

/// Compute vertex normal by averaging adjacent face normals.
///
/// `face_normals`: iterator of normals of faces incident to the vertex.
/// Returns the normalized average, or `None` if all normals cancel out.
pub fn average_normal<'a>(face_normals: impl Iterator<Item = &'a Vector3r>) -> Option<Vector3r> {
    let mut sum = Vector3r::zeros();
    let mut count = 0usize;
    for n in face_normals {
        sum += n;
        count += 1;
    }
    if count == 0 {
        return None;
    }
    let len = sum.norm();
    if len < TOLERANCE {
        return None;
    }
    Some(sum / len)
}

/// Angle-weighted vertex normal: weight each face normal by the angle
/// subtended at the vertex.
///
/// `faces`: iterator of `(face_normal, angle_at_vertex)` pairs.
pub fn angle_weighted_normal(faces: impl Iterator<Item = (Vector3r, Real)>) -> Option<Vector3r> {
    let mut sum = Vector3r::zeros();
    for (normal, angle) in faces {
        sum += normal * angle;
    }
    let len = sum.norm();
    if len < TOLERANCE {
        return None;
    }
    Some(sum / len)
}
