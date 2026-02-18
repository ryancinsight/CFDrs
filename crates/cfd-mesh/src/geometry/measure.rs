//! Area and volume measurements.

use crate::core::scalar::{Real, Point3r};

/// Area of a triangle.
pub fn triangle_area(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    let ab = b - a;
    let ac = c - a;
    ab.cross(&ac).norm() * 0.5
}

/// Signed volume contribution of a triangle (for divergence theorem).
///
/// For a closed mesh, summing over all faces gives the enclosed volume.
/// Uses the formula: V_i = (1/6) * (v0 · (v1 × v2)).
pub fn signed_triangle_volume(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    a.coords.dot(&b.coords.cross(&c.coords)) / 6.0
}

/// Total surface area from an iterator of triangle vertex triples.
pub fn total_surface_area<'a>(
    triangles: impl Iterator<Item = (&'a Point3r, &'a Point3r, &'a Point3r)>,
) -> Real {
    triangles.map(|(a, b, c)| triangle_area(a, b, c)).sum()
}

/// Total signed volume from an iterator of triangle vertex triples.
///
/// For an outward-oriented closed mesh, this returns the positive enclosed volume.
pub fn total_signed_volume<'a>(
    triangles: impl Iterator<Item = (&'a Point3r, &'a Point3r, &'a Point3r)>,
) -> Real {
    triangles
        .map(|(a, b, c)| signed_triangle_volume(a, b, c))
        .sum()
}
