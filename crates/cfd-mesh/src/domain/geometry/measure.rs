//! Area and volume measurements — generic over scalar precision `T: Scalar`.

use crate::domain::core::scalar::Scalar;
use nalgebra::Point3;

/// Area of a triangle.
pub fn triangle_area<T: Scalar>(a: &Point3<T>, b: &Point3<T>, c: &Point3<T>) -> T {
    let ab = b - a;
    let ac = c - a;
    ab.cross(&ac).norm() * <T as Scalar>::from_f64(0.5)
}

/// Signed volume contribution of one triangle (divergence theorem).
///
/// Summed over all faces of a closed mesh, this yields the enclosed volume.
/// Formula: `V_i = (1/6) · v0 · (v1 × v2)`.
pub fn signed_triangle_volume<T: Scalar>(a: &Point3<T>, b: &Point3<T>, c: &Point3<T>) -> T {
    a.coords.dot(&b.coords.cross(&c.coords)) / <T as Scalar>::from_f64(6.0)
}

/// Total surface area from an iterator of triangle vertex triples.
pub fn total_surface_area<'a, T: Scalar + 'a>(
    triangles: impl Iterator<Item = (&'a Point3<T>, &'a Point3<T>, &'a Point3<T>)>,
) -> T {
    triangles
        .map(|(a, b, c)| triangle_area(a, b, c))
        .fold(T::zero(), |acc, v| acc + v)
}

/// Total signed volume from an iterator of triangle vertex triples.
///
/// Positive for outward-oriented closed meshes.
pub fn total_signed_volume<'a, T: Scalar + 'a>(
    triangles: impl Iterator<Item = (&'a Point3<T>, &'a Point3<T>, &'a Point3<T>)>,
) -> T {
    triangles
        .map(|(a, b, c)| signed_triangle_volume(a, b, c))
        .fold(T::zero(), |acc, v| acc + v)
}
