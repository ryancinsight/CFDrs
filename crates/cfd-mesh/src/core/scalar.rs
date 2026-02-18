//! Scalar type abstraction — compile-time selection of `f64` or `f32`.

use nalgebra::{Point3, Vector3, Matrix4};

/// The floating-point type used throughout the crate.
#[cfg(feature = "f64")]
pub type Real = f64;

/// The floating-point type used throughout the crate.
#[cfg(feature = "f32")]
pub type Real = f32;

/// 3D point alias.
pub type Point3r = Point3<Real>;

/// 3D vector alias.
pub type Vector3r = Vector3<Real>;

/// 4×4 transform matrix alias.
pub type Matrix4r = Matrix4<Real>;

/// Tolerance for floating-point comparisons.
///
/// Defaults based on precision:
/// - `f64`: 1e-9 (sub-nanometer for millifluidic mm-scale)
/// - `f32`: 1e-5
#[cfg(feature = "f64")]
pub const TOLERANCE: Real = 1e-9;

/// Tolerance for floating-point comparisons.
#[cfg(feature = "f32")]
pub const TOLERANCE: Real = 1e-5;

/// Squared tolerance (avoids sqrt in distance checks).
pub const TOLERANCE_SQ: Real = TOLERANCE * TOLERANCE;

/// Sanitize a scalar: replace NaN/Inf with zero.
#[inline]
pub fn sanitize(v: Real) -> Real {
    if v.is_finite() { v } else { 0.0 as Real }
}

/// Sanitize a point: replace NaN/Inf components with zero.
#[inline]
pub fn sanitize_point(p: &Point3r) -> Point3r {
    Point3r::new(sanitize(p.x), sanitize(p.y), sanitize(p.z))
}

/// Sanitize a vector: replace NaN/Inf components with zero.
#[inline]
pub fn sanitize_vector(v: &Vector3r) -> Vector3r {
    Vector3r::new(sanitize(v.x), sanitize(v.y), sanitize(v.z))
}
