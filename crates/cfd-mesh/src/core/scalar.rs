//! Scalar type abstraction — zero-cost generic floating-point precision.
//!
//! # Design
//!
//! Instead of a compile-time feature flag that forces a single precision
//! across the whole crate, every mesh type is generic over `T: Scalar`.
//! Monomorphisation generates optimal machine code per instantiation —
//! identical to a hand-written `f64`-only implementation — while letting
//! callers freely mix `IndexedMesh<f32>` and `IndexedMesh<f64>` in one
//! binary without recompilation.
//!
//! # Theorem: Scalar completeness
//!
//! **Statement**: `Scalar` covers exactly the floating-point types that
//! support all mesh-geometry operations required by `cfd-mesh`.
//!
//! **Proof sketch**: `nalgebra::RealField` provides a complete ordered field
//! with algebraic/trigonometric operations needed for nalgebra vector maths.
//! `num_traits::Float` adds `infinity`, `neg_infinity`, `is_finite`, `floor`,
//! `sqrt`, and `min`/`max`.  `num_traits::ToPrimitive` enables lossless
//! conversion to `f64` for human-readable outputs.  The sealed super-trait
//! restricts the impl set to `{f32, f64}`, matching IEEE 754 hardware support.

use nalgebra::{Point3, Vector3, Matrix4};

// ── Sealed trait ──────────────────────────────────────────────────────────────
// Prevents downstream crates from implementing `Scalar` for arbitrary types.
mod private {
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

/// Zero-cost generic floating-point scalar for all CFD mesh operations.
///
/// Implemented only by `f32` and `f64`.  Mesh types parameterised by
/// `T: Scalar` monomorphise to zero-overhead code, equivalent to writing a
/// separate `f32` and `f64` implementation by hand.
///
/// # Choosing between precisions
///
/// | Precision | Tolerance | Use case |
/// |-----------|-----------|----------|
/// | `f64` (default) | 1 nm | High-fidelity CFD, validation, export to OpenFOAM |
/// | `f32` | 10 µm | GPU-side geometry staging where bandwidth matters |
///
/// # Example
///
/// ```rust,ignore
/// // Both coexist in the same binary — no feature flag, no recompilation:
/// let hi: IndexedMesh<f64> = IndexedMesh::new();
/// let lo: IndexedMesh<f32> = IndexedMesh::new();
/// ```
pub trait Scalar:
    nalgebra::RealField
    + num_traits::Float
    + num_traits::ToPrimitive
    + Copy
    + Default
    + std::fmt::Debug
    + std::fmt::Display
    + Send
    + Sync
    + 'static
    + private::Sealed
{
    /// Absolute geometry tolerance appropriate for this precision.
    ///
    /// - `f64` → `1 × 10⁻⁹` m  (sub-nanometer; millifluidic mm-scale geometry)
    /// - `f32` → `1 × 10⁻⁵` m  (10 µm; bounded by single-precision rounding)
    fn tolerance() -> Self;

    /// Convert an `f64` literal to this scalar type.
    ///
    /// Zero-cost identity for `f64`; one `as` cast for `f32`.
    /// Enables generic code to write `T::from_f64(0.5)` instead of `0.5_T`.
    fn from_f64(v: f64) -> Self;

    /// Squared tolerance — avoids `sqrt` in distance comparisons.
    #[inline]
    fn tolerance_sq() -> Self {
        let t = Self::tolerance();
        t * t
    }
}

impl Scalar for f64 {
    #[inline(always)]
    fn tolerance() -> Self { 1e-9 }
    #[inline(always)]
    fn from_f64(v: f64) -> Self { v }
}

impl Scalar for f32 {
    #[inline(always)]
    fn tolerance() -> Self { 1e-5_f32 }
    #[inline(always)]
    fn from_f64(v: f64) -> Self { v as f32 }
}

// ── Default-precision convenience aliases ─────────────────────────────────────

/// Default scalar precision — `f64` for sub-nanometer millifluidic accuracy.
///
/// All existing code that references `Real` continues to compile unchanged.
/// New code should prefer the generic `T: Scalar` pattern.
pub type Real = f64;

/// 3-D point at default (`f64`) precision.
pub type Point3r = Point3<Real>;

/// 3-D vector at default (`f64`) precision.
pub type Vector3r = Vector3<Real>;

/// 4×4 homogeneous transform at default (`f64`) precision.
pub type Matrix4r = Matrix4<Real>;

/// Absolute geometry tolerance at default precision (1 nm).
pub const TOLERANCE: Real = 1e-9;

/// Squared tolerance at default precision — avoids `sqrt` in distance checks.
pub const TOLERANCE_SQ: Real = TOLERANCE * TOLERANCE;

// ── Generic sanitisation helpers ──────────────────────────────────────────────

/// Replace NaN or ±Inf with zero — generic over any `T: Scalar`.
#[inline]
pub fn sanitize<T: Scalar>(v: T) -> T {
    if <T as num_traits::Float>::is_finite(v) { v } else { T::zero() }
}

/// Replace NaN / ±Inf components of a point with zero — generic.
#[inline]
pub fn sanitize_point<T: Scalar>(p: &Point3<T>) -> Point3<T> {
    Point3::new(sanitize(p.x), sanitize(p.y), sanitize(p.z))
}

/// Replace NaN / ±Inf components of a vector with zero — generic.
#[inline]
pub fn sanitize_vector<T: Scalar>(v: &Vector3<T>) -> Vector3<T> {
    Vector3::new(sanitize(v.x), sanitize(v.y), sanitize(v.z))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerance_ordering() {
        assert!(f32::tolerance() > f64::tolerance() as f32,
            "f32 tolerance must be coarser than f64");
    }

    #[test]
    fn from_f64_identity_f64() {
        assert_eq!(f64::from_f64(1.0_f64), 1.0_f64);
    }

    #[test]
    fn from_f64_cast_f32() {
        let v: f32 = f32::from_f64(0.5_f64);
        assert!((v - 0.5_f32).abs() < 1e-7, "f32 cast must be accurate");
    }

    #[test]
    fn sanitize_finite_passthrough() {
        assert_eq!(sanitize(1.5_f64), 1.5_f64);
        assert_eq!(sanitize(1.5_f32), 1.5_f32);
    }

    #[test]
    fn sanitize_nan_to_zero() {
        assert_eq!(sanitize(f64::NAN), 0.0_f64);
        assert_eq!(sanitize(f32::NAN), 0.0_f32);
    }

    #[test]
    fn sanitize_inf_to_zero() {
        assert_eq!(sanitize(f64::INFINITY), 0.0_f64);
        assert_eq!(sanitize(f32::NEG_INFINITY), 0.0_f32);
    }
}
