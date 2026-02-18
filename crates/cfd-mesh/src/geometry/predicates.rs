//! # Exact Geometric Predicates
//!
//! Thin wrappers around [Shewchuk's adaptive-precision arithmetic][shewchuk]
//! as provided by the [`geometry-predicates`] crate.  All functions in this
//! module run in **exact arithmetic** — they never return a wrong sign due to
//! floating-point rounding, even for nearly-degenerate configurations.
//!
//! ## Theorem — Shewchuk Adaptive Arithmetic
//!
//! Let `fl(e)` denote the floating-point evaluation of an expression `e`.
//! Shewchuk's method computes `sign(e)` exactly by maintaining an *expansion*
//! — a non-overlapping sequence of floating-point numbers whose sum equals `e`
//! exactly.  The adaptive stages refine the expansion only until the sign is
//! determined, achieving near-floating-point speed for well-separated inputs
//! while falling back to full multi-precision for near-degenerate cases.
//!
//! ## Available Predicates
//!
//! | Function | Returns | Meaning |
//! |---|---|---|
//! | [`orient_2d`] | [`Orientation`] | Sign of the 2-D cross product |
//! | [`orient_3d`] | [`Orientation`] | Sign of the 3×3 tet volume determinant |
//! | [`incircle`] | [`Orientation`] | Point inside/on/outside circumcircle |
//! | [`insphere`] | [`Orientation`] | Point inside/on/outside circumsphere |
//!
//! ## Diagram
//!
//! ```text
//! orient_2d(a, b, c):
//!
//!      c
//!     /
//!    /  CCW (+) → Positive
//!   a ──── b
//!
//!   det = | ax  ay  1 |
//!         | bx  by  1 |  > 0 ↔ CCW
//!         | cx  cy  1 |
//! ```
//!
//! [shewchuk]: https://www.cs.cmu.edu/~quake/robust.html

use geometry_predicates as gp;
use nalgebra::Point2;

use crate::core::scalar::Real;

/// The sign of an orientation determinant.
///
/// Returned by all predicate functions.  The `Degenerate` variant corresponds
/// to an exact zero determinant (collinear / coplanar / co-spherical inputs).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orientation {
    /// Determinant is strictly positive (CCW / above / inside).
    Positive,
    /// Determinant is exactly zero (degenerate configuration).
    Degenerate,
    /// Determinant is strictly negative (CW / below / outside).
    Negative,
}

impl Orientation {
    /// Convert a raw determinant value to an `Orientation`.
    #[inline]
    pub fn from_det(d: f64) -> Self {
        if d > 0.0 {
            Orientation::Positive
        } else if d < 0.0 {
            Orientation::Negative
        } else {
            Orientation::Degenerate
        }
    }

    /// Returns `true` if the orientation is [`Positive`](Orientation::Positive).
    #[inline]
    pub fn is_positive(self) -> bool {
        self == Orientation::Positive
    }

    /// Returns `true` if the orientation is [`Negative`](Orientation::Negative).
    #[inline]
    pub fn is_negative(self) -> bool {
        self == Orientation::Negative
    }

    /// Returns `true` for a degenerate (zero) determinant.
    #[inline]
    pub fn is_degenerate(self) -> bool {
        self == Orientation::Degenerate
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Convert a [`Real`] to `f64` for the predicates crate.
///
/// `geometry-predicates` always operates in `f64`.  If the crate is compiled
/// with `f32`, we upcast — this is safe because `f32` is a strict subset of
/// `f64`'s representable values.
#[inline]
fn r(v: Real) -> f64 {
    v as f64
}

// ── 2-D predicates ────────────────────────────────────────────────────────────

/// **Exact 2-D orientation test.**
///
/// Returns the sign of the 2×2 determinant:
///
/// ```text
/// | ax - cx   ay - cy |
/// | bx - cx   by - cy |
/// ```
///
/// - [`Positive`][Orientation::Positive] — `a`, `b`, `c` are in CCW order.
/// - [`Negative`][Orientation::Negative] — `a`, `b`, `c` are in CW order.
/// - [`Degenerate`][Orientation::Degenerate] — `a`, `b`, `c` are collinear.
///
/// # Example
/// ```rust
/// use cfd_mesh::geometry::predicates::{orient_2d, Orientation};
/// use nalgebra::Point2;
///
/// let a = Point2::new(0.0_f64, 0.0);
/// let b = Point2::new(1.0, 0.0);
/// let c = Point2::new(0.0, 1.0);
/// assert_eq!(orient_2d(&a, &b, &c), Orientation::Positive); // CCW
/// ```
pub fn orient_2d(a: &Point2<Real>, b: &Point2<Real>, c: &Point2<Real>) -> Orientation {
    let det = gp::orient2d(
        [r(a.x), r(a.y)],
        [r(b.x), r(b.y)],
        [r(c.x), r(c.y)],
    );
    Orientation::from_det(det)
}

/// **Exact 2-D orientation test from raw arrays.**
///
/// Convenience overload accepting `[Real; 2]` arrays.
#[inline]
pub fn orient_2d_arr(a: [Real; 2], b: [Real; 2], c: [Real; 2]) -> Orientation {
    let det = gp::orient2d(
        [r(a[0]), r(a[1])],
        [r(b[0]), r(b[1])],
        [r(c[0]), r(c[1])],
    );
    Orientation::from_det(det)
}

// ── 3-D predicates ────────────────────────────────────────────────────────────

/// **Exact 3-D orientation test.**
///
/// Returns the sign of the **signed volume** of the tetrahedron `(a, b, c, d)`:
///
/// - [`Positive`][Orientation::Positive] — `d` lies *above* the plane `abc`
///   when `a→b→c` is counter-clockwise (right-hand rule); positive signed volume.
/// - [`Negative`][Orientation::Negative] — `d` lies *below* the plane.
/// - [`Degenerate`][Orientation::Degenerate] — all four points are coplanar.
///
/// Used in the BVH Boolean narrow phase to determine on which side of a
/// supporting plane a vertex lies.
///
/// # Example
/// ```rust
/// use cfd_mesh::geometry::predicates::{orient_3d, Orientation};
///
/// // Tetrahedron with positive volume
/// let o  = [0.0_f64, 0.0, 0.0];
/// let ex = [1.0, 0.0, 0.0];
/// let ey = [0.0, 1.0, 0.0];
/// let ez = [0.0, 0.0, 1.0];
/// assert_eq!(orient_3d(o, ex, ey, ez), Orientation::Positive);
/// ```
pub fn orient_3d(
    a: [Real; 3],
    b: [Real; 3],
    c: [Real; 3],
    d: [Real; 3],
) -> Orientation {
    let to64 = |v: [Real; 3]| [r(v[0]), r(v[1]), r(v[2])];
    // gp::orient3d returns positive when d is *below* the abc plane (Shewchuk
    // convention).  We negate so that our API convention is "d above = Positive",
    // which matches standard signed-volume / right-hand-rule intuition.
    let det = -gp::orient3d(to64(a), to64(b), to64(c), to64(d));
    Orientation::from_det(det)
}

/// Convenience wrapper: accept `nalgebra::Point3<Real>`.
pub fn orient_3d_pts(
    a: &nalgebra::Point3<Real>,
    b: &nalgebra::Point3<Real>,
    c: &nalgebra::Point3<Real>,
    d: &nalgebra::Point3<Real>,
) -> Orientation {
    orient_3d(
        [a.x, a.y, a.z],
        [b.x, b.y, b.z],
        [c.x, c.y, c.z],
        [d.x, d.y, d.z],
    )
}

// ── In-circle / In-sphere predicates ─────────────────────────────────────────

/// **Exact in-circle test.**
///
/// Returns the sign of the 3×3 determinant that determines whether point `d`
/// lies inside, on, or outside the circumcircle of triangle `abc` (with `abc`
/// in CCW order):
///
/// ```text
/// | ax - dx   ay - dy   (ax²+ay²) - (dx²+dy²) |
/// | bx - dx   by - dy   (bx²+by²) - (dx²+dy²) |
/// | cx - dx   cy - dy   (cx²+cy²) - (dx²+dy²) |
/// ```
///
/// - [`Positive`][Orientation::Positive] — `d` is strictly inside the circumcircle.
/// - [`Negative`][Orientation::Negative] — `d` is strictly outside.
/// - [`Degenerate`][Orientation::Degenerate] — `d` lies exactly on the circle.
///
/// Used in Delaunay mesh refinement to enforce the empty-circumcircle property.
pub fn incircle(
    a: &Point2<Real>,
    b: &Point2<Real>,
    c: &Point2<Real>,
    d: &Point2<Real>,
) -> Orientation {
    let det = gp::incircle(
        [r(a.x), r(a.y)],
        [r(b.x), r(b.y)],
        [r(c.x), r(c.y)],
        [r(d.x), r(d.y)],
    );
    Orientation::from_det(det)
}

/// **Exact in-sphere test.**
///
/// Returns the sign of the 4×4 determinant that determines whether point `e`
/// lies inside, on, or outside the circumsphere of tetrahedron `abcd` (with
/// `abcd` in positive orientation):
///
/// - [`Positive`][Orientation::Positive] — `e` is strictly inside the sphere.
/// - [`Negative`][Orientation::Negative] — `e` is strictly outside.
/// - [`Degenerate`][Orientation::Degenerate] — `e` lies on the sphere.
///
/// Used in 3-D Delaunay mesh generation to enforce the Delaunay property.
pub fn insphere(
    a: [Real; 3],
    b: [Real; 3],
    c: [Real; 3],
    d: [Real; 3],
    e: [Real; 3],
) -> Orientation {
    let to64 = |v: [Real; 3]| [r(v[0]), r(v[1]), r(v[2])];
    let det = gp::insphere(to64(a), to64(b), to64(c), to64(d), to64(e));
    Orientation::from_det(det)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point2;

    // ── orient_2d ──────────────────────────────────────────────────────────

    #[test]
    fn orient_2d_ccw() {
        // Standard CCW triangle
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.0, 1.0);
        assert_eq!(orient_2d(&a, &b, &c), Orientation::Positive);
    }

    #[test]
    fn orient_2d_cw() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(0.0, 1.0);
        let c = Point2::new(1.0, 0.0);
        assert_eq!(orient_2d(&a, &b, &c), Orientation::Negative);
    }

    #[test]
    fn orient_2d_collinear() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(2.0, 0.0);
        assert_eq!(orient_2d(&a, &b, &c), Orientation::Degenerate);
    }

    // ── orient_3d ──────────────────────────────────────────────────────────

    #[test]
    fn orient_3d_positive_tet() {
        let o  = [0.0, 0.0, 0.0];
        let ex = [1.0, 0.0, 0.0];
        let ey = [0.0, 1.0, 0.0];
        let ez = [0.0, 0.0, 1.0];
        assert_eq!(orient_3d(o, ex, ey, ez), Orientation::Positive);
    }

    #[test]
    fn orient_3d_negative_tet() {
        let o  = [0.0, 0.0, 0.0];
        let ex = [1.0, 0.0, 0.0];
        let ey = [0.0, 1.0, 0.0];
        // Swap ex and ey to flip sign
        assert_eq!(orient_3d(o, ey, ex, [0.0, 0.0, 1.0]), Orientation::Negative);
    }

    #[test]
    fn orient_3d_coplanar() {
        // All four points in z=0 plane
        assert_eq!(
            orient_3d([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]),
            Orientation::Degenerate
        );
    }

    // ── incircle ───────────────────────────────────────────────────────────

    #[test]
    fn incircle_inside() {
        // Unit circle: points on the circle, test point at origin (inside)
        let a = Point2::new(1.0, 0.0);
        let b = Point2::new(0.0, 1.0);
        let c = Point2::new(-1.0, 0.0);
        let d = Point2::new(0.0, 0.0); // origin — strictly inside
        assert_eq!(incircle(&a, &b, &c, &d), Orientation::Positive);
    }

    #[test]
    fn incircle_outside() {
        let a = Point2::new(1.0, 0.0);
        let b = Point2::new(0.0, 1.0);
        let c = Point2::new(-1.0, 0.0);
        let d = Point2::new(0.0, 2.0); // above the circle — outside
        assert_eq!(incircle(&a, &b, &c, &d), Orientation::Negative);
    }

    // ── insphere ───────────────────────────────────────────────────────────

    #[test]
    fn insphere_inside() {
        // Unit sphere vertices, test point at origin
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = [0.0, 0.0, 1.0];
        let d = [-1.0, 0.0, 0.0];
        let e = [0.0, 0.0, 0.0]; // origin — inside unit sphere
        assert_eq!(insphere(a, b, c, d, e), Orientation::Positive);
    }

    #[test]
    fn insphere_outside() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = [0.0, 0.0, 1.0];
        let d = [-1.0, 0.0, 0.0];
        let e = [2.0, 0.0, 0.0]; // outside unit sphere
        assert_eq!(insphere(a, b, c, d, e), Orientation::Negative);
    }

    // ── Orientation helpers ────────────────────────────────────────────────

    #[test]
    fn orientation_helpers() {
        assert!(Orientation::Positive.is_positive());
        assert!(!Orientation::Positive.is_negative());
        assert!(!Orientation::Positive.is_degenerate());

        assert!(Orientation::Negative.is_negative());
        assert!(Orientation::Degenerate.is_degenerate());
    }
}
