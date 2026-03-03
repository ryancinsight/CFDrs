//! Shared exact 3-D geometric predicates for CSG pipelines.
//!
//! ## Algorithm
//!
//! This module centralizes exact helper predicates used by arrangement,
//! coplanar, and seam-propagation stages:
//!
//! 1. `triangle_is_degenerate_exact`:
//!    exact 3-D triangle degeneracy via three projected `orient_2d_arr` tests.
//! 2. `collinear_3d_exact`:
//!    exact collinearity of point `p` with segment endpoints `(a,b)`.
//! 3. `point_on_segment_exact`:
//!    strict interior point-on-segment (`0 < t < 1`) after exact collinearity.
//!
//! ## Theorem — Projection-Minor Equivalence in 3-D
//!
//! Let vectors `u = b-a` and `v = p-a`.
//! `u` and `v` are linearly dependent in `R^3` iff all three 2x2 minors
//! (`xy`, `xz`, `yz`) vanish. These minors are exactly the determinants tested
//! by `orient_2d_arr` on projected coordinates.
//!
//! **Proof sketch.**
//! Linear dependence in `R^3` is equivalent to rank 1 of matrix `[u v]`.
//! Rank 1 iff every 2x2 minor is zero. Each projected orientation determinant
//! computes one of these minors exactly (Shewchuk predicates). ∎

use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};

/// Exact 3-D triangle degeneracy check.
#[inline]
pub(crate) fn triangle_is_degenerate_exact(a: &Point3r, b: &Point3r, c: &Point3r) -> bool {
    orient_2d_arr([a.x, a.y], [b.x, b.y], [c.x, c.y]) == Orientation::Degenerate
        && orient_2d_arr([a.x, a.z], [b.x, b.z], [c.x, c.z]) == Orientation::Degenerate
        && orient_2d_arr([a.y, a.z], [b.y, b.z], [c.y, c.z]) == Orientation::Degenerate
}

/// Exact 3-D collinearity check for point `p` with segment `(a,b)`.
#[inline]
pub(crate) fn collinear_3d_exact(a: &Point3r, b: &Point3r, p: &Point3r) -> bool {
    orient_2d_arr([a.x, a.y], [b.x, b.y], [p.x, p.y]) == Orientation::Degenerate
        && orient_2d_arr([a.x, a.z], [b.x, b.z], [p.x, p.z]) == Orientation::Degenerate
        && orient_2d_arr([a.y, a.z], [b.y, b.z], [p.y, p.z]) == Orientation::Degenerate
}

/// Exact strict-interior point-on-segment parameter.
///
/// Returns `Some(t)` only when `p = a + t*(b-a)` with `0 < t < 1`.
#[inline]
pub(crate) fn point_on_segment_exact(a: &Point3r, b: &Point3r, p: &Point3r) -> Option<Real> {
    if !collinear_3d_exact(a, b, p) {
        return None;
    }
    let edge = *b - *a;
    let edge_len_sq = edge.dot(&edge);
    if edge_len_sq <= 0.0 {
        return None;
    }
    let t = (*p - *a).dot(&edge) / edge_len_sq;
    if t > 0.0 && t < 1.0 {
        Some(t)
    } else {
        None
    }
}

#[derive(Clone, Copy)]
enum DropAxis {
    X,
    Y,
    Z,
}

#[inline]
fn choose_drop_axis(n: &Vector3r) -> Option<DropAxis> {
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    let max_c = ax.max(ay).max(az);
    if max_c <= 0.0 {
        return None;
    }
    if ax >= ay && ax >= az {
        Some(DropAxis::X)
    } else if ay >= az {
        Some(DropAxis::Y)
    } else {
        Some(DropAxis::Z)
    }
}

#[inline]
fn project_drop_axis(p: &Point3r, axis: DropAxis) -> [Real; 2] {
    match axis {
        DropAxis::X => [p.y, p.z],
        DropAxis::Y => [p.x, p.z],
        DropAxis::Z => [p.x, p.y],
    }
}

/// Proper interior segment intersection via exact projected orientation tests.
///
/// Returns parameters `(t, u)` such that:
/// - `ab(t) = a + t*(b-a)`
/// - `cd(u) = c + u*(d-c)`
/// - `0 < t < 1`, `0 < u < 1` (proper interior intersection)
///
/// The projection plane is chosen by dropping the dominant component of
/// `plane_n`, maximizing projected area and numeric conditioning.
///
/// # Theorem — Dominant-Axis Projection Preserves Coplanar Segment Ordering
///
/// For coplanar points on plane normal `n != 0`, dropping the dominant axis
/// of `n` defines a one-to-one affine map from that plane to 2-D. Affine maps
/// preserve orientation signs up to a global sign, so exact 2-D orientation
/// predicates on projected points correctly determine segment straddling.
///
/// **Proof sketch.**
/// Restrict projection to the plane; its Jacobian determinant equals a nonzero
/// component of `n` (by dominant-axis choice), so the map is invertible on the
/// plane. Invertible affine transforms preserve incidence and orientation sign
/// consistency, thus proper 2-D straddle implies and is implied by proper 3-D
/// coplanar segment intersection. ∎
pub(crate) fn proper_segment_intersection_params_projected_exact(
    a: &Point3r,
    b: &Point3r,
    c: &Point3r,
    d: &Point3r,
    plane_n: &Vector3r,
) -> Option<(Real, Real)> {
    let axis = choose_drop_axis(plane_n)?;
    let a2 = project_drop_axis(a, axis);
    let b2 = project_drop_axis(b, axis);
    let c2 = project_drop_axis(c, axis);
    let d2 = project_drop_axis(d, axis);

    let o1 = orient_2d_arr(a2, b2, c2);
    let o2 = orient_2d_arr(a2, b2, d2);
    let o3 = orient_2d_arr(c2, d2, a2);
    let o4 = orient_2d_arr(c2, d2, b2);

    // Proper interior intersection only: no collinear/endpoint-touch cases.
    if o1 == Orientation::Degenerate
        || o2 == Orientation::Degenerate
        || o3 == Orientation::Degenerate
        || o4 == Orientation::Degenerate
    {
        return None;
    }
    let ab_straddle = o1 != o2;
    let cd_straddle = o3 != o4;
    if !ab_straddle || !cd_straddle {
        return None;
    }

    let r = [b2[0] - a2[0], b2[1] - a2[1]];
    let s = [d2[0] - c2[0], d2[1] - c2[1]];
    let qp = [c2[0] - a2[0], c2[1] - a2[1]];
    let denom = r[0] * s[1] - r[1] * s[0];
    if denom == 0.0 {
        return None;
    }
    let t = (qp[0] * s[1] - qp[1] * s[0]) / denom;
    let u = (qp[0] * r[1] - qp[1] * r[0]) / denom;
    Some((t, u))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: Real, y: Real, z: Real) -> Point3r {
        Point3r::new(x, y, z)
    }

    #[test]
    fn triangle_is_degenerate_exact_detects_collinear_triangle() {
        assert!(triangle_is_degenerate_exact(
            &p(0.0, 0.0, 0.0),
            &p(1.0, 0.0, 0.0),
            &p(2.0, 0.0, 0.0)
        ));
    }

    #[test]
    fn point_on_segment_exact_accepts_true_collinear_midpoint() {
        let t = point_on_segment_exact(&p(0.0, 0.0, 0.0), &p(1.0, 2.0, 3.0), &p(0.5, 1.0, 1.5));
        assert!(t.is_some());
    }

    #[test]
    fn point_on_segment_exact_rejects_projection_false_positive() {
        let t = point_on_segment_exact(&p(0.0, 0.0, 0.0), &p(1.0, 1.0, 1.0), &p(0.5, 0.5, 0.6));
        assert!(t.is_none());
    }

    #[test]
    fn proper_segment_intersection_detects_near_parallel_crossing() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 1.0e-10, 0.0);
        let c = p(0.5, -1.0e-10, 0.0);
        let d = p(0.5, 1.0e-10, 0.0);
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let hit = proper_segment_intersection_params_projected_exact(&a, &b, &c, &d, &n);
        assert!(
            hit.is_some(),
            "near-parallel but proper crossing must be detected"
        );
        let (t, u) = hit.expect("intersection params");
        assert!(t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0);
    }

    #[test]
    fn proper_segment_intersection_rejects_endpoint_touch() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 0.0, 1.0);
        let c = p(1.0, 0.0, 1.0);
        let d = p(1.0, 1.0, 1.0);
        let n = Vector3r::new(1.0, 0.0, -1.0);
        let hit = proper_segment_intersection_params_projected_exact(&a, &b, &c, &d, &n);
        assert!(
            hit.is_none(),
            "endpoint-touch must not be classified as proper interior crossing"
        );
    }

    #[test]
    fn proper_segment_intersection_detects_tilted_plane_crossing() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 0.0, 1.0);
        let c = p(0.5, -1.0, 0.5);
        let d = p(0.5, 1.0, 0.5);
        let n = Vector3r::new(1.0, 0.0, -1.0);
        let hit = proper_segment_intersection_params_projected_exact(&a, &b, &c, &d, &n);
        assert!(hit.is_some(), "crossing in tilted plane must be detected");
        let (t, u) = hit.expect("intersection params");
        assert!((t - 0.5).abs() < 1e-12);
        assert!((u - 0.5).abs() < 1e-12);
    }
}
