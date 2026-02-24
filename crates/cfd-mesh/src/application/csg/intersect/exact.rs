//! Exact intersection predicates.

use crate::domain::core::scalar::Point3r;
use crate::domain::geometry::predicates::{orient_3d, Orientation};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

use super::types::IntersectionType;
use super::segment::compute_segment;

/// Test whether two triangles intersect using exact orientation predicates.
///
/// # Arguments
///
/// * `fa`, `pool_a` — Triangle and vertex pool from mesh A.
/// * `fb`, `pool_b` — Triangle and vertex pool from mesh B.  May equal `pool_a`.
///
/// # Returns
///
/// - [`IntersectionType::None`] when the triangles are disjoint.
/// - [`IntersectionType::Coplanar`] when they lie in the same plane.
/// - [`IntersectionType::Segment`] when they intersect along a segment.
pub fn intersect_triangles(
    fa: &FaceData,
    pool_a: &VertexPool,
    fb: &FaceData,
    pool_b: &VertexPool,
) -> IntersectionType {
    let a = pool_a.position(fa.vertices[0]);
    let b = pool_a.position(fa.vertices[1]);
    let c = pool_a.position(fa.vertices[2]);
    let d = pool_b.position(fb.vertices[0]);
    let e = pool_b.position(fb.vertices[1]);
    let f = pool_b.position(fb.vertices[2]);

    let arr = |p: &Point3r| [p.x, p.y, p.z];
    let dp = arr(d);
    let ep = arr(e);
    let fp = arr(f);

    // ── Step 1: classify T1 vertices against T2's plane ──────────────────
    let signs_t1 = [
        orient_3d(dp, ep, fp, arr(a)),
        orient_3d(dp, ep, fp, arr(b)),
        orient_3d(dp, ep, fp, arr(c)),
    ];

    if signs_t1.iter().all(|s| *s == Orientation::Degenerate) {
        return IntersectionType::Coplanar;
    }
    if not_straddling(&signs_t1) {
        return IntersectionType::None;
    }

    // ── Step 2: classify T2 vertices against T1's plane ──────────────────
    let ap = arr(a);
    let bp = arr(b);
    let cp = arr(c);
    let signs_t2 = [
        orient_3d(ap, bp, cp, arr(d)),
        orient_3d(ap, bp, cp, arr(e)),
        orient_3d(ap, bp, cp, arr(f)),
    ];

    if signs_t2.iter().all(|s| *s == Orientation::Degenerate) {
        return IntersectionType::Coplanar;
    }
    if not_straddling(&signs_t2) {
        return IntersectionType::None;
    }

    // ── Step 3: compute the intersection segment ──────────────────────────
    compute_segment(a, b, c, d, e, f)
}

/// Returns `true` when the triangle does **not** straddle the plane.
///
/// A triangle straddles a plane iff at least one vertex is `Positive` AND at
/// least one is `Negative`.  A `Degenerate` (on-plane) vertex counts as
/// neither side.
fn not_straddling(signs: &[Orientation; 3]) -> bool {
    let any_pos = signs.iter().any(|s| *s == Orientation::Positive);
    let any_neg = signs.iter().any(|s| *s == Orientation::Negative);
    !(any_pos && any_neg)
}
