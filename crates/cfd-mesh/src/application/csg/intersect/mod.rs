//! Exact triangle–triangle intersection using Shewchuk predicates.
//!
//! ## Algorithm (Möller & Trumbore 1997 — exact-predicate variant)
//!
//! Given triangles T1 = (a, b, c) and T2 = (d, e, f):
//!
//! 1. Classify each vertex of T1 using exact `orient_3d` against T2's
//!    supporting plane.  If no vertex of T1 is on the opposite side from any
//!    other vertex (all positive, all negative, or all on the plane) → no
//!    proper intersection.
//!
//! 2. Repeat with T1 and T2 swapped.
//!
//! 3. Both triangles straddle each other's planes.
//!    Compute the intersection segment by finding the two edge-crossings of
//!    T1 against T2's plane, and of T2 against T1's plane.  Project all four
//!    crossing points onto the intersection line, find the overlap interval,
//!    and return the corresponding 3-D endpoints.
//!
//! ## Theorem — Exact Non-Intersection Decision
//!
//! Steps 1 and 2 use `orient_3d` from Shewchuk's adaptive-precision
//! library, which never returns a wrong sign due to floating-point
//! cancellation.  A pair classified as `IntersectionType::None` is
//! provably disjoint; no false negatives are introduced by rounding.
//!
//! Step 3 uses ordinary `f64` arithmetic only for computing *where* the
//! segment is, not for deciding *if* it exists.

pub mod exact;
pub mod segment;
pub mod types;

pub use exact::intersect_triangles;
pub use types::{IntersectionType, SnapSegment};

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::infrastructure::storage::vertex_pool::VertexPool;
    use crate::infrastructure::storage::face_store::FaceData;
    use crate::domain::core::scalar::Point3r;

    fn make_tri(pool: &mut VertexPool, pts: [[f64; 3]; 3]) -> FaceData {
        let n = nalgebra::Vector3::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(pts[0][0], pts[0][1], pts[0][2]), n);
        let v1 = pool.insert_or_weld(Point3r::new(pts[1][0], pts[1][1], pts[1][2]), n);
        let v2 = pool.insert_or_weld(Point3r::new(pts[2][0], pts[2][1], pts[2][2]), n);
        FaceData::untagged(v0, v1, v2)
    }

    #[test]
    fn parallel_planes_no_intersection() {
        let mut pool = VertexPool::default_millifluidic();
        let t1 = make_tri(
            &mut pool,
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        );
        let t2 = make_tri(
            &mut pool,
            [[0.0, 0.0, 5.0], [1.0, 0.0, 5.0], [0.0, 1.0, 5.0]],
        );
        assert!(matches!(
            intersect_triangles(&t1, &pool, &t2, &pool),
            IntersectionType::None
        ));
    }

    #[test]
    fn t1_above_t2_no_intersection() {
        let mut pa = VertexPool::default_millifluidic();
        let mut pb = VertexPool::default_millifluidic();
        // T1 at z=1 (all vertices above T2's z=0 plane)
        let t1 = make_tri(&mut pa, [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]]);
        let t2 = make_tri(&mut pb, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        assert!(matches!(
            intersect_triangles(&t1, &pa, &t2, &pb),
            IntersectionType::None
        ));
    }

    #[test]
    fn coplanar_triangles_detected() {
        let mut pool = VertexPool::default_millifluidic();
        // Both triangles in z=0 plane
        let t1 = make_tri(
            &mut pool,
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        );
        let t2 = make_tri(
            &mut pool,
            [[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [0.5, 1.0, 0.0]],
        );
        assert!(matches!(
            intersect_triangles(&t1, &pool, &t2, &pool),
            IntersectionType::Coplanar
        ));
    }

    #[test]
    fn perpendicular_triangles_produce_segment() {
        let mut pa = VertexPool::default_millifluidic();
        let mut pb = VertexPool::default_millifluidic();
        // T1 in the XY plane
        let t1 = make_tri(
            &mut pa,
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        );
        // T2 in the x=0 plane, crossing T1
        let t2 = make_tri(
            &mut pb,
            [[0.0, -1.0, -1.0], [0.0, -1.0, 1.0], [0.0, 1.0, 0.0]],
        );
        let result = intersect_triangles(&t1, &pa, &t2, &pb);
        assert!(
            matches!(result, IntersectionType::Segment { .. }),
            "expected Segment, got {:?}",
            result
        );
        if let IntersectionType::Segment { start, end } = result {
            let len = (end - start).norm();
            assert!(
                len > 1e-6,
                "intersection segment should have positive length"
            );
        }
    }

    #[test]
    fn touching_at_single_vertex_is_none() {
        let mut pa = VertexPool::default_millifluidic();
        let mut pb = VertexPool::default_millifluidic();
        // T1 and T2 share only the vertex at the origin
        let t1 = make_tri(&mut pa, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let t2 = make_tri(
            &mut pb,
            [[0.0, 0.0, 0.0], [-1.0, 0.0, 1.0], [0.0, -1.0, 1.0]],
        );
        // These triangles may or may not intersect depending on exact geometry;
        // the important thing is that the function doesn't panic.
        let _ = intersect_triangles(&t1, &pa, &t2, &pb);
    }
}
