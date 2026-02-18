//! Polygon classification against a splitting plane.

use crate::core::index::VertexId;
use crate::core::scalar::Point3r;
use crate::geometry::plane::{Plane, PointClassification};
use crate::geometry::predicates::{orient_3d, Orientation};
use crate::storage::vertex_pool::VertexPool;
use crate::csg::bsp::BSP_PLANE_TOLERANCE;

/// Classification of a polygon relative to a plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolygonClassification {
    /// All vertices are in front of (or coplanar with) the plane.
    Front,
    /// All vertices are behind (or coplanar with) the plane.
    Back,
    /// All vertices are coplanar with the plane.
    Coplanar,
    /// Polygon straddles the plane.
    Spanning,
}

/// Classify a triangle relative to a plane using a floating-point epsilon.
///
/// Returns the classification and the individual vertex classifications.
///
/// **Note:** Uses `BSP_PLANE_TOLERANCE = 1e-5` for robustness against
/// near-coplanar vertices in the BSP pipeline.  For exact classification
/// without epsilon (Phase 8 Boolean ops), use [`classify_triangle_exact`].
pub fn classify_triangle(
    verts: [VertexId; 3],
    pool: &VertexPool,
    plane: &Plane,
) -> (PolygonClassification, [PointClassification; 3]) {
    let classifications = [
        plane.classify_point_with_eps(&pool.position(verts[0]), BSP_PLANE_TOLERANCE),
        plane.classify_point_with_eps(&pool.position(verts[1]), BSP_PLANE_TOLERANCE),
        plane.classify_point_with_eps(&pool.position(verts[2]), BSP_PLANE_TOLERANCE),
    ];

    let mut front = false;
    let mut back = false;
    for &c in &classifications {
        match c {
            PointClassification::Front => front = true,
            PointClassification::Back => back = true,
            PointClassification::Coplanar => {}
        }
    }

    let poly_class = match (front, back) {
        (true, true) => PolygonClassification::Spanning,
        (true, false) => PolygonClassification::Front,
        (false, true) => PolygonClassification::Back,
        (false, false) => PolygonClassification::Coplanar,
    };

    (poly_class, classifications)
}

// ── Exact classification (Phase 8) ───────────────────────────────────────────

/// Classify a triangle against a plane **without any floating-point epsilon**,
/// using Shewchuk's exact [`orient_3d`] predicate.
///
/// The plane is defined by three points `pa`, `pb`, `pc` in CCW order.
/// A vertex P is classified:
/// - `Front`   iff `orient_3d(pa, pb, pc, P) == Positive`
/// - `Back`    iff `orient_3d(pa, pb, pc, P) == Negative`
/// - `Coplanar` iff `orient_3d(pa, pb, pc, P) == Degenerate`
///
/// ## Theorem — No False Spanning Triangle
///
/// Because `orient_3d` is evaluated with adaptive-precision arithmetic
/// (Shewchuk 1997), the sign it returns is always *exactly* correct.
/// A triangle classified as non-spanning (all Front, all Back, or all
/// Coplanar) provably does not straddle the plane; no epsilon-proximity
/// band can misclassify a vertex on the "wrong" side.
///
/// The BSP-tolerance epsilon `1e-5` used by [`classify_triangle`] can
/// misclassify nearly-coplanar vertices in double-precision arithmetic;
/// `classify_triangle_exact` is immune to this.
pub fn classify_triangle_exact(
    verts: [VertexId; 3],
    pool: &VertexPool,
    pa: &Point3r,
    pb: &Point3r,
    pc: &Point3r,
) -> PolygonClassification {
    let arr = |p: &Point3r| [p.x, p.y, p.z];
    let pap = arr(pa); let pbp = arr(pb); let pcp = arr(pc);

    let orientations = [
        orient_3d(pap, pbp, pcp, arr(pool.position(verts[0]))),
        orient_3d(pap, pbp, pcp, arr(pool.position(verts[1]))),
        orient_3d(pap, pbp, pcp, arr(pool.position(verts[2]))),
    ];

    let any_pos = orientations.iter().any(|o| *o == Orientation::Positive);
    let any_neg = orientations.iter().any(|o| *o == Orientation::Negative);

    match (any_pos, any_neg) {
        (true, true)   => PolygonClassification::Spanning,
        (true, false)  => PolygonClassification::Front,
        (false, true)  => PolygonClassification::Back,
        (false, false) => PolygonClassification::Coplanar,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::scalar::Vector3r;
    use crate::storage::vertex_pool::VertexPool;

    fn make_pool_with_triangle(pts: [[f64; 3]; 3]) -> (VertexPool, [VertexId; 3]) {
        let mut pool = VertexPool::default_millifluidic();
        let zero = Vector3r::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(pts[0][0], pts[0][1], pts[0][2]), zero);
        let v1 = pool.insert_or_weld(Point3r::new(pts[1][0], pts[1][1], pts[1][2]), zero);
        let v2 = pool.insert_or_weld(Point3r::new(pts[2][0], pts[2][1], pts[2][2]), zero);
        (pool, [v0, v1, v2])
    }

    #[test]
    fn exact_all_front() {
        let (pool, verts) = make_pool_with_triangle([
            [1.0, 0.0, 1.0], [2.0, 0.0, 1.0], [1.0, 1.0, 1.0],
        ]);
        // Plane z=0, normal +z: pa=(0,0,0), pb=(1,0,0), pc=(0,1,0)
        let pa = Point3r::new(0.0, 0.0, 0.0);
        let pb = Point3r::new(1.0, 0.0, 0.0);
        let pc = Point3r::new(0.0, 1.0, 0.0);
        let class = classify_triangle_exact(verts, &pool, &pa, &pb, &pc);
        assert_eq!(class, PolygonClassification::Front);
    }

    #[test]
    fn exact_all_back() {
        let (pool, verts) = make_pool_with_triangle([
            [1.0, 0.0, -1.0], [2.0, 0.0, -1.0], [1.0, 1.0, -1.0],
        ]);
        let pa = Point3r::new(0.0, 0.0, 0.0);
        let pb = Point3r::new(1.0, 0.0, 0.0);
        let pc = Point3r::new(0.0, 1.0, 0.0);
        let class = classify_triangle_exact(verts, &pool, &pa, &pb, &pc);
        assert_eq!(class, PolygonClassification::Back);
    }

    #[test]
    fn exact_coplanar() {
        let (pool, verts) = make_pool_with_triangle([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        ]);
        let pa = Point3r::new(0.0, 0.0, 0.0);
        let pb = Point3r::new(1.0, 0.0, 0.0);
        let pc = Point3r::new(0.0, 1.0, 0.0);
        let class = classify_triangle_exact(verts, &pool, &pa, &pb, &pc);
        assert_eq!(class, PolygonClassification::Coplanar);
    }

    #[test]
    fn exact_spanning() {
        let (pool, verts) = make_pool_with_triangle([
            [1.0, 0.0, 1.0], [2.0, 0.0, -1.0], [1.0, 1.0, 0.0],
        ]);
        let pa = Point3r::new(0.0, 0.0, 0.0);
        let pb = Point3r::new(1.0, 0.0, 0.0);
        let pc = Point3r::new(0.0, 1.0, 0.0);
        let class = classify_triangle_exact(verts, &pool, &pa, &pb, &pc);
        assert_eq!(class, PolygonClassification::Spanning);
    }
}
