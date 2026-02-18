//! BVH broad-phase for CSG Boolean operations.
//!
//! Finds candidate triangle pairs from two meshes whose axis-aligned bounding
//! boxes overlap.  Only overlapping pairs need exact narrow-phase testing.
//!
//! ## Complexity
//!
//! This O(n·m) implementation is sufficient for millifluidic meshes
//! (typically < 10 000 triangles).  A BVH-accelerated variant using the `bvh`
//! crate can be substituted here without changing the public API.
//!
//! ## Diagram
//!
//! ```text
//! mesh_A triangles:  [T0, T1, T2, … Tn]
//!                           │
//!                    AABB(Ti) overlaps?
//!                           │
//!                           ▼
//! mesh_B triangles:  [T0, T1, T2, … Tm]
//!                           │
//!                   → candidate pair (i, j)
//! ```

use crate::core::scalar::Real;
use crate::geometry::aabb::Aabb;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;

// ── Public types ──────────────────────────────────────────────────────────────

/// A candidate triangle pair `(face_a, face_b)` whose AABBs overlap.
///
/// Both indices are zero-based positions within their respective face slices.
/// A candidate pair is not guaranteed to geometrically intersect; it must be
/// passed to the exact narrow-phase test in [`crate::csg::intersect`].
///
/// # Invariant
/// `face_a < faces_a.len()` and `face_b < faces_b.len()` for the slices passed
/// to [`broad_phase_pairs`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CandidatePair {
    /// Zero-based index of the face in mesh A's face slice.
    pub face_a: usize,
    /// Zero-based index of the face in mesh B's face slice.
    pub face_b: usize,
}

// ── AABB ──────────────────────────────────────────────────────────────────────

/// AABB growth applied to each triangle to absorb floating-point rounding.
///
/// 1 nm — smaller than any millifluidic feature but large enough to avoid
/// dropping truly-intersecting pairs due to representational error.
const AABB_EPSILON: Real = 1e-9;

/// Compute the axis-aligned bounding box of a triangle.
///
/// The box is grown by [`AABB_EPSILON`] on each side to guard against
/// floating-point rounding when comparing AABBs of nearly-touching triangles.
pub fn triangle_aabb(face: &FaceData, pool: &VertexPool) -> Aabb {
    let a = *pool.position(face.vertices[0]);
    let b = *pool.position(face.vertices[1]);
    let c = *pool.position(face.vertices[2]);

    let mut aabb = Aabb::empty();
    aabb.expand(&a);
    aabb.expand(&b);
    aabb.expand(&c);

    aabb.min.x -= AABB_EPSILON;
    aabb.min.y -= AABB_EPSILON;
    aabb.min.z -= AABB_EPSILON;
    aabb.max.x += AABB_EPSILON;
    aabb.max.y += AABB_EPSILON;
    aabb.max.z += AABB_EPSILON;

    aabb
}

// ── Broad phase ───────────────────────────────────────────────────────────────

/// Find all pairs `(i, j)` where the AABB of `faces_a[i]` overlaps the
/// AABB of `faces_b[j]`.
///
/// Pre-computes all B-side AABBs to avoid recomputation for every A-face.
/// Returned pairs are in `(face_a, face_b)` order and must be validated by
/// the exact narrow-phase test before being used for mesh splitting.
///
/// # Complexity
/// `O(n·m)` where n = `faces_a.len()`, m = `faces_b.len()`.
///
/// # Theorem: Completeness
/// Every pair of triangles that geometrically intersect must have overlapping
/// AABBs (since the AABB contains the triangle).  Therefore no intersecting
/// pair is missing from the output.  ∎
pub fn broad_phase_pairs(
    faces_a: &[FaceData],
    pool_a: &VertexPool,
    faces_b: &[FaceData],
    pool_b: &VertexPool,
) -> Vec<CandidatePair> {
    // Pre-compute all B-side AABBs once.
    let aabbs_b: Vec<Aabb> = faces_b
        .iter()
        .map(|f| triangle_aabb(f, pool_b))
        .collect();

    let mut pairs = Vec::new();
    for (i, fa) in faces_a.iter().enumerate() {
        let aabb_a = triangle_aabb(fa, pool_a);
        for (j, aabb_b) in aabbs_b.iter().enumerate() {
            if aabb_a.intersects(aabb_b) {
                pairs.push(CandidatePair { face_a: i, face_b: j });
            }
        }
    }
    pairs
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::scalar::Point3r;

    fn make_pool_and_face(pts: [[f64; 3]; 3]) -> (VertexPool, FaceData) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(pts[0][0], pts[0][1], pts[0][2]), n);
        let v1 = pool.insert_or_weld(Point3r::new(pts[1][0], pts[1][1], pts[1][2]), n);
        let v2 = pool.insert_or_weld(Point3r::new(pts[2][0], pts[2][1], pts[2][2]), n);
        (pool, FaceData::untagged(v0, v1, v2))
    }

    #[test]
    fn overlapping_aabbs_yields_pair() {
        let (pa, fa) =
            make_pool_and_face([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let (pb, fb) =
            make_pool_and_face([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [0.5, 1.0, 0.0]]);
        let pairs = broad_phase_pairs(&[fa], &pa, &[fb], &pb);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], CandidatePair { face_a: 0, face_b: 0 });
    }

    #[test]
    fn non_overlapping_aabbs_no_pair() {
        let (pa, fa) =
            make_pool_and_face([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let (pb, fb) =
            make_pool_and_face([[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]]);
        let pairs = broad_phase_pairs(&[fa], &pa, &[fb], &pb);
        assert!(pairs.is_empty());
    }

    #[test]
    fn multiple_b_triangles() {
        let (pa, fa) =
            make_pool_and_face([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]);
        // Two B triangles: one overlapping, one not
        let (pb1, fb1) =
            make_pool_and_face([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [0.5, 1.0, 0.0]]);
        let (pb2, fb2) =
            make_pool_and_face([[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]]);

        // Build a synthetic combined pool and face list for B
        let mut combined_b = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();

        let mut map_b1 = [crate::core::index::VertexId::default(); 3];
        for (k, &vid) in fb1.vertices.iter().enumerate() {
            map_b1[k] = combined_b.insert_or_weld(*pb1.position(vid), n);
        }
        let mut map_b2 = [crate::core::index::VertexId::default(); 3];
        for (k, &vid) in fb2.vertices.iter().enumerate() {
            map_b2[k] = combined_b.insert_or_weld(*pb2.position(vid), n);
        }
        let faces_b = vec![
            FaceData::untagged(map_b1[0], map_b1[1], map_b1[2]),
            FaceData::untagged(map_b2[0], map_b2[1], map_b2[2]),
        ];

        let pairs = broad_phase_pairs(&[fa], &pa, &faces_b, &combined_b);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].face_a, 0);
        assert_eq!(pairs[0].face_b, 0); // only the first B triangle overlaps
    }

    #[test]
    fn triangle_aabb_bounds_are_correct() {
        let (p, f) =
            make_pool_and_face([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]]);
        let aabb = triangle_aabb(&f, &p);
        assert!(aabb.min.x <= 0.0);
        assert!(aabb.min.y <= 0.0);
        assert!(aabb.max.x >= 2.0);
        assert!(aabb.max.y >= 3.0);
        // Epsilon growth
        assert!(aabb.min.x < 0.0);
        assert!(aabb.max.x > 2.0);
    }
}
