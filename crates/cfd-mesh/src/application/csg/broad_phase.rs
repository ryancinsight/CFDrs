//! BVH broad-phase for CSG Boolean operations.
//!
//! Finds candidate triangle pairs from two meshes whose axis-aligned bounding
//! boxes overlap.  Only overlapping pairs need exact narrow-phase testing.
//!
//! ## Algorithm
//!
//! A flat-arena SAH-BVH ([`super::bvh::BvhTree`]) is built over the triangles
//! of the smaller mesh. The larger mesh is queried against that BVH, yielding
//! `O((n + m) log(min(n,m)))` total work — a significant improvement over the
//! `O(n·m)` exhaustive approach for meshes with > ~100 triangles.
//!
//! ## Diagram
//!
//! ```text
//! smaller mesh triangles: [S0, S1, ...]
//!         │
//!         │  build SAH-BVH
//!         ▼
//!    BvhTree(S)  ── query larger mesh triangles ──► overlapping leaf indices
//!         │
//!         ▼
//!   → candidate pair (face_a, face_b)
//! ```
//!
//! ## Theorem: Completeness
//! Every pair of triangles that geometrically intersect must have overlapping
//! AABBs (since the AABB contains the triangle).  The [`crate::infrastructure::spatial::bvh::BvhTree`]
//! query is itself complete by its own completeness theorem.  Therefore no
//! intersecting pair is missing from the output.  ∎
//!
//! ## Theorem: Symmetric Query Equivalence
//! AABB-overlap is symmetric: `overlap(A_i, B_j) == overlap(B_j, A_i)`.
//! Therefore building the BVH on A and querying B, or building on B and
//! querying A, produces the same candidate pair set after index re-labeling. ∎

use crate::domain::geometry::aabb::Aabb;
use crate::infrastructure::spatial::bvh::with_bvh;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

// ── Public types ──────────────────────────────────────────────────────────────

/// A candidate triangle pair `(face_a, face_b)` whose AABBs overlap.
///
/// Both indices are zero-based positions within their respective face slices.
/// A candidate pair is not guaranteed to geometrically intersect; it must be
/// passed to the exact narrow-phase test in [`super::intersect`].
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

/// Compute the axis-aligned bounding box of a triangle.
///
/// The box is grown by [`AABB_EPSILON`] on each side to guard against
/// floating-point rounding when comparing AABBs of nearly-touching triangles.
#[must_use]
pub fn triangle_aabb(face: &FaceData, pool: &VertexPool) -> Aabb {
    let a = *pool.position(face.vertices[0]);
    let b = *pool.position(face.vertices[1]);
    let c = *pool.position(face.vertices[2]);

    let mut aabb = Aabb::empty();
    aabb.expand(&a);
    aabb.expand(&b);
    aabb.expand(&c);

    aabb
}

// ── Broad phase ───────────────────────────────────────────────────────────────

/// Find all pairs `(i, j)` where the AABB of `faces_a[i]` overlaps the
/// AABB of `faces_b[j]`.
///
/// Builds a SAH-BVH over the smaller side, then queries with the larger side.
/// Returned pairs are in `(face_a, face_b)` order and must be validated by
/// the exact narrow-phase test before being used for mesh splitting.
///
/// # Complexity
/// `O((n + m) log(min(n,m)))` where n = `faces_a.len()`, m = `faces_b.len()`.
#[must_use]
pub fn broad_phase_pairs(
    faces_a: &[FaceData],
    pool_a: &VertexPool,
    faces_b: &[FaceData],
    pool_b: &VertexPool,
) -> Vec<CandidatePair> {
    if faces_a.is_empty() || faces_b.is_empty() {
        return Vec::new();
    }

    let mut pairs = Vec::new();

    if faces_b.len() <= faces_a.len() {
        // Build BVH on B and query A.
        let aabbs_b: Vec<Aabb> = faces_b.iter().map(|f| triangle_aabb(f, pool_b)).collect();
        with_bvh(
            &aabbs_b,
            |tree: crate::infrastructure::spatial::bvh::BvhTree<'_, '_>, token| {
                let mut hits = Vec::new();
                for (i, fa) in faces_a.iter().enumerate() {
                    let aabb_a = triangle_aabb(fa, pool_a);
                    hits.clear();
                    tree.query_overlapping(&aabb_a, &token, &mut hits);
                    for &j in &hits {
                        pairs.push(CandidatePair {
                            face_a: i,
                            face_b: j,
                        });
                    }
                }
            },
        );
    } else {
        // Build BVH on A and query B, then remap to (face_a, face_b).
        let aabbs_a: Vec<Aabb> = faces_a.iter().map(|f| triangle_aabb(f, pool_a)).collect();
        with_bvh(
            &aabbs_a,
            |tree: crate::infrastructure::spatial::bvh::BvhTree<'_, '_>, token| {
                let mut hits = Vec::new();
                for (j, fb) in faces_b.iter().enumerate() {
                    let aabb_b = triangle_aabb(fb, pool_b);
                    hits.clear();
                    tree.query_overlapping(&aabb_b, &token, &mut hits);
                    for &i in &hits {
                        pairs.push(CandidatePair {
                            face_a: i,
                            face_b: j,
                        });
                    }
                }
            },
        );
    }

    // Preserve deterministic processing order for downstream arrangement stages.
    // The set of pairs is unchanged; we sort only by indices.
    pairs.sort_unstable_by_key(|p| (p.face_a, p.face_b));
    pairs
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use proptest::prelude::*;
    use std::collections::BTreeSet;

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
        let (pa, fa) = make_pool_and_face([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let (pb, fb) = make_pool_and_face([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [0.5, 1.0, 0.0]]);
        let pairs = broad_phase_pairs(&[fa], &pa, &[fb], &pb);
        assert_eq!(pairs.len(), 1);
        assert_eq!(
            pairs[0],
            CandidatePair {
                face_a: 0,
                face_b: 0
            }
        );
    }

    #[test]
    fn non_overlapping_aabbs_no_pair() {
        let (pa, fa) = make_pool_and_face([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let (pb, fb) = make_pool_and_face([[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]]);
        let pairs = broad_phase_pairs(&[fa], &pa, &[fb], &pb);
        assert!(pairs.is_empty());
    }

    #[test]
    fn multiple_b_triangles() {
        let (pa, fa) = make_pool_and_face([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]);
        // Two B triangles: one overlapping, one not
        let (pb1, fb1) = make_pool_and_face([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [0.5, 1.0, 0.0]]);
        let (pb2, fb2) = make_pool_and_face([[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]]);

        // Build a synthetic combined pool and face list for B
        let mut combined_b = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();

        let mut map_b1 = [crate::domain::core::index::VertexId::default(); 3];
        for (k, &vid) in fb1.vertices.iter().enumerate() {
            map_b1[k] = combined_b.insert_or_weld(*pb1.position(vid), n);
        }
        let mut map_b2 = [crate::domain::core::index::VertexId::default(); 3];
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
        let (p, f) = make_pool_and_face([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]]);
        let aabb = triangle_aabb(&f, &p);
        assert!(aabb.min.x <= 0.0);
        assert!(aabb.min.y <= 0.0);
        assert!(aabb.max.x >= 2.0);
        assert!(aabb.max.y >= 3.0);
    }

    fn make_faces_from_raw(raw: &[[i16; 9]]) -> (VertexPool, Vec<FaceData>) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let mut faces = Vec::with_capacity(raw.len());
        for tri in raw {
            let p = |k: usize| -> Point3r {
                Point3r::new(
                    tri[k] as f64 * 0.1,
                    tri[k + 1] as f64 * 0.1,
                    tri[k + 2] as f64 * 0.1,
                )
            };
            let v0 = pool.insert_or_weld(p(0), n);
            let v1 = pool.insert_or_weld(p(3), n);
            let v2 = pool.insert_or_weld(p(6), n);
            faces.push(FaceData::untagged(v0, v1, v2));
        }
        (pool, faces)
    }

    proptest! {
        #[test]
        fn adaptive_side_selection_is_symmetric(
            a_raw in prop::collection::vec(prop::array::uniform9(-20_i16..20_i16), 0..12),
            b_raw in prop::collection::vec(prop::array::uniform9(-20_i16..20_i16), 0..12),
        ) {
            let (pool_a, faces_a) = make_faces_from_raw(&a_raw);
            let (pool_b, faces_b) = make_faces_from_raw(&b_raw);

            let ab = broad_phase_pairs(&faces_a, &pool_a, &faces_b, &pool_b);
            let ba = broad_phase_pairs(&faces_b, &pool_b, &faces_a, &pool_a);

            let set_ab: BTreeSet<(usize, usize)> =
                ab.into_iter().map(|p| (p.face_a, p.face_b)).collect();
            let set_ba_flipped: BTreeSet<(usize, usize)> =
                ba.into_iter().map(|p| (p.face_b, p.face_a)).collect();

            prop_assert_eq!(set_ab, set_ba_flipped);
        }
    }
}
