//! Fragment classification for CSG mesh arrangement.
//!
//! This module composes the GWN computation ([`super::gwn`]) and tiebreaker
//! predicates ([`super::tiebreaker`]) into the public classification API.
//!
//! ## `classify_fragment` Decision Flow
//!
//! ```text
//! classify_fragment(centroid, frag_normal, other_faces, pool)
//!         │
//!  ┌──────▼──────────────────────────────────────────────┐
//!  │  GWN(centroid, other_faces)                          │
//!  │  |wn| > GWN_INSIDE_THRESHOLD  ──────► Inside        │
//!  │  |wn| < GWN_OUTSIDE_THRESHOLD ──────► Outside       │
//!  └──────┬──────────────────────────────────────────────┘
//!         │ band: GWN_OUTSIDE_THRESHOLD ≤ |wn| ≤ GWN_INSIDE_THRESHOLD
//!  ┌──────▼──────────────────────────────────────────────┐
//!  │  coplanarity_tiebreak (orient3d majority vote)       │
//!  │  → CoplanarSame / CoplanarOpposite on majority       │
//!  └──────┬──────────────────────────────────────────────┘
//!         │ tied / no coplanar faces
//!  ┌──────▼──────────────────────────────────────────────┐
//!  │  nearest_face_tiebreak (signed distance fallback)    │
//!  │  negative sign → Inside; positive → Outside          │
//!  │  near-zero → CoplanarSame (conservative)             │
//!  └──────────────────────────────────────────────────────┘
//! ```
//!
//! ## Complexity
//!
//! | Function | Complexity | Notes |
//! |----------|------------|-------|
//! | `gwn()` | O(n) | n = face count of reference mesh |
//! | `classify_fragment()` | O(n) | single GWN + O(n) tiebreakers |
//!
//! ## References
//!
//! - Jacobson et al. (2013), *Robust Inside-Outside Segmentation using
//!   Generalized Winding Numbers*, ACM SIGGRAPH.

// ── Public re-exports (stable API for fragment classification and callers) ───

pub use super::gwn::{gwn, prepare_classification_faces, wnnc_score, PreparedFace};
pub use super::gwn_bvh::{gwn_bvh, prepare_bvh_mesh, PreparedBvhMesh};
pub use super::tiebreaker::FragmentClass;

use super::gwn::{gwn_bounded_prepared, gwn_prepared};
use super::tiebreaker::{
    coplanarity_tiebreak_pool, coplanarity_tiebreak_prepared, nearest_face_tiebreak_pool,
    nearest_face_tiebreak_prepared,
};
use crate::domain::core::constants::{GWN_INSIDE_THRESHOLD, GWN_OUTSIDE_THRESHOLD};
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

// ── Fragment record ───────────────────────────────────────────────────────────

/// One subdivision fragment of a parent face produced by CDT co-refinement.
pub struct FragRecord {
    /// The triangulated sub-face with pool-registered vertex IDs.
    pub face: FaceData,
    /// Index of the parent face in the originating face slice (A or B).
    pub parent_idx: usize,
    /// True if this fragment originated from `mesh_a`, false if from `mesh_b`.
    pub from_a: bool,
}

// ── Classification API ────────────────────────────────────────────────────────

/// Classify whether a fragment's centroid is inside the opposing mesh.
///
/// Calls GWN, then tiebreakers in order — see module-level decision flow.
#[must_use]
pub fn classify_fragment(
    centroid: &Point3r,
    frag_normal: &Vector3r,
    other_faces: &[FaceData],
    pool: &VertexPool<f64>,
) -> FragmentClass {
    let wn_abs = gwn::<f64>(centroid, other_faces, pool).abs();
    if wn_abs > GWN_INSIDE_THRESHOLD {
        return FragmentClass::Inside;
    }
    if wn_abs < GWN_OUTSIDE_THRESHOLD {
        return FragmentClass::Outside;
    }
    if let Some(cls) = coplanarity_tiebreak_pool(centroid, frag_normal, other_faces, pool) {
        return cls;
    }
    nearest_face_tiebreak_pool(centroid, other_faces, pool)
}

/// Classify whether a fragment's centroid is inside the opposing mesh using
/// precomputed reference-face geometry.
///
/// Drop-in equivalent to [`classify_fragment`] for hot loops.
///
/// When the standard GWN falls in the band, a bounded GWN refinement is
/// attempted first — this resolves many near-surface cases without needing
/// the tiebreaker chain.
#[must_use]
pub fn classify_fragment_prepared(
    centroid: &Point3r,
    frag_normal: &Vector3r,
    other_faces: &[PreparedFace],
) -> FragmentClass {
    let wn_abs = gwn_prepared(centroid, other_faces).abs();
    if wn_abs > GWN_INSIDE_THRESHOLD {
        return FragmentClass::Inside;
    }
    if wn_abs < GWN_OUTSIDE_THRESHOLD {
        return FragmentClass::Outside;
    }
    // Bounded GWN refinement: the per-triangle solid angle clip stabilises
    // the result for near-surface queries.  This often resolves the band case
    // without needing orient3d / signed-distance tiebreakers.
    let wn_bounded_abs = gwn_bounded_prepared(centroid, other_faces).abs();
    if wn_bounded_abs > GWN_INSIDE_THRESHOLD {
        return FragmentClass::Inside;
    }
    if wn_bounded_abs < GWN_OUTSIDE_THRESHOLD {
        return FragmentClass::Outside;
    }
    if let Some(cls) = coplanarity_tiebreak_prepared(centroid, frag_normal, other_faces) {
        return cls;
    }
    nearest_face_tiebreak_prepared(centroid, other_faces)
}

// ── Geometric utilities ───────────────────────────────────────────────────────

/// Triangle centroid.
#[inline]
#[must_use]
pub fn centroid(tri: &[Point3r; 3]) -> Point3r {
    Point3r::new(
        (tri[0].x + tri[1].x + tri[2].x) / 3.0,
        (tri[0].y + tri[1].y + tri[2].y) / 3.0,
        (tri[0].z + tri[1].z + tri[2].z) / 3.0,
    )
}

/// Geometric normal of a triangle (not normalised).
#[inline]
#[must_use]
pub fn tri_normal(tri: &[Point3r; 3]) -> Vector3r {
    let ab = Vector3r::new(
        tri[1].x - tri[0].x,
        tri[1].y - tri[0].y,
        tri[1].z - tri[0].z,
    );
    let ac = Vector3r::new(
        tri[2].x - tri[0].x,
        tri[2].y - tri[0].y,
        tri[2].z - tri[0].z,
    );
    ab.cross(&ac)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::infrastructure::storage::face_store::FaceData;
    use crate::infrastructure::storage::vertex_pool::VertexPool;

    fn unit_cube_mesh() -> (VertexPool, Vec<FaceData>) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let s = 0.5_f64;
        let mut v = |x, y, z| pool.insert_or_weld(Point3r::new(x, y, z), n);
        let c000 = v(-s, -s, -s);
        let c100 = v(s, -s, -s);
        let c010 = v(-s, s, -s);
        let c110 = v(s, s, -s);
        let c001 = v(-s, -s, s);
        let c101 = v(s, -s, s);
        let c011 = v(-s, s, s);
        let c111 = v(s, s, s);
        let f = FaceData::untagged;
        let faces = vec![
            f(c000, c010, c110),
            f(c000, c110, c100),
            f(c001, c101, c111),
            f(c001, c111, c011),
            f(c000, c001, c011),
            f(c000, c011, c010),
            f(c100, c110, c111),
            f(c100, c111, c101),
            f(c000, c100, c101),
            f(c000, c101, c001),
            f(c010, c011, c111),
            f(c010, c111, c110),
        ];
        (pool, faces)
    }

    #[test]
    fn classify_inside_fragment() {
        let (pool, faces) = unit_cube_mesh();
        let c = Point3r::new(0.0, 0.0, 0.0);
        let n = Vector3r::new(0.0, 0.0, 1.0);
        assert_eq!(
            classify_fragment(&c, &n, &faces, &pool),
            FragmentClass::Inside
        );
    }

    #[test]
    fn classify_outside_fragment() {
        let (pool, faces) = unit_cube_mesh();
        let c = Point3r::new(5.0, 0.0, 0.0);
        let n = Vector3r::new(1.0, 0.0, 0.0);
        assert_eq!(
            classify_fragment(&c, &n, &faces, &pool),
            FragmentClass::Outside
        );
    }

    #[test]
    fn classify_coplanar_same_fragment() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(-0.5, -0.5, 0.5), n);
        let v1 = pool.insert_or_weld(Point3r::new(0.5, -0.5, 0.5), n);
        let v2 = pool.insert_or_weld(Point3r::new(-0.5, 0.5, 0.5), n);
        let v3 = pool.insert_or_weld(Point3r::new(0.5, 0.5, 0.5), n);
        let faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v1, v3, v2),
        ];
        let c = Point3r::new(0.0, 0.0, 0.5);
        let frag_n = Vector3r::new(0.0, 0.0, 1.0);
        let cls = classify_fragment(&c, &frag_n, &faces, &pool);
        assert!(
            matches!(cls, FragmentClass::CoplanarSame | FragmentClass::Outside),
            "outward-coplanar fragment should be CoplanarSame or Outside, got {cls:?}"
        );
    }

    #[test]
    fn centroid_helper_correct() {
        let tri = [
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(3.0, 0.0, 0.0),
            Point3r::new(0.0, 3.0, 0.0),
        ];
        let c = centroid(&tri);
        assert!((c.x - 1.0).abs() < 1e-10);
        assert!((c.y - 1.0).abs() < 1e-10);
        assert!(c.z.abs() < 1e-10);
    }

    #[test]
    fn tri_normal_ccw_xy_is_positive_z() {
        let tri = [
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(1.0, 0.0, 0.0),
            Point3r::new(0.0, 1.0, 0.0),
        ];
        let n = tri_normal(&tri);
        assert!(n.z > 0.0, "CCW XY triangle normal should be +Z, got {n:?}");
    }

    /// Prepared and non-prepared classification paths must produce identical results.
    #[test]
    fn classify_prepared_matches_unprepared() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        let samples = [
            (Point3r::new(0.0, 0.0, 0.0), Vector3r::new(0.0, 0.0, 1.0)),
            (Point3r::new(4.0, 0.0, 0.0), Vector3r::new(1.0, 0.0, 0.0)),
            (Point3r::new(0.0, 0.0, 0.5), Vector3r::new(0.0, 0.0, 1.0)),
            (Point3r::new(0.5, 0.0, 0.0), Vector3r::new(1.0, 0.0, 0.0)),
        ];
        for (c, n) in samples {
            let a = classify_fragment(&c, &n, &faces, &pool);
            let b = classify_fragment_prepared(&c, &n, &prepared);
            assert_eq!(a, b, "prepared/unprepared diverged at centroid={c:?}");
        }
    }
}
