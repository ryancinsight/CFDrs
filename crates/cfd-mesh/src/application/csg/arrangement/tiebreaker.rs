//! Tiebreaker predicates for GWN boundary fragments in CSG classification.
//!
//! When GWN(q, M) lies in the band `[GWN_OUTSIDE_THRESHOLD, GWN_INSIDE_THRESHOLD]`
//! (nominally 0.35 ≤ |wn| ≤ 0.65), the fragment centroid lies on or very near the
//! opposing mesh surface.  Two tiebreakers resolve ambiguity in order:
//!
//! 1. **Coplanarity + normal comparison** (`coplanarity_tiebreak_*`):
//!    For each reference face whose plane exactly contains the centroid
//!    (`orient3d = 0`), vote on whether the fragment and reference normals are
//!    co-directed (exterior / `CoplanarSame`) or counter-directed (interior /
//!    `CoplanarOpposite`).  A majority vote is returned.
//!
//! 2. **Nearest-face signed distance** (`nearest_face_tiebreak_*`):
//!    Find the reference face whose centroid is closest (L² norm) to the
//!    fragment centroid and compute the signed distance of the fragment centroid
//!    from that face's plane.  Negative → behind face → Inside; positive →
//!    in front → Outside.
//!
//! Conservative fallback (when all tiebreakers are inconclusive): `CoplanarSame`.

use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::topology::predicates::{orient3d, Sign};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

use super::gwn::PreparedFace;

// ── Fragment classification result ────────────────────────────────────────────

/// CSG fragment classification relative to the opposing mesh.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FragmentClass {
    /// Centroid lies strictly inside the opposing mesh.
    Inside,
    /// Centroid lies strictly outside the opposing mesh.
    Outside,
    /// Centroid lies on a face plane with co-directed normal (shared boundary, exterior side).
    CoplanarSame,
    /// Centroid lies on a face plane with counter-directed normal (shared boundary, interior side).
    CoplanarOpposite,
}

// ── Tiebreaker 1: coplanarity ─────────────────────────────────────────────────

/// Coplanarity vote against pool-based reference faces.
///
/// Returns `Some(class)` when a majority vote exists; `None` when tied or
/// when no reference face is coplanar with `centroid`.
pub(super) fn coplanarity_tiebreak_pool(
    centroid: &Point3r,
    frag_normal: &Vector3r,
    faces: &[FaceData],
    pool: &VertexPool<f64>,
) -> Option<FragmentClass> {
    let mut interior = 0i32;
    let mut exterior = 0i32;
    for oface in faces {
        let pa = pool.position(oface.vertices[0]);
        let pb = pool.position(oface.vertices[1]);
        let pc = pool.position(oface.vertices[2]);
        if orient3d(pa, pb, pc, centroid) != Sign::Zero {
            continue;
        }
        let n_face = Vector3r::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z).cross(&Vector3r::new(
            pc.x - pa.x,
            pc.y - pa.y,
            pc.z - pa.z,
        ));
        let dot = n_face.dot(frag_normal);
        if dot > 0.0 {
            exterior += 1;
        } else if dot < 0.0 {
            interior += 1;
        }
    }
    majority_vote(interior, exterior)
}

/// Coplanarity vote against precomputed `PreparedFace` geometry.
pub(super) fn coplanarity_tiebreak_prepared(
    centroid: &Point3r,
    frag_normal: &Vector3r,
    faces: &[PreparedFace],
) -> Option<FragmentClass> {
    let mut interior = 0i32;
    let mut exterior = 0i32;
    for face in faces {
        if orient3d(&face.a, &face.b, &face.c, centroid) != Sign::Zero {
            continue;
        }
        let dot = face.normal.dot(frag_normal);
        if dot > 0.0 {
            exterior += 1;
        } else if dot < 0.0 {
            interior += 1;
        }
    }
    majority_vote(interior, exterior)
}

#[inline]
fn majority_vote(interior: i32, exterior: i32) -> Option<FragmentClass> {
    if interior > exterior {
        Some(FragmentClass::CoplanarOpposite)
    } else if exterior > interior {
        Some(FragmentClass::CoplanarSame)
    } else {
        None
    }
}

// ── Tiebreaker 2: nearest-face signed distance ────────────────────────────────

/// Nearest-face signed distance tiebreak against pool-based reference faces.
///
/// Returns `Inside` / `Outside` when the signed distance magnitude exceeds 1e-9;
/// falls back to `CoplanarSame` (conservative) when the centroid is on-boundary.
pub(super) fn nearest_face_tiebreak_pool(
    centroid: &Point3r,
    faces: &[FaceData],
    pool: &VertexPool<f64>,
) -> FragmentClass {
    let mut best_dist_sq = f64::MAX;
    let mut best_sign = 0.0_f64;
    for oface in faces {
        let pa = pool.position(oface.vertices[0]);
        let pb = pool.position(oface.vertices[1]);
        let pc = pool.position(oface.vertices[2]);
        let fc_x = (pa.x + pb.x + pc.x) / 3.0;
        let fc_y = (pa.y + pb.y + pc.y) / 3.0;
        let fc_z = (pa.z + pb.z + pc.z) / 3.0;
        let d = (centroid.x - fc_x, centroid.y - fc_y, centroid.z - fc_z);
        let dist_sq = d.0 * d.0 + d.1 * d.1 + d.2 * d.2;
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            let n = Vector3r::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z).cross(&Vector3r::new(
                pc.x - pa.x,
                pc.y - pa.y,
                pc.z - pa.z,
            ));
            let cp = Vector3r::new(centroid.x - pa.x, centroid.y - pa.y, centroid.z - pa.z);
            best_sign = cp.dot(&n);
        }
    }
    classify_by_sign(best_sign)
}

/// Nearest-face signed distance tiebreak against precomputed `PreparedFace` geometry.
pub(super) fn nearest_face_tiebreak_prepared(
    centroid: &Point3r,
    faces: &[PreparedFace],
) -> FragmentClass {
    let mut best_dist_sq = f64::MAX;
    let mut best_sign = 0.0_f64;
    for face in faces {
        let d = (
            centroid.x - face.centroid.x,
            centroid.y - face.centroid.y,
            centroid.z - face.centroid.z,
        );
        let dist_sq = d.0 * d.0 + d.1 * d.1 + d.2 * d.2;
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            let cp = Vector3r::new(
                centroid.x - face.a.x,
                centroid.y - face.a.y,
                centroid.z - face.a.z,
            );
            best_sign = cp.dot(&face.normal);
        }
    }
    classify_by_sign(best_sign)
}

#[inline]
fn classify_by_sign(sign: f64) -> FragmentClass {
    if sign.abs() > 1e-9 {
        if sign < 0.0 {
            FragmentClass::Inside
        } else {
            FragmentClass::Outside
        }
    } else {
        // Conservative fallback: keeps coplanar seam faces in the output.
        FragmentClass::CoplanarSame
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::infrastructure::storage::face_store::FaceData;
    use crate::infrastructure::storage::vertex_pool::VertexPool;

    /// Majority vote: more interior → CoplanarOpposite.
    #[test]
    fn majority_vote_interior_wins() {
        assert_eq!(majority_vote(2, 1), Some(FragmentClass::CoplanarOpposite));
    }

    /// Majority vote: more exterior → CoplanarSame.
    #[test]
    fn majority_vote_exterior_wins() {
        assert_eq!(majority_vote(1, 3), Some(FragmentClass::CoplanarSame));
    }

    /// Majority vote: tied → None.
    #[test]
    fn majority_vote_tie_gives_none() {
        assert_eq!(majority_vote(2, 2), None);
    }

    /// classify_by_sign: large negative → Inside.
    #[test]
    fn classify_by_sign_negative_is_inside() {
        assert_eq!(classify_by_sign(-0.01), FragmentClass::Inside);
    }

    /// classify_by_sign: large positive → Outside.
    #[test]
    fn classify_by_sign_positive_is_outside() {
        assert_eq!(classify_by_sign(0.01), FragmentClass::Outside);
    }

    /// classify_by_sign: near-zero → CoplanarSame (conservative fallback).
    #[test]
    fn classify_by_sign_near_zero_is_coplanar_same() {
        assert_eq!(classify_by_sign(1e-12), FragmentClass::CoplanarSame);
    }

    /// Coplanarity tiebreak: centroid on the +Z face of a unit square → exterior vote.
    #[test]
    fn coplanarity_pool_outward_normal_votes_exterior() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(-0.5, -0.5, 0.5), n);
        let v1 = pool.insert_or_weld(Point3r::new(0.5, -0.5, 0.5), n);
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 0.5, 0.5), n);
        let faces = vec![FaceData::untagged(v0, v1, v2)];

        // Centroid on the face plane z = 0.5, frag_normal = +Z (same as CCW face normal).
        let centroid = Point3r::new(0.0, 0.0, 0.5);
        let frag_normal = Vector3r::new(0.0, 0.0, 1.0);
        let result = coplanarity_tiebreak_pool(&centroid, &frag_normal, &faces, &pool);
        assert!(
            matches!(result, Some(FragmentClass::CoplanarSame)),
            "co-directed normals → exterior (CoplanarSame), got {result:?}"
        );
    }
}
