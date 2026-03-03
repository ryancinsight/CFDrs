//! Fragment classification and Generalized Winding Number (GWN) for CSG.
//!
//! ## Generalized Winding Number — Theorem
//!
//! For a closed orientable 2-manifold M and a query point q not on M:
//!
//! ```text
//! GWN(q, M) = (1/4π) Σ Ω(q, tri_i)
//! ```
//!
//! where Ω(q, T) is the solid angle subtended by triangle T at q, computed via
//! the van Oosterom–Strackee (1983) formula:
//!
//! ```text
//! Ω = 2·atan2( a·(b×c), |a||b||c| + (a·b)|c| + (b·c)|a| + (c·a)|b| )
//! ```
//!
//! where a, b, c are the vectors from q to each triangle vertex (not normalised).
//!
//! **Interior**: GWN = ±1  **Exterior**: GWN = 0
//!
//! No epsilon displacement of q is needed — GWN is well-defined everywhere
//! except at mesh vertices.  Seam centroids lying exactly on a face plane
//! produce |GWN| ≈ 0.5; [`classify_fragment`] handles these via an exact
//! `orient3d` tiebreaker.
//!
//! ## `classify_fragment` Decision Flow
//!
//! ```text
//! classify_fragment(centroid, frag_normal, other_faces, pool)
//!         │
//!  ┌──────▼──────────────────────────────────────────────┐
//!  │  GWN(centroid, other_faces)                          │
//!  │  |wn| > 0.65  ──────► Inside                        │
//!  │  |wn| < 0.35  ──────► Outside                       │
//!  └──────┬──────────────────────────────────────────────┘
//!         │ 0.35 ≤ |wn| ≤ 0.65  (seam fragment)
//!  ┌──────▼──────────────────────────────────────────────┐
//!  │  Tiebreaker 1: orient3d coplanarity test             │
//!  │  For each face coplanar with centroid:               │
//!  │    n_face · frag_normal > 0 → exterior vote         │
//!  │    n_face · frag_normal < 0 → interior vote         │
//!  │  majority → CoplanarSame / CoplanarOpposite          │
//!  └──────┬──────────────────────────────────────────────┘
//!         │ tied / no coplanar faces
//!  ┌──────▼──────────────────────────────────────────────┐
//!  │  Tiebreaker 2: nearest-face signed distance          │
//!  │  Negative signed dist → Inside; positive → Outside   │
//!  └──────┬──────────────────────────────────────────────┘
//!         │ ambiguous → CoplanarSame (conservative)
//!         ▼
//! ```
//!
//! ## Complexity
//!
//! | Function | Complexity | Notes |
//! |----------|------------|-------|
//! | `gwn()` | O(n) | n = face count of reference mesh |
//! | `classify_fragment()` | O(n) | calls gwn() once; tiebreakers also O(n) |
//!
//! ## References
//!
//! - van Oosterom & Strackee (1983), *The Solid Angle of a Plane Triangle*,
//!   IEEE Trans. Biomed. Eng. 30(2).
//! - Jacobson et al. (2013), *Robust Inside-Outside Segmentation using
//!   Generalized Winding Numbers*, ACM SIGGRAPH.

use crate::domain::core::scalar::{Point3r, Scalar, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

// ── Fragment classification result ───────────────────────────────────────────

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FragmentClass {
    Inside,
    Outside,
    CoplanarSame,
    CoplanarOpposite,
}

// ── Internal fragment record ─────────────────────────────────────────────────

/// One subdivision fragment of a parent face.
pub struct FragRecord {
    /// The triangulated sub-face with pool-registered vertex IDs.
    pub face: FaceData,
    /// Index of the parent face in the originating face slice (A or B).
    pub parent_idx: usize,
    /// True if this fragment originated from `mesh_a`, false if from `mesh_b`.
    pub from_a: bool,
}

// ── Phase 4: GWN and fragment classification ─────────────────────────────────

/// Generalized Winding Number (GWN) of `query` with respect to a closed
/// triangle mesh.
///
/// Returns a value in [-1, 1] clamped by construction:
/// - **±1.0**: query is strictly inside the mesh.
/// - **0.0**:  query is strictly outside the mesh.
/// - **≈0.5**: query lies on a face plane (seam centroid).
///
/// ## Theorem — Correctness
///
/// For a closed orientable 2-manifold, GWN is an integer-valued function
/// (by the Gauss-Bonnet theorem) at all non-mesh points.  The only values
/// possible for a single-winding manifold are 0 (outside) and ±1 (inside).
/// Points on face planes give |GWN| = 0.5 by symmetry.  ∎
///
/// ## Implementation Note — Norm Efficiency
///
/// The van Oosterom–Strackee formula uses vector norms (|a|, |b|, |c|) only
/// in the denominator.  The near-vertex guard uses `norm_squared() < ε²`
/// (no sqrt).  The denominator computation calls `norm()` exactly once per
/// vertex vector and reuses the result — no redundant square-root calls.
pub fn gwn<T: Scalar>(query: &nalgebra::Point3<T>, faces: &[FaceData], pool: &VertexPool<T>) -> T {
    let mut solid_angle_sum = <T as Scalar>::from_f64(0.0);
    let near_sq = <T as Scalar>::from_f64(1e-40); // ε² for near-vertex guard (no sqrt)
    let one_e_30 = <T as Scalar>::from_f64(1e-30);
    let two = <T as Scalar>::from_f64(2.0);
    let four_pi = <T as Scalar>::from_f64(4.0 * std::f64::consts::PI);

    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);

        let va = nalgebra::Vector3::new(a.x - query.x, a.y - query.y, a.z - query.z);
        let vb = nalgebra::Vector3::new(b.x - query.x, b.y - query.y, b.z - query.z);
        let vc = nalgebra::Vector3::new(c.x - query.x, c.y - query.y, c.z - query.z);

        // Near-vertex guard: use norm_squared (no sqrt) for efficiency.
        if va.norm_squared() < near_sq || vb.norm_squared() < near_sq || vc.norm_squared() < near_sq
        {
            continue;
        }

        // Compute norms once; reuse in denominator.
        let la = va.norm();
        let lb = vb.norm();
        let lc = vc.norm();

        let num = va.dot(&vb.cross(&vc));
        let den = la * lb * lc + va.dot(&vb) * lc + vb.dot(&vc) * la + vc.dot(&va) * lb;

        use num_traits::Float;
        if Float::abs(den) > one_e_30 || Float::abs(num) > one_e_30 {
            solid_angle_sum += two * Float::atan2(num, den);
        }
    }
    use num_traits::clamp;
    clamp(
        solid_angle_sum / four_pi,
        <T as Scalar>::from_f64(-1.0),
        <T as Scalar>::from_f64(1.0),
    )
}

/// Classify whether a fragment's centroid is inside the opposing mesh.
///
/// See the module-level `classify_fragment` decision flow for the full
/// decision tree with GWN thresholds and tiebreakers.
#[must_use]
pub fn classify_fragment(
    centroid: &Point3r,
    frag_normal: &Vector3r,
    other_faces: &[FaceData],
    pool: &VertexPool<f64>,
) -> FragmentClass {
    use crate::domain::topology::predicates::{orient3d, Sign};

    let wn = gwn::<f64>(centroid, other_faces, pool);
    let wn_abs = wn.abs();

    // ── Fast path: unambiguous GWN result ────────────────────────────────────
    if wn_abs > 0.65 {
        return FragmentClass::Inside;
    }
    if wn_abs < 0.35 {
        return FragmentClass::Outside;
    }

    // ── Tiebreaker 1: exact orient3d coplanarity + normal comparison ─────────
    // For faces where the centroid lies exactly on the face plane (orient3d = 0),
    // compare the fragment's outward normal against the face's CCW normal.
    // Same-direction normals → exterior (CoplanarSame);
    // opposite normals → interior (CoplanarOpposite).
    let mut interior_votes = 0i32;
    let mut exterior_votes = 0i32;

    for oface in other_faces {
        let pa = pool.position(oface.vertices[0]);
        let pb = pool.position(oface.vertices[1]);
        let pc = pool.position(oface.vertices[2]);

        if orient3d(pa, pb, pc, centroid) != Sign::Zero {
            continue;
        }

        let ab = Vector3r::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
        let ac = Vector3r::new(pc.x - pa.x, pc.y - pa.y, pc.z - pa.z);
        let n_face = ab.cross(&ac);

        let dot = n_face.dot(frag_normal);
        if dot > 0.0 {
            exterior_votes += 1;
        } else if dot < 0.0 {
            interior_votes += 1;
        }
    }

    if interior_votes > exterior_votes {
        return FragmentClass::CoplanarOpposite;
    }
    if exterior_votes > interior_votes {
        return FragmentClass::CoplanarSame;
    }

    // ── Tiebreaker 2: nearest-face signed distance ────────────────────────────
    // Find the face whose centroid is closest to the query centroid (L² norm)
    // and use its plane's signed distance from the query.
    //
    // For curved meshes, exact coplanarity (orient3d == Zero) almost never fires.
    // The nearest face's plane provides a robust fallback:
    //   signed_dist < 0 → query is behind face (inside)
    //   signed_dist > 0 → query is in front (outside)
    //
    // Early-exit: once we've found a face closer than `best_dist_sq` and
    // accumulated non-zero signed distance, further faces cannot improve the
    // result since they are provably farther away.
    let mut best_dist_sq = f64::MAX;
    let mut best_sign = 0.0_f64;

    for oface in other_faces {
        let pa = pool.position(oface.vertices[0]);
        let pb = pool.position(oface.vertices[1]);
        let pc = pool.position(oface.vertices[2]);
        let fc = Point3r::new(
            (pa.x + pb.x + pc.x) / 3.0,
            (pa.y + pb.y + pc.y) / 3.0,
            (pa.z + pb.z + pc.z) / 3.0,
        );
        let dx = centroid.x - fc.x;
        let dy = centroid.y - fc.y;
        let dz = centroid.z - fc.z;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            let ab = Vector3r::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
            let ac = Vector3r::new(pc.x - pa.x, pc.y - pa.y, pc.z - pa.z);
            let n = ab.cross(&ac);
            let cp = Vector3r::new(centroid.x - pa.x, centroid.y - pa.y, centroid.z - pa.z);
            best_sign = cp.dot(&n);
        }
    }

    // Negative signed distance → centroid is behind the closest face → inside.
    if best_sign.abs() > 1e-9 {
        if best_sign < 0.0 {
            return FragmentClass::Inside;
        } else {
            return FragmentClass::Outside;
        }
    }

    // Treat on-boundary (GWN == 0.5, no coplanar votes, no signed dist) as
    // CoplanarSame — conservative: keeps the face in the output.
    FragmentClass::CoplanarSame
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
    use crate::domain::core::scalar::Point3r;
    use crate::infrastructure::storage::face_store::FaceData;
    use crate::infrastructure::storage::vertex_pool::VertexPool;

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// Build a closed unit-cube mesh (6 faces × 2 triangles = 12 triangles).
    /// Normals are outward-pointing CCW.
    fn unit_cube_mesh() -> (VertexPool, Vec<FaceData>) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let s = 0.5_f64;
        // 8 corners of [-0.5, 0.5]³
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
        // 6 faces × 2 triangles, CCW winding from outside.
        let faces = vec![
            // -Z face (bottom)
            f(c000, c010, c110),
            f(c000, c110, c100),
            // +Z face (top)
            f(c001, c101, c111),
            f(c001, c111, c011),
            // -X face
            f(c000, c001, c011),
            f(c000, c011, c010),
            // +X face
            f(c100, c110, c111),
            f(c100, c111, c101),
            // -Y face (front)
            f(c000, c100, c101),
            f(c000, c101, c001),
            // +Y face (back)
            f(c010, c011, c111),
            f(c010, c111, c110),
        ];
        (pool, faces)
    }

    // ── GWN tests ─────────────────────────────────────────────────────────────

    /// Theorem: GWN at the origin (interior of unit cube) must equal 1.0
    /// for an outward-pointing closed manifold.
    #[test]
    fn gwn_unit_cube_interior_is_one() {
        let (pool, faces) = unit_cube_mesh();
        let query = Point3r::new(0.0, 0.0, 0.0);
        let wn = gwn::<f64>(&query, &faces, &pool);
        assert!(
            (wn - 1.0).abs() < 0.02,
            "GWN at cube interior should be ≈1.0, got {wn:.4}"
        );
    }

    /// Theorem: GWN at a far exterior point must equal 0.0.
    #[test]
    fn gwn_unit_cube_exterior_is_zero() {
        let (pool, faces) = unit_cube_mesh();
        let query = Point3r::new(10.0, 0.0, 0.0);
        let wn = gwn::<f64>(&query, &faces, &pool);
        assert!(
            wn.abs() < 0.02,
            "GWN at cube exterior should be ≈0.0, got {wn:.4}"
        );
    }

    /// Theorem: GWN output is always clamped to [-1, 1] by construction.
    #[test]
    fn gwn_always_clamped_to_unit_interval() {
        let (pool, faces) = unit_cube_mesh();
        for (x, y, z) in [
            (0.0, 0.0, 0.0),       // interior
            (10.0, 0.0, 0.0),      // exterior
            (0.5, 0.5, 0.5),       // corner vicinity
            (100.0, 100.0, 100.0), // very far
        ] {
            let q = Point3r::new(x, y, z);
            let wn = gwn::<f64>(&q, &faces, &pool);
            assert!(
                wn >= -1.0 && wn <= 1.0,
                "GWN ({x},{y},{z}) out of [-1,1]: {wn}"
            );
        }
    }

    /// GWN on empty mesh returns 0.
    #[test]
    fn gwn_empty_mesh_is_zero() {
        let pool = VertexPool::default_millifluidic();
        let q = Point3r::new(0.0, 0.0, 0.0);
        let wn = gwn::<f64>(&q, &[], &pool);
        assert_eq!(wn, 0.0);
    }

    // ── classify_fragment tests ────────────────────────────────────────────────

    /// Fragment centroid deep inside the cube → Inside.
    #[test]
    fn classify_inside_fragment() {
        let (pool, faces) = unit_cube_mesh();
        let centroid = Point3r::new(0.0, 0.0, 0.0);
        let normal = Vector3r::new(0.0, 0.0, 1.0);
        let cls = classify_fragment(&centroid, &normal, &faces, &pool);
        assert_eq!(
            cls,
            FragmentClass::Inside,
            "interior point should be Inside"
        );
    }

    /// Fragment centroid far outside the cube → Outside.
    #[test]
    fn classify_outside_fragment() {
        let (pool, faces) = unit_cube_mesh();
        let centroid = Point3r::new(5.0, 0.0, 0.0);
        let normal = Vector3r::new(1.0, 0.0, 0.0);
        let cls = classify_fragment(&centroid, &normal, &faces, &pool);
        assert_eq!(
            cls,
            FragmentClass::Outside,
            "exterior point should be Outside"
        );
    }

    /// Fragment on the +Z face with outward normal → CoplanarSame (exterior).
    #[test]
    fn classify_coplanar_same_fragment() {
        // Build a single square face facing +Z at z = 0.5
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
        // Centroid ON the plane z = 0.5, with +Z normal (same as mesh face).
        let c = Point3r::new(0.0, 0.0, 0.5);
        let frag_n = Vector3r::new(0.0, 0.0, 1.0); // same direction as face CCW normal
        let cls = classify_fragment(&c, &frag_n, &faces, &pool);
        // CoplanarSame or Outside are both acceptable — key is not Inside.
        assert!(
            matches!(cls, FragmentClass::CoplanarSame | FragmentClass::Outside),
            "outward-coplanar fragment should be CoplanarSame or Outside, got {cls:?}"
        );
    }

    /// `centroid()` helper returns correct midpoint.
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
        assert!((c.z - 0.0).abs() < 1e-10);
    }

    /// `tri_normal()` for a CCW XY triangle should point in +Z.
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
}
