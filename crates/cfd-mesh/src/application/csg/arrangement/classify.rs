use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
// Ã¢â€â‚¬Ã¢â€â‚¬ Internal types Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬

/// One subdivision fragment of a parent face.
pub struct FragRecord {
    /// The triangulated sub-face with pool-registered vertex IDs.
    pub face: FaceData,
    /// Index of the parent face in the originating face slice (A or B).
    pub parent_idx: usize,
    /// True if this fragment originated from mesh_a, false if from mesh_b.
    pub from_a: bool,
}

// Ã¢â€â‚¬Ã¢â€â‚¬ Phase 4: classify and keep fragments Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬

/// Generalized Winding Number (GWN) of `query` with respect to a closed triangle mesh.
///
/// ## Theorem Ã¢â‚¬â€ GWN Inside/Outside
///
/// For a closed orientable 2-manifold, the GWN at any interior point equals Ã‚Â±1
/// and at any exterior point equals 0.  This is computed via the van OosteromÃ¢â‚¬â€œ
/// Strackee (1983) solid-angle formula applied per triangle:
///   ÃŽÂ© = 2Ã‚Â·atan2(aÃ‚Â·(bÃƒâ€”c), |a||b||c| + (aÃ‚Â·b)|c| + (bÃ‚Â·c)|a| + (cÃ‚Â·a)|b|)
/// where a,b,c are the unit vectors from query to the triangle vertices.
///
/// No epsilon displacement is required Ã¢â‚¬â€ GWN is well-defined everywhere except
/// at mesh vertices.  Query points lying exactly on a face plane produce |GWN| Ã¢â€°Ë† 0.5;
/// the `classify_fragment` function uses an exact `orient3d` tiebreaker for these.
///
/// Reference: Jacobson et al. (2013), "Robust Inside-Outside Segmentation using
/// Generalized Winding Numbers", ACM SIGGRAPH.
pub fn gwn(query: &Point3r, faces: &[FaceData], pool: &VertexPool) -> Real {
    let mut solid_angle_sum = 0.0_f64;
    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);

        let va = Vector3r::new(a.x - query.x, a.y - query.y, a.z - query.z);
        let vb = Vector3r::new(b.x - query.x, b.y - query.y, b.z - query.z);
        let vc = Vector3r::new(c.x - query.x, c.y - query.y, c.z - query.z);

        let la = va.norm();
        let lb = vb.norm();
        let lc = vc.norm();

        // Skip faces where query coincides with a vertex.
        if la < 1e-20 || lb < 1e-20 || lc < 1e-20 {
            continue;
        }

        let num = va.dot(&vb.cross(&vc));
        let den = la * lb * lc + va.dot(&vb) * lc + vb.dot(&vc) * la + vc.dot(&va) * lb;

        if den.abs() > 1e-30 || num.abs() > 1e-30 {
            solid_angle_sum += 2.0 * num.atan2(den);
        }
    }
    (solid_angle_sum / (4.0 * std::f64::consts::PI)).clamp(-1.0, 1.0)
}

/// Classify whether a fragment's centroid is inside the opposing mesh.
///
/// ## Algorithm
///
/// 1. Evaluate GWN at the exact centroid (no nudge).
///    - Clearly inside  (|GWN| > 0.65): return `true`.
///    - Clearly outside (|GWN| < 0.35): return `false`.
/// 2. **Exact-predicate tiebreaker** for seam fragments (0.35 Ã¢â€°Â¤ |GWN| Ã¢â€°Â¤ 0.65):
///    The centroid lies on an opposing face plane.  Compare the fragment's
///    outward normal against each coplanar opposing face's CCW normal via
///    exact `orient3d`.  Same-direction normals Ã¢â€ â€™ exterior; opposite Ã¢â€ â€™ interior.
pub fn classify_fragment(
    centroid: &Point3r,
    frag_normal: &Vector3r,
    other_faces: &[FaceData],
    pool: &VertexPool,
) -> bool {
    use crate::domain::topology::predicates::{orient3d, Sign};

    let wn = gwn(centroid, other_faces, pool);
    let wn_abs = wn.abs();
    // Fast path: unambiguous.
    if wn_abs > 0.65 {
        return true;
    }
    if wn_abs < 0.35 {
        return false;
    }

    // Tiebreaker 1: exact normal comparison for seam fragments.
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
        return true;
    }
    if exterior_votes > interior_votes {
        return false;
    }

    // Tiebreaker 2: nearest-face-centroid signed distance.
    // For curved meshes, exact coplanarity (orient3d == Zero) almost never
    // fires.  Instead, find the face whose centroid is closest (L2) to the
    // query centroid and use its plane's signed distance.
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
    // Negative signed distance Ã¢â€¡â€™ centroid is behind the closest face Ã¢â€¡â€™ inside.
    if best_sign.abs() > 1e-15 {
        return best_sign < 0.0;
    }
    // Ultimate fallback when everything is degenerate.
    wn_abs > 0.5
}

/// Triangle centroid.
#[inline]
pub fn centroid(tri: &[Point3r; 3]) -> Point3r {
    Point3r::new(
        (tri[0].x + tri[1].x + tri[2].x) / 3.0,
        (tri[0].y + tri[1].y + tri[2].y) / 3.0,
        (tri[0].z + tri[1].z + tri[2].z) / 3.0,
    )
}

/// Geometric normal of a triangle (not normalised).
#[inline]
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

