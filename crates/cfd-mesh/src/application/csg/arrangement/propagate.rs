use crate::application::csg::intersect::SnapSegment;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
use std::collections::HashMap;

/// Collinearity tolerance for point-on-edge detection in seam propagation.
///
/// A point P is collinear with edge [Va, Vb] if:
///   `|cross(Vb-Va, P-Va)|² < COLLINEAR_TOL_SQ * |Vb-Va|²`
///
/// Equivalently: `sin(angle) < sqrt(COLLINEAR_TOL_SQ) ≈ 1e-3` (0.06°).
///
/// Widened from 1e-8 to 1e-6 to handle shallow-angle sliver faces from
/// elbow-cylinder junctions where edge_len_sq is very small.  Used by both
/// `propagate_seam_vertices` and `inject_cap_seam_into_barrels` to ensure
/// consistent seam vertex detection across the pipeline.
const COLLINEAR_TOL_SQ: Real = 1e-6;

#[inline]
fn collinear_3d_exact(a: &Point3r, b: &Point3r, p: &Point3r) -> bool {
    orient_2d_arr([a.x, a.y], [b.x, b.y], [p.x, p.y]) == Orientation::Degenerate
        && orient_2d_arr([a.x, a.z], [b.x, b.z], [p.x, p.z]) == Orientation::Degenerate
        && orient_2d_arr([a.y, a.z], [b.y, b.z], [p.y, p.z]) == Orientation::Degenerate
}

/// Exact-first point-on-segment parameter.
///
/// ## Algorithm
///
/// 1. Verify exact 3-D collinearity by requiring all three projection
///    orientations (`xy`, `xz`, `yz`) to be degenerate.
/// 2. Compute segment parameter `t = dot(p-a, b-a)/|b-a|²`.
/// 3. Accept only strict interior points (`0 < t < 1`).
///
/// ## Theorem — 3-D Collinearity by Projection
///
/// A point `p` is collinear with segment endpoints `a,b` in 3-D iff the three
/// independent 2-D projections `(xy, xz, yz)` are all collinear.
///
/// **Proof sketch.**
/// `p-a` and `b-a` are linearly dependent in R³ iff all pairwise 2x2 minors of
/// the matrix `[b-a, p-a]` vanish. Those minors are exactly the determinants
/// tested by projected orientation predicates on `xy`, `xz`, and `yz`. ∎
#[inline]
fn point_on_segment_exact(a: &Point3r, b: &Point3r, p: &Point3r) -> Option<Real> {
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

/// Ensure that every seam vertex created by CDT co-refinement is injected into
/// all faces that share the face edge on which the seam vertex lies.
///
/// ## Problem — T-junctions at shared edges
///
/// When `intersect_triangles` produces a snap-segment endpoint P that lies
/// on the boundary edge `[Va, Vb]` of a triangle face `f1`, the CDT of `f1`
/// inserts a Steiner vertex at P and produces sub-edges `Va→P` and `P→Vb`.
/// However, the adjacent face `f2` (which shares the undirected edge `{Va,Vb}`)
/// has no snap segment touching P, so its CDT leaves edge `Va→Vb` unsplit.
///
/// In the final mesh, sub-edge `Va→P` appears once (from `f1`'s CDT) with no
/// counterpart from `f2` → open boundary edge → non-manifold output.
///
/// ## Algorithm
///
/// 1. Build an undirected edge adjacency map: `{Va, Vb} → [face_idx, …]`.
/// 2. For each face `f` that has snap segments, collect all seam-endpoint
///    positions from those segments.
/// 3. For each endpoint P and each edge `[Va, Vb]` of face `f`, check if P
///    lies strictly between Va and Vb (collinearity + parameter check).
/// 4. If yes, inject snap segments `Va→P` and `P→Vb` into every OTHER face
///    that shares edge `{Va, Vb}` — propagating the Steiner vertex across the
///    shared edge.
///
/// ## Collinearity threshold
///
/// A point P is considered to lie on edge [Va, Vb] if:
/// - `|(Vb−Va) × (P−Va)|² / |Vb−Va|² < COLLINEAR_TOL_SQ` (1e-6; sub-millimetre)
/// - parameter `t = (P−Va)·(Vb−Va) / |Vb−Va|² ∈ (1e-7, 1−1e-7)`
pub fn propagate_seam_vertices(
    faces: &[FaceData],
    segs: &mut HashMap<usize, Vec<SnapSegment>>,
    pool: &VertexPool,
) {
    use crate::domain::core::index::VertexId;

    if segs.is_empty() {
        return;
    }

    // Build undirected edge → face-index adjacency.
    type EdgeKey = (VertexId, VertexId);
    let mut edge_to_faces: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
    for (fi, face) in faces.iter().enumerate() {
        let v = face.vertices;
        for i in 0..3_usize {
            let va = v[i];
            let vb = v[(i + 1) % 3];
            let key = if va < vb { (va, vb) } else { (vb, va) };
            edge_to_faces.entry(key).or_default().push(fi);
        }
    }

    let mut injections: Vec<(usize, SnapSegment)> = Vec::new();

    for (&fi, snap_segs) in segs.iter() {
        if fi >= faces.len() {
            continue;
        }
        let face = &faces[fi];
        let v = face.vertices;

        // For each face edge, collect all positions where snap segments touch or cross it.
        // These include:
        //   (A) snap-segment endpoints that lie on the edge
        //   (B) crossing points where a snap segment crosses the edge interior
        // Both types create Steiner vertices in the CDT that must be propagated
        // to the adjacent face sharing that edge.
        for i in 0..3_usize {
            let va_id = v[i];
            let vb_id = v[(i + 1) % 3];
            let pa = *pool.position(va_id);
            let pb = *pool.position(vb_id);

            let edge_vec = pb - pa;
            let edge_len_sq = edge_vec.dot(&edge_vec);
            if edge_len_sq < 1e-20 {
                continue;
            }

            let edge_key = if va_id < vb_id {
                (va_id, vb_id)
            } else {
                (vb_id, va_id)
            };
            let adj_faces = match edge_to_faces.get(&edge_key) {
                Some(f) => f,
                None => continue,
            };
            if !adj_faces.iter().any(|&f| f != fi) {
                continue;
            }

            // Collect t-parameters (on edge [pa,pb], t ∈ (0,1)) for all contacts.
            let mut t_params: Vec<Real> = Vec::new();
            const MARGIN: Real = 1e-7;

            for seg in snap_segs {
                // (A) Endpoint on edge.
                for &p in &[seg.start, seg.end] {
                    if let Some(t_exact) = point_on_segment_exact(&pa, &pb, &p) {
                        if t_exact > MARGIN && t_exact < 1.0 - MARGIN {
                            t_params.push(t_exact);
                        }
                        continue;
                    }

                    // Fallback: tolerance-based on-edge check for residual drift.
                    let sp: nalgebra::Vector3<f64> = p - pa;
                    let cross_v = edge_vec.cross(&sp);
                    if cross_v.norm_squared() <= COLLINEAR_TOL_SQ * edge_len_sq {
                        let t = sp.dot(&edge_vec) / edge_len_sq;
                        if t > MARGIN && t < 1.0 - MARGIN {
                            t_params.push(t);
                        }
                    }
                }

                // (B) Segment-edge crossing in 3-D.
                // Solve pa + t*(pb-pa) = seg.start + s*(seg.end-seg.start):
                // Use the two axis-equations with the largest determinant.
                let sv = seg.end - seg.start;
                let r_vec = seg.start - pa;
                let pairs: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];
                let mut best_det_abs = 0.0_f64;
                let mut best_t = 0.0_f64;
                let mut best_s = 0.0_f64;
                for &(ax, ay) in &pairs {
                    // e[0]*t - sv[0]*s = r[0]
                    // e[1]*t - sv[1]*s = r[1]
                    let e0 = edge_vec[ax];
                    let e1 = edge_vec[ay];
                    let s0 = sv[ax];
                    let s1 = sv[ay];
                    let r0 = r_vec[ax];
                    let r1 = r_vec[ay];
                    let det = e0 * (-s1) - e1 * (-s0); // det([e,-sv])
                    if det.abs() > best_det_abs {
                        best_det_abs = det.abs();
                        best_t = (r0 * (-s1) - r1 * (-s0)) / det;
                        best_s = (e0 * r1 - e1 * r0) / det;
                    }
                }
                let min_det = 1e-14 * (edge_len_sq + sv.norm_squared()).sqrt();
                if best_det_abs < min_det {
                    continue;
                }
                if best_t <= MARGIN || best_t >= 1.0 - MARGIN {
                    continue;
                }
                if best_s <= MARGIN || best_s >= 1.0 - MARGIN {
                    continue;
                }
                // Verify that the crossing is consistent (lines actually meet).
                let x_edge = pa + edge_vec * best_t;
                let x_seg = seg.start + sv * best_s;
                // Widened from 1e-8 to 1e-6 to match WELD_TOL_SQ in corefine.rs.
                if nalgebra::distance_squared(&x_edge, &x_seg) > 1e-6 * edge_len_sq {
                    continue;
                }
                t_params.push(best_t);
            }

            if t_params.is_empty() {
                continue;
            }

            // Deduplicate and sort t-parameters.
            t_params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            t_params.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

            // Build sub-interval snap segments and inject into adjacent faces.
            let pts: Vec<Point3r> = std::iter::once(pa)
                .chain(t_params.iter().map(|&t| pa + edge_vec * t))
                .chain(std::iter::once(pb))
                .collect();

            for &adj_fi in adj_faces {
                if adj_fi == fi {
                    continue;
                }
                for w in pts.windows(2) {
                    if (w[1] - w[0]).norm_squared() < 1e-20 {
                        continue;
                    }
                    injections.push((
                        adj_fi,
                        SnapSegment {
                            start: w[0],
                            end: w[1],
                        },
                    ));
                }
            }
        }
    }

    for (fi, seg) in injections {
        segs.entry(fi).or_default().push(seg);
    }
}

// ══ Phase 2d helper: seam vertex injection into barrel rim faces ═══════════════

/// Inject snap segments into barrel rim faces so they are corefined at every
/// seam vertex produced by `boolean_coplanar`.
///
/// ## Problem
///
/// `boolean_coplanar` clips cap triangles against each other in 2-D, producing
/// NEW intersection vertices (e.g., where an A-cap interior edge crosses a B-cap
/// interior edge).  These vertices appear on the *cap* side of the mesh but NOT
/// on the adjacent barrel face's rim edge, creating T-junctions → boundary edges
/// in the output.
///
/// ## Algorithm
///
/// For each barrel face with exactly 2 on-plane vertices (the "rim edge" `pa→pb`),
/// iterate over every seam position `s` produced by `boolean_coplanar`.  If `s`
/// lies strictly between `pa` and `pb` on the rim edge (collinearity + parameter
/// check), inject a `SnapSegment` covering the full rim edge through `s`.
///
/// The injected `SnapSegment`s cause `corefine_face` (Phase 3) to subdivide the
/// barrel face at exactly the same vertex position as the cap triangle, so both
/// sides produce matching edges → watertight seam.
///
/// ## Collinearity test
///
/// Point `s` lies on segment `pa→pb` iff the cross-product `(pb−pa)×(s−pa)` is
/// the zero vector (collinear) and the dot-product parameter
/// `t = (s−pa)·(pb−pa) / |pb−pa|²` lies in `(MARGIN, 1−MARGIN)`.
///
/// To handle floating-point imprecision from the 2-D clipping step we use
/// `COLLINEAR_TOL_SQ` (1e-6 on cross²/edge², i.e., |cross|/|edge| < 1e-3),
/// matching `propagate_seam_vertices` to ensure consistent seam detection.
pub fn inject_cap_seam_into_barrels(
    barrel_faces: &[FaceData],
    coplanar_used: &std::collections::HashSet<usize>,
    plane_pt: &Point3r,
    plane_n: &Vector3r,
    seam_positions: &[Point3r],
    segs_out: &mut HashMap<usize, Vec<SnapSegment>>,
    pool: &VertexPool,
) {
    let plane_n_len_sq = plane_n.dot(plane_n);
    if plane_n_len_sq < 1e-20 || seam_positions.is_empty() {
        return;
    }
    let plane_n_len = plane_n_len_sq.sqrt();

    const ON_TOL: Real = 1e-7; // signed-distance tolerance (relative to normal length)
    const SEG_MARGIN: Real = 1e-7; // parameter margin for "strictly interior"

    for (face_idx, face) in barrel_faces.iter().enumerate() {
        if coplanar_used.contains(&face_idx) {
            continue;
        }

        let v0 = *pool.position(face.vertices[0]);
        let v1 = *pool.position(face.vertices[1]);
        let v2 = *pool.position(face.vertices[2]);

        // Signed distance of each vertex to the cap plane (unnormalised).
        let d0 = (v0 - plane_pt).dot(plane_n);
        let d1 = (v1 - plane_pt).dot(plane_n);
        let d2 = (v2 - plane_pt).dot(plane_n);

        let tol = ON_TOL * plane_n_len;
        let on0 = d0.abs() < tol;
        let on1 = d1.abs() < tol;
        let on2 = d2.abs() < tol;

        // A rim barrel face has exactly 2 on-plane vertices.
        let on_count = u8::from(on0) + u8::from(on1) + u8::from(on2);
        if on_count != 2 {
            continue;
        }

        let (pa, pb) = match (on0, on1, on2) {
            (true, true, false) => (v0, v1),
            (true, false, true) => (v0, v2),
            (false, true, true) => (v1, v2),
            _ => continue,
        };

        // Rim edge vector and squared length.
        let edge = pb - pa;
        let edge_len_sq = edge.dot(&edge);
        if edge_len_sq < 1e-20 {
            continue;
        }

        // For each seam position, check if it lies strictly on the rim edge pa→pb.
        let mut cut_params: Vec<Real> = Vec::new();

        for s in seam_positions {
            // (1) Check that s is on the cap plane (should always be true, but guard anyway).
            let ds = (*s - plane_pt).dot(plane_n);
            if ds.abs() > tol * 10.0 {
                continue;
            }

            // (2)+(3): exact-first on-segment detection.
            if let Some(t_exact) = point_on_segment_exact(&pa, &pb, s) {
                if t_exact > SEG_MARGIN && t_exact < 1.0 - SEG_MARGIN {
                    cut_params.push(t_exact);
                }
                continue;
            }

            // Fallback: tolerance-based collinearity + parameter check.
            let sp = *s - pa;
            let cross = edge.cross(&sp);
            let cross_len_sq = cross.dot(&cross);
            if cross_len_sq <= COLLINEAR_TOL_SQ * edge_len_sq {
                let t = sp.dot(&edge) / edge_len_sq;
                if t > SEG_MARGIN && t < 1.0 - SEG_MARGIN {
                    cut_params.push(t);
                }
            }
        }

        if cut_params.is_empty() {
            continue;
        }

        // Deduplicate and sort cut parameters.
        cut_params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        cut_params.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

        // Build sub-intervals [0, t0, t1, …, 1] and emit each as a SnapSegment.
        // corefine_face will insert the interior break-points as constrained CDT
        // vertices, forcing the mesh to split at those positions.
        let mut params: Vec<Real> = Vec::with_capacity(cut_params.len() + 2);
        params.push(0.0);
        params.extend_from_slice(&cut_params);
        params.push(1.0);

        for w in params.windows(2) {
            let (t0, t1) = (w[0], w[1]);
            if (t1 - t0).abs() < 1e-12 {
                continue;
            }
            let start_3d = pa + edge * t0;
            let end_3d = pa + edge * t1;
            if (end_3d - start_3d).norm_squared() < 1e-20 {
                continue;
            }
            segs_out.entry(face_idx).or_default().push(SnapSegment {
                start: start_3d,
                end: end_3d,
            });
        }
    }
}

// ══ Public entry point ═════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: Real, y: Real, z: Real) -> Point3r {
        Point3r::new(x, y, z)
    }

    #[test]
    fn point_on_segment_exact_accepts_true_3d_collinear_point() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 2.0, 3.0);
        let mid = p(0.5, 1.0, 1.5);
        let t = point_on_segment_exact(&a, &b, &mid);
        assert!(t.is_some(), "true 3D collinear midpoint must be accepted");
    }

    #[test]
    fn point_on_segment_exact_rejects_projection_false_positive() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 1.0, 1.0);
        let q = p(0.5, 0.5, 0.6); // collinear in XY only
        assert!(
            point_on_segment_exact(&a, &b, &q).is_none(),
            "non-collinear 3D point must be rejected"
        );
    }
}
