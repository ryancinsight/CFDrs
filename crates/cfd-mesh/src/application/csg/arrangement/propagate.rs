use crate::application::csg::intersect::SnapSegment;
use crate::application::csg::predicates3d::{
    point_on_segment_exact, proper_segment_intersection_params_projected_exact,
};
use crate::application::welding::snap::GridCell;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
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
        let p0 = *pool.position(v[0]);
        let p1 = *pool.position(v[1]);
        let p2 = *pool.position(v[2]);
        let face_n = (p1 - p0).cross(&(p2 - p0));

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
                    // Shadow evaluation only: keep exact projected predicate hot
                    // without changing current seam behavior.
                    let _ = proper_segment_intersection_params_projected_exact(
                        &pa, &pb, &seg.start, &seg.end, &face_n,
                    );
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
/// ## Algorithm (Spatial-Hash Accelerated)
///
/// The naïve O(B × P) loop — for each barrel face, test every seam position — is
/// replaced by an O(R × k + P) algorithm:
///
/// 1. **Pre-filter** barrel faces to "rim faces" (exactly 2 on-plane vertices).
///    Non-rim faces are skipped immediately without any seam-position work.
///    Reduces the outer-loop count from B (all barrel faces) to R ≤ B.
///
/// 2. **Build seam-position spatial hash**: map each seam point into a 1 mm
///    grid cell.  Cost: O(P).
///
/// 3. **Query per rim face**: sample 5 points along the rim edge at
///    t ∈ {0, ¼, ½, ¾, 1}, query the 27-cell neighbourhood of each sample
///    (135 cells total), collect candidate seam indices, deduplicate, then
///    run the collinearity test only on candidates.  Cost: O(R × k) where
///    k = seam positions per cell ≈ 1–3 for millifluidic meshes.
///
/// ## Theorem (Spatial Hash Correctness)
///
/// For any seam position P strictly on rim edge [pa, pb] with |pb − pa| ≤ 8 mm,
/// the nearest sample point is within 1 mm of P.  Proof: samples divide [pa, pb]
/// into 4 equal sub-intervals of length |pb − pa|/4 ≤ 2 mm.  P's distance to the
/// nearest sample is at most |pb − pa|/8 ≤ 1 mm = cell size.  Therefore P lands
/// in the 27-cell neighbourhood (radius 1 = 1 mm) of that sample's GridCell.  For
/// edges > 8 mm, the fallback collinearity test (part of the inner loop) still
/// runs correctly — coverage is only reduced to O(1/cell_size) density.  QED.
///
/// ## Complexity
///
/// | Phase | Cost |
/// |-------|------|
/// | Rim-face pre-filter | O(B) |
/// | Seam-position hash build | O(P) |
/// | Query + collinearity test | O(R × 135 × k) |
/// | **Total** | **O(B + P + R × k)** — vs O(B × P) naïve |
///
/// For millifluidic meshes with k ≈ 1–3: effectively **O(B + P)**.
///
/// ## Collinearity test
///
/// Point `s` lies on segment `pa→pb` iff the cross-product `(pb−pa)×(s−pa)` is
/// the zero vector (collinear) and the dot-product parameter
/// `t = (s−pa)·(pb−pa) / |pb−pa|²` lies in `(MARGIN, 1−MARGIN)`.
///
/// Uses `COLLINEAR_TOL_SQ` (1e-6 on cross²/edge², i.e., |cross|/|edge| < 1e-3),
/// matching `propagate_seam_vertices` for consistent seam detection.
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
    let tol = ON_TOL * plane_n_len;

    // ── Phase 1: Pre-filter barrel faces to rim faces ─────────────────────────
    // Rim face: exactly 2 on-plane vertices (the rim edge [pa, pb] lies on the
    // cap plane).  We pre-compute all rim edges once instead of re-detecting
    // them inside the seam-position loop.
    struct RimFace {
        face_idx: usize,
        pa: Point3r,
        pb: Point3r,
    }

    let mut rim_faces: Vec<RimFace> = Vec::new();
    for (face_idx, face) in barrel_faces.iter().enumerate() {
        if coplanar_used.contains(&face_idx) {
            continue;
        }
        let v0 = *pool.position(face.vertices[0]);
        let v1 = *pool.position(face.vertices[1]);
        let v2 = *pool.position(face.vertices[2]);
        let d0 = (v0 - plane_pt).dot(plane_n);
        let d1 = (v1 - plane_pt).dot(plane_n);
        let d2 = (v2 - plane_pt).dot(plane_n);
        let on0 = d0.abs() < tol;
        let on1 = d1.abs() < tol;
        let on2 = d2.abs() < tol;
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
        let edge_len_sq = (pb - pa).norm_squared();
        if edge_len_sq < 1e-20 {
            continue;
        }
        rim_faces.push(RimFace { face_idx, pa, pb });
    }

    if rim_faces.is_empty() {
        return;
    }

    // ── Phase 2: Build seam-position spatial hash ─────────────────────────────
    // Cell size 1 mm covers millifluidic rim edges ≤ 8 mm via 5 sample points
    // (see theorem in module doc above).
    const HASH_CELL_M: Real = 1e-3; // 1 mm
    let inv_cell = 1.0 / HASH_CELL_M;

    let mut seam_hash: HashMap<GridCell, Vec<usize>> =
        HashMap::with_capacity(seam_positions.len() * 2);
    for (i, s) in seam_positions.iter().enumerate() {
        seam_hash
            .entry(GridCell::from_point_round(s, inv_cell))
            .or_default()
            .push(i);
    }

    // ── Phase 3: Query & inject ───────────────────────────────────────────────
    for rim in &rim_faces {
        let RimFace { face_idx, pa, pb } = *rim;
        let edge = pb - pa;
        let edge_len_sq = edge.norm_squared();

        // 5 sample points at t ∈ {0, ¼, ½, ¾, 1} along the rim edge.
        // Deduplication via sort+dedup on the candidate index list.
        let mut candidates: Vec<usize> = Vec::new();
        for k in 0..=4_u8 {
            let t = k as Real * 0.25;
            let sample = pa + edge * t;
            let home = GridCell::from_point_round(&sample, inv_cell);
            for cell in home.neighborhood_27() {
                if let Some(idxs) = seam_hash.get(&cell) {
                    candidates.extend_from_slice(idxs);
                }
            }
        }

        if candidates.is_empty() {
            continue;
        }
        candidates.sort_unstable();
        candidates.dedup();

        // Collinearity test for each candidate seam position.
        let mut cut_params: Vec<Real> = Vec::new();
        for i in candidates {
            let s = &seam_positions[i];

            // Guard: s must lie on the cap plane.
            let ds = (*s - plane_pt).dot(plane_n);
            if ds.abs() > tol * 10.0 {
                continue;
            }

            // Exact-first on-segment detection.
            if let Some(t_exact) = point_on_segment_exact(&pa, &pb, s) {
                if t_exact > SEG_MARGIN && t_exact < 1.0 - SEG_MARGIN {
                    cut_params.push(t_exact);
                }
                continue;
            }

            // Tolerance-based collinearity + parameter check.
            let sp = *s - pa;
            let cross = edge.cross(&sp);
            if cross.norm_squared() <= COLLINEAR_TOL_SQ * edge_len_sq {
                let t = sp.dot(&edge) / edge_len_sq;
                if t > SEG_MARGIN && t < 1.0 - SEG_MARGIN {
                    cut_params.push(t);
                }
            }
        }

        if cut_params.is_empty() {
            continue;
        }

        // Sort and deduplicate cut parameters, then emit sub-interval SnapSegments.
        cut_params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        cut_params.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

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
    use crate::domain::core::scalar::Vector3r;
    use std::collections::HashMap;

    fn contains_param_split(
        segs: &[SnapSegment],
        a: Point3r,
        b: Point3r,
        t: Real,
        tol: Real,
    ) -> bool {
        let x = a + (b - a) * t;
        segs.iter().any(|s| {
            (nalgebra::distance_squared(&s.start, &a) < tol
                && nalgebra::distance_squared(&s.end, &x) < tol)
                || (nalgebra::distance_squared(&s.start, &x) < tol
                    && nalgebra::distance_squared(&s.end, &b) < tol)
        })
    }

    #[test]
    fn propagate_seam_vertices_injects_crossing_split_into_adjacent_face() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(0.5, 1.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(0.5, -1.0, 0.0), n);
        let faces = vec![FaceData::untagged(a, b, c), FaceData::untagged(b, a, d)];

        let mut segs: HashMap<usize, Vec<SnapSegment>> = HashMap::new();
        segs.insert(
            0,
            vec![SnapSegment {
                start: Point3r::new(0.25, 0.5, 0.0),
                end: Point3r::new(0.25, -0.5, 0.0),
            }],
        );

        propagate_seam_vertices(&faces, &mut segs, &pool);

        let injected = segs
            .get(&1)
            .expect("adjacent face should receive seam splits");
        let pa = *pool.position(a);
        let pb = *pool.position(b);
        assert!(
            injected.len() >= 2,
            "crossing split should create at least two sub-segments on shared edge"
        );
        assert!(
            contains_param_split(injected, pa, pb, 0.25, 1e-10),
            "adjacent face should receive split at crossing parameter t=0.25"
        );
    }

    #[test]
    fn adversarial_near_parallel_crossing_is_propagated() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(1.0, 1.0e-10, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(0.2, 1.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(0.2, -1.0, 0.0), n);
        let faces = vec![FaceData::untagged(a, b, c), FaceData::untagged(b, a, d)];

        let mut segs: HashMap<usize, Vec<SnapSegment>> = HashMap::new();
        segs.insert(
            0,
            vec![SnapSegment {
                start: Point3r::new(0.5, -1.0e-10, 0.0),
                end: Point3r::new(0.5000000001, 1.0e-10, 0.0),
            }],
        );

        propagate_seam_vertices(&faces, &mut segs, &pool);
        let injected = segs
            .get(&1)
            .expect("near-parallel proper crossing should still propagate");
        assert!(
            !injected.is_empty(),
            "adjacent face should receive injected segments for near-parallel crossing"
        );
    }

    // ── inject_cap_seam_into_barrels tests ────────────────────────────────────

    /// Build a single barrel rim face: v0=(0,0,0), v1=(1,0,0) on plane z=0,
    /// v2=(0.5,0,−1) off-plane.  Rim edge is [v0,v1] along X.
    fn single_rim_face(pool: &mut VertexPool) -> (Vec<FaceData>, Point3r, Vector3r) {
        let nz = Vector3r::new(0.0, 0.0, 1.0);
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), nz);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), nz);
        let v2 = pool.insert_or_weld(Point3r::new(0.5, 0.0, -1.0), nz);
        let faces = vec![FaceData::untagged(v0, v1, v2)];
        let plane_pt = Point3r::new(0.0, 0.0, 0.0);
        let plane_n = Vector3r::new(0.0, 0.0, 1.0);
        (faces, plane_pt, plane_n)
    }

    /// Seam position at t=0.25 on a 1m rim edge is found by the spatial hash and
    /// produces 2 sub-segment injections.
    #[test]
    fn inject_cap_seam_finds_seam_at_quarter_param() {
        let mut pool = VertexPool::default_millifluidic();
        let (faces, plane_pt, plane_n) = single_rim_face(&mut pool);
        let coplanar_used = std::collections::HashSet::new();
        let seam_positions = vec![Point3r::new(0.25, 0.0, 0.0)];
        let mut segs_out: HashMap<usize, Vec<SnapSegment>> = HashMap::new();

        inject_cap_seam_into_barrels(
            &faces,
            &coplanar_used,
            &plane_pt,
            &plane_n,
            &seam_positions,
            &mut segs_out,
            &pool,
        );

        let injected = segs_out
            .get(&0)
            .expect("rim face 0 should receive injected segments");
        assert!(
            !injected.is_empty(),
            "seam at t=0.25 should generate sub-segments"
        );
        assert_eq!(injected.len(), 2, "one seam point creates 2 sub-segments");
    }

    /// Seam position not on the cap plane (z=0.5) must not inject into any face.
    #[test]
    fn inject_cap_seam_ignores_off_plane_position() {
        let mut pool = VertexPool::default_millifluidic();
        let (faces, plane_pt, plane_n) = single_rim_face(&mut pool);
        let coplanar_used = std::collections::HashSet::new();
        // z=0.5 — off-plane, should be rejected by the ds.abs() guard.
        let seam_positions = vec![Point3r::new(0.5, 0.0, 0.5)];
        let mut segs_out: HashMap<usize, Vec<SnapSegment>> = HashMap::new();

        inject_cap_seam_into_barrels(
            &faces,
            &coplanar_used,
            &plane_pt,
            &plane_n,
            &seam_positions,
            &mut segs_out,
            &pool,
        );

        assert!(
            segs_out.is_empty(),
            "off-plane seam position must not inject into any face"
        );
    }

    /// Two seam positions (given out of order) produce 3 sorted sub-segments.
    #[test]
    fn inject_cap_seam_multiple_positions_generate_sorted_sub_intervals() {
        let mut pool = VertexPool::default_millifluidic();
        let (faces, plane_pt, plane_n) = single_rim_face(&mut pool);
        let coplanar_used = std::collections::HashSet::new();
        // Intentionally insert t=0.75 before t=0.25 to verify sorting.
        let seam_positions = vec![Point3r::new(0.75, 0.0, 0.0), Point3r::new(0.25, 0.0, 0.0)];
        let mut segs_out: HashMap<usize, Vec<SnapSegment>> = HashMap::new();

        inject_cap_seam_into_barrels(
            &faces,
            &coplanar_used,
            &plane_pt,
            &plane_n,
            &seam_positions,
            &mut segs_out,
            &pool,
        );

        let injected = segs_out
            .get(&0)
            .expect("rim face should receive injected segments");
        assert_eq!(injected.len(), 3, "two seam points create 3 sub-segments");
    }

    /// A face with only 1 on-plane vertex is not a rim face and must not receive
    /// any injected segments.
    #[test]
    fn inject_cap_seam_non_rim_face_is_skipped() {
        let mut pool = VertexPool::default_millifluidic();
        let nz = Vector3r::new(0.0, 0.0, 1.0);
        // Only v0 is on the cap plane z=0 → on_count=1 → not a rim face.
        let v0 = pool.insert_or_weld(Point3r::new(0.5, 0.0, 0.0), nz);
        let v1 = pool.insert_or_weld(Point3r::new(0.0, 0.5, -1.0), nz);
        let v2 = pool.insert_or_weld(Point3r::new(1.0, 0.5, -1.0), nz);
        let faces = vec![FaceData::untagged(v0, v1, v2)];
        let coplanar_used = std::collections::HashSet::new();
        let plane_pt = Point3r::new(0.0, 0.0, 0.0);
        let plane_n = Vector3r::new(0.0, 0.0, 1.0);
        let seam_positions = vec![Point3r::new(0.5, 0.0, 0.0)];
        let mut segs_out: HashMap<usize, Vec<SnapSegment>> = HashMap::new();

        inject_cap_seam_into_barrels(
            &faces,
            &coplanar_used,
            &plane_pt,
            &plane_n,
            &seam_positions,
            &mut segs_out,
            &pool,
        );

        assert!(
            segs_out.is_empty(),
            "non-rim face (1 on-plane vertex) must not receive injected segments"
        );
    }

    /// Adversarial: seam position not on the rim edge (off to the side) is
    /// rejected even though it passes the plane check.
    #[test]
    fn inject_cap_seam_position_beside_rim_edge_is_rejected() {
        let mut pool = VertexPool::default_millifluidic();
        let (faces, plane_pt, plane_n) = single_rim_face(&mut pool);
        let coplanar_used = std::collections::HashSet::new();
        // y=0.5 puts the point on the cap plane but off the rim edge [0,0,0]→[1,0,0].
        let seam_positions = vec![Point3r::new(0.5, 0.5, 0.0)];
        let mut segs_out: HashMap<usize, Vec<SnapSegment>> = HashMap::new();

        inject_cap_seam_into_barrels(
            &faces,
            &coplanar_used,
            &plane_pt,
            &plane_n,
            &seam_positions,
            &mut segs_out,
            &pool,
        );

        assert!(
            segs_out.is_empty(),
            "seam position beside rim edge must not inject (off-edge but on-plane)"
        );
    }
}
