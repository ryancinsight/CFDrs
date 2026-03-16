//! Boundary-hole patching for arrangement output.
//!
//! After co-refinement, fragment classification, and seam stitching, the
//! output mesh may still contain small boundary holes at barrel–barrel
//! junctions or where excluded coplanar faces leave dangling edges.
//! This module traces those boundary loops and fills them with patch
//! triangles to restore a watertight surface.
//!
//! ## Algorithm — Boundary Hole Patching (7-step pipeline)
//!
//! 1. **Degenerate face removal** — Remove zero-area and extreme-sliver
//!    faces (aspect ratio below threshold).
//!
//! 2. **Face deduplication** — Remove coincident faces with opposing
//!    orientation (internal membrane cancellation).
//!
//! 3. **Conservative seam stitch** — Merge near-coincident boundary
//!    vertices using spatial hashing with tolerance welding.
//!
//! 4. **Snap-round T-junction repair** — Inject boundary vertices that
//!    lie on edges of adjacent faces, closing T-junctions.
//!
//! 5. **Boundary loop tracing** — Trace connected boundary half-edges
//!    into closed loops using a successor map.
//!
//! 6. **Collinear loop filtering** — Skip degenerate loops whose vertices
//!    are nearly collinear (area²/diameter⁴ below threshold).
//!
//! 7. **CDT / ear-clip fill** — Fill each valid loop with constrained
//!    Delaunay triangulation (preferred) or ear-clipping (fallback).
//!
//! ## Theorem — Boundary Loop Tracing Correctness
//!
//! A boundary half-edge `(a → b)` is one whose reverse `(b → a)` does not
//! appear in any face.  The successor map `next: a → b` for boundary
//! edges forms a set of disjoint directed cycles if and only if every
//! boundary vertex has exactly one incoming and one outgoing boundary edge.
//!
//! *Proof.*  Each boundary vertex `v` has `in-degree(v)` = number of
//! boundary edges ending at `v` and `out-degree(v)` = number starting at
//! `v`.  In a manifold mesh, each boundary vertex has exactly one incoming
//! and one outgoing boundary edge (Euler characteristic argument), so the
//! successor map is a permutation, which decomposes into disjoint cycles.  ∎
//!
//! ## References
//!
//! - Liepa, P. (2003). "Filling holes in meshes." *Eurographics Symposium
//!   on Geometry Processing*, 200–205.
//! - Held, M. (2001). "FIST: Fast industrial-strength triangulation of
//!   polygons." *Algorithmica*, 30(4), 563–596.

use hashbrown::{HashMap, HashSet};

use super::mesh_ops::{apply_vertex_merge, boundary_half_edges, dedup_faces_unordered, merge_root};
use super::seam::stitch_boundary_seams_conservative;
use super::snap_round;
use super::stitch;
#[cfg(test)]
use crate::application::csg::diagnostics::trace_enabled;
use crate::application::csg::predicates3d::triangle_is_degenerate_exact;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Real;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Detect small boundary loops (`<= MAX_PATCH_LOOP` edges) in `faces` and
/// fill each with patch triangles.
///
/// A boundary half-edge is one that appears in exactly one face's directed edge
/// set and has no matching reverse. Small boundary loops arise at barrel-barrel
/// junction boundaries where an excluded face's grid edge is not shared by the
/// other mesh's kept fragments.
pub(crate) fn patch_small_boundary_holes(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    const MAX_PATCH_LOOP: usize = 256;
    // (2e-3)^2 -- spatial tolerance for near-duplicate boundary vertex merging.
    // Widened from 4e-8 (=2e-4 mm) to 4e-6 (=2e-3 mm) to match corefine.rs
    // WELD_TOL_SQ, ensuring arithmetic-drift Steiner vertices from shallow-angle
    // elbow-cylinder junctions are welded during the patch pass.
    const BOUNDARY_MERGE_TOL_SQ: Real = 4e-6;

    // Collinear loop threshold: area^2 / diameter^4 < this -> degenerate.
    const COLLINEAR_THRESH: Real = 1e-8;

    // -- Step 1: Remove degenerate faces (zero-area or extreme slivers). ------
    faces.retain(|f| {
        let p0 = pool.position(f.vertices[0]);
        let p1 = pool.position(f.vertices[1]);
        let p2 = pool.position(f.vertices[2]);
        let e01 = p1 - p0;
        let e02 = p2 - p0;
        let e12 = p2 - p1;
        let area_sq = e01.cross(&e02).norm_squared();
        let max_edge_sq = e01
            .norm_squared()
            .max(e02.norm_squared())
            .max(e12.norm_squared());
        area_sq > 1e-8 * max_edge_sq
    });

    // -- Step 2: Remove duplicate faces (same 3 vertex IDs in any order). ----
    dedup_faces_unordered(faces);

    // -- Step 3: Non-manifold edge repair -- deterministic. -------------------
    // For each directed half-edge shared by multiple faces, keep only the
    // face with the largest area (smallest is likely a sliver from CDT).
    {
        let mut he_faces: HashMap<(VertexId, VertexId), Vec<usize>> = HashMap::new();
        for (fi, face) in faces.iter().enumerate() {
            let v = face.vertices;
            for i in 0..3 {
                let j = (i + 1) % 3;
                he_faces.entry((v[i], v[j])).or_default().push(fi);
            }
        }

        let mut remove_set: HashSet<usize> = HashSet::new();
        let mut nm_edges: Vec<(VertexId, VertexId)> = he_faces
            .iter()
            .filter(|(_, v)| v.len() > 1)
            .map(|(&k, _)| k)
            .collect();
        nm_edges.sort();

        for he in &nm_edges {
            let face_indices = &he_faces[he];
            let mut best_fi = face_indices[0];
            let mut best_area = {
                let f = &faces[best_fi];
                let p0 = pool.position(f.vertices[0]);
                let p1 = pool.position(f.vertices[1]);
                let p2 = pool.position(f.vertices[2]);
                (p1 - p0).cross(&(p2 - p0)).norm_squared()
            };
            for &fi in &face_indices[1..] {
                let f = &faces[fi];
                let p0 = pool.position(f.vertices[0]);
                let p1 = pool.position(f.vertices[1]);
                let p2 = pool.position(f.vertices[2]);
                let area_sq = (p1 - p0).cross(&(p2 - p0)).norm_squared();
                if area_sq > best_area {
                    remove_set.insert(best_fi);
                    best_fi = fi;
                    best_area = area_sq;
                } else {
                    remove_set.insert(fi);
                }
            }
        }

        if !remove_set.is_empty() {
            let mut idx = 0;
            faces.retain(|_| {
                let keep = !remove_set.contains(&idx);
                idx += 1;
                keep
            });
        }
    }

    // -- Step 3.5: repair boundary edges exposed by cleanup. ------------------
    // Steps 1-3 remove slivers, duplicates, and non-manifold faces that were
    // masking boundary edges. Now the true boundary topology is visible.
    //
    // Two-pass repair:
    //  (a) Fill closed boundary loops with ear-clipping triangulation.
    //  (b) Conservative short-edge collapse for residual seam gaps.
    //  (c) Another fill pass for loops created by the collapse.
    stitch::fill_boundary_loops(faces, pool);
    stitch_boundary_seams_conservative(faces, pool);
    stitch::fill_boundary_loops(faces, pool);

    // -- Helper: build sorted boundary-edge list from current face set. -------
    let build_boundary = |faces: &[FaceData]| -> Vec<(VertexId, VertexId)> {
        let mut bnd = boundary_half_edges(faces);
        bnd.sort_unstable();
        bnd
    };

    // -- Helper: check if a polygon (by vertex id list) is collinear. ---------
    let is_collinear = |poly: &[VertexId]| -> bool {
        let n = poly.len();
        if n < 3 {
            return true;
        }
        let p0 = *pool.position(poly[0]);
        let mut exact_collinear = true;
        'exact: for i in 1..n {
            let pi = *pool.position(poly[i]);
            for j in (i + 1)..n {
                let pj = *pool.position(poly[j]);
                if !triangle_is_degenerate_exact(&p0, &pi, &pj) {
                    exact_collinear = false;
                    break 'exact;
                }
            }
        }
        if exact_collinear {
            return true;
        }

        // Near-collinear fallback for small residual numerical drift.
        let mut max_area_sq: Real = 0.0;
        let mut diameter_sq: Real = 0.0;
        for &vi in &poly[1..] {
            let pi = pool.position(vi);
            let d = (pi - p0).norm_squared();
            if d > diameter_sq {
                diameter_sq = d;
            }
        }
        for i in 1..n {
            for j in (i + 1)..n {
                let pi = pool.position(poly[i]);
                let pj = pool.position(poly[j]);
                let a = (pi - p0).cross(&(pj - p0)).norm_squared();
                if a > max_area_sq {
                    max_area_sq = a;
                }
            }
        }
        max_area_sq < COLLINEAR_THRESH * diameter_sq * diameter_sq
    };

    // -- Helper: trace closed boundary loops <= MAX_PATCH_LOOP. ---------------
    //
    // Handles figure-8 boundary graphs where a single vertex has boundary
    // out-degree > 1. When the DFS trace reaches a vertex already in the path:
    // 1) extract the inner cycle, 2) truncate path to that vertex, 3) continue.
    let trace_loops = |boundary_edges: &[(VertexId, VertexId)]| -> Vec<Vec<VertexId>> {
        let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
        for &(vi, vj) in boundary_edges {
            adj.entry(vi).or_default().push(vj);
        }
        for v in adj.values_mut() {
            v.sort();
        }

        let mut used_edges: HashSet<(VertexId, VertexId)> = HashSet::new();
        let mut loops: Vec<Vec<VertexId>> = Vec::new();
        let mut starts: Vec<VertexId> = adj.keys().copied().collect();
        starts.sort();

        for start in starts {
            let Some(successors) = adj.get(&start) else {
                continue;
            };
            for &next in successors {
                if used_edges.contains(&(start, next)) {
                    continue;
                }
                let mut loop_verts: Vec<VertexId> = vec![start, next];
                used_edges.insert((start, next));
                let mut cur = next;
                let mut closed = false;
                'trace: loop {
                    if loop_verts.len() > MAX_PATCH_LOOP * 4 {
                        break;
                    }
                    let nexts = match adj.get(&cur) {
                        Some(s) => s,
                        None => break,
                    };
                    let mut found = false;
                    for &n in nexts {
                        if used_edges.contains(&(cur, n)) {
                            continue;
                        }
                        used_edges.insert((cur, n));
                        if n == start {
                            closed = true;
                            break 'trace;
                        }
                        if let Some(pos) = loop_verts.iter().position(|&v| v == n) {
                            let inner = loop_verts[pos..].to_vec(); // [n, ..., cur]
                            if inner.len() >= 3 && inner.len() <= MAX_PATCH_LOOP {
                                loops.push(inner);
                            }
                            loop_verts.truncate(pos + 1);
                            cur = n;
                            found = true;
                            break;
                        }
                        loop_verts.push(n);
                        cur = n;
                        found = true;
                        break;
                    }
                    if !found {
                        break;
                    }
                }
                if closed && loop_verts.len() >= 3 && loop_verts.len() <= MAX_PATCH_LOOP {
                    loops.push(loop_verts);
                }
            }
        }
        loops
    };

    // -- Iterative patching loop. ---------------------------------------------
    // Each iteration:
    //   (a) Build boundary edges.
    //   (b) Merge near-duplicate boundary vertices (Step 5).
    //   (c) Rebuild boundary, trace loops.
    //   (d) Collapse collinear degenerate loops (Step 6).
    //   (e) Fill remaining non-degenerate loops (Step 7).
    const MAX_ITERS: usize = 16;
    for _iter in 0..MAX_ITERS {
        let boundary_edges = build_boundary(faces);
        if boundary_edges.is_empty() {
            break;
        }

        // -- (b) Step 5: exact constrained insertion first, tolerance fallback.
        let before_split = faces.len();
        snap_round::snap_round_tjunctions(faces, pool);
        if faces.len() == before_split {
            let mut bnd_verts: Vec<VertexId> =
                boundary_edges.iter().flat_map(|&(a, b)| [a, b]).collect();
            bnd_verts.sort();
            bnd_verts.dedup();

            let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
            for i in 0..bnd_verts.len() {
                if merge_map.contains_key(&bnd_verts[i]) {
                    continue;
                }
                let pi = pool.position(bnd_verts[i]);
                for j in (i + 1)..bnd_verts.len() {
                    if merge_map.contains_key(&bnd_verts[j]) {
                        continue;
                    }
                    let pj = pool.position(bnd_verts[j]);
                    if (pj - pi).norm_squared() < BOUNDARY_MERGE_TOL_SQ {
                        merge_map.insert(bnd_verts[j], bnd_verts[i]);
                    }
                }
            }
            apply_vertex_merge(faces, &merge_map, pool);
        }

        // -- (c) Rebuild boundary, trace loops. -------------------------------
        let boundary_edges = build_boundary(faces);
        if boundary_edges.is_empty() {
            break;
        }

        #[cfg(test)]
        if trace_enabled() {
            eprintln!("[patch-iter] {} boundary edges", boundary_edges.len());
        }

        let loops = trace_loops(&boundary_edges);

        // -- (d) Step 6: collapse collinear degenerate loops. -----------------
        {
            let mut all_face_verts: Vec<VertexId> = faces
                .iter()
                .flat_map(|f| f.vertices.iter().copied())
                .collect();
            all_face_verts.sort();
            all_face_verts.dedup();

            let mut global_merge: HashMap<VertexId, VertexId> = HashMap::new();

            for poly in &loops {
                let n = poly.len();
                if n < 3 {
                    continue;
                }
                if !is_collinear(poly) {
                    continue;
                }

                let mut best_dist_sq: Real = 0.0;
                let mut endpoint_a = poly[0];
                let mut endpoint_b = poly[1];
                for i in 0..n {
                    let pi = pool.position(poly[i]);
                    for j in (i + 1)..n {
                        let pj = pool.position(poly[j]);
                        let d = (pj - pi).norm_squared();
                        if d > best_dist_sq {
                            best_dist_sq = d;
                            endpoint_a = poly[i];
                            endpoint_b = poly[j];
                        }
                    }
                }
                let pa = pool.position(endpoint_a);
                let pb = pool.position(endpoint_b);

                for &vi in poly {
                    if vi == endpoint_a || vi == endpoint_b {
                        continue;
                    }
                    let pv = pool.position(vi);
                    let da = (pv - pa).norm_squared();
                    let db = (pv - pb).norm_squared();
                    let target = if da <= db { endpoint_a } else { endpoint_b };
                    let final_target = merge_root(&global_merge, target);
                    global_merge.entry(vi).or_insert(final_target);

                    let pvi = pool.position(vi);
                    for &vw in &all_face_verts {
                        if vw == vi || vw == endpoint_a || vw == endpoint_b {
                            continue;
                        }
                        if global_merge.contains_key(&vw) {
                            continue;
                        }
                        let pw = pool.position(vw);
                        if (pw - pvi).norm_squared() < BOUNDARY_MERGE_TOL_SQ {
                            global_merge.insert(vw, final_target);
                        }
                    }
                }
            }
            apply_vertex_merge(faces, &global_merge, pool);
        }

        // -- (e) Step 7: fill non-degenerate loops. ---------------------------
        let boundary_edges_after_collapse = build_boundary(faces);
        if boundary_edges_after_collapse.is_empty() {
            break;
        }
        let loops_after = trace_loops(&boundary_edges_after_collapse);

        let mut valence = stitch::build_canonical_valence(faces);
        let mut any_patch = false;
        for poly in &loops_after {
            if poly.len() < 3 {
                continue;
            }
            if is_collinear(poly) {
                continue;
            }
            // Prefer CDT loop fill (exact predicates) and fall back to ear clip.
            let added = {
                let cdt_added = stitch::cdt_fill_loop(poly, pool, faces, &mut valence);
                if cdt_added > 0 {
                    cdt_added
                } else {
                    stitch::ear_clip_fill(poly, pool, faces, &mut valence)
                }
            };
            if added > 0 {
                any_patch = true;
            }
        }

        if !any_patch {
            break;
        }
    }

    // -- Final cleanup: duplicate removal. ------------------------------------
    dedup_faces_unordered(faces);

    #[cfg(test)]
    {
        let remain = build_boundary(faces);
        if trace_enabled() {
            eprintln!("[patch-done] {} boundary edges remain", remain.len());
            for (a, b) in &remain {
                let pa = pool.position(*a);
                let pb = pool.position(*b);
                eprintln!(
                    "  {:?}->{:?}  ({:.9},{:.9},{:.9})->({:.9},{:.9},{:.9})",
                    a, b, pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};

    /// Helper: build a VertexPool and insert vertices at given positions.
    fn pool_with_positions(pts: &[Point3r]) -> (VertexPool, Vec<VertexId>) {
        let mut pool = VertexPool::new(1e-6_f64);
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let ids: Vec<VertexId> = pts.iter().map(|&p| pool.insert_or_weld(p, n)).collect();
        (pool, ids)
    }

    fn p(x: Real, y: Real, z: Real) -> Point3r {
        Point3r::new(x, y, z)
    }

    // ── Empty / no-op cases ───────────────────────────────────────────────

    /// Patching an empty face list must not panic and must leave faces empty.
    #[test]
    fn patch_empty_faces() {
        let pool = VertexPool::new(1e-6_f64);
        let mut faces: Vec<FaceData> = Vec::new();
        patch_small_boundary_holes(&mut faces, &pool);
        assert!(faces.is_empty());
    }

    /// A single isolated triangle has three boundary edges. The patcher should
    /// not panic but also cannot close such a loop (it needs matching reverse
    /// edges from adjacent faces). The face list must remain non-empty.
    #[test]
    fn patch_single_triangle_survives() {
        let (pool, v) = pool_with_positions(&[p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)]);
        let mut faces = vec![FaceData::untagged(v[0], v[1], v[2])];
        patch_small_boundary_holes(&mut faces, &pool);
        assert!(!faces.is_empty(), "single triangle must not be deleted");
    }

    // ── Already-closed mesh ──────────────────────────────────────────────

    /// A watertight tetrahedron has no boundary edges. Patching should be a
    /// no-op and preserve all 4 faces.
    #[test]
    fn patch_closed_tetrahedron_is_noop() {
        let (pool, v) = pool_with_positions(&[
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.5, 1.0, 0.0),
            p(0.5, 0.5, 1.0),
        ]);
        // 4-face closed tetrahedron (consistent CCW winding from outside)
        let mut faces = vec![
            FaceData::untagged(v[0], v[2], v[1]), // bottom (viewed from -Z)
            FaceData::untagged(v[0], v[1], v[3]), // front
            FaceData::untagged(v[1], v[2], v[3]), // right
            FaceData::untagged(v[2], v[0], v[3]), // left
        ];
        let before = faces.len();
        patch_small_boundary_holes(&mut faces, &pool);
        assert_eq!(
            faces.len(),
            before,
            "closed tetrahedron must not gain or lose faces"
        );
        // Verify still closed
        let boundary = boundary_half_edges(&faces);
        assert!(boundary.is_empty(), "tetrahedron must remain watertight");
    }

    // ── Degenerate face removal ──────────────────────────────────────────

    /// Step 1 of patch removes degenerate (zero-area) faces. A face with two
    /// coincident vertices should be cleaned out.
    #[test]
    fn patch_removes_degenerate_faces() {
        let (pool, v) = pool_with_positions(&[
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
        ]);
        let mut faces = vec![
            FaceData::untagged(v[0], v[1], v[2]),       // good triangle
            FaceData::untagged(v[0], v[0], v[1]),       // degenerate: v0 == v0
        ];
        patch_small_boundary_holes(&mut faces, &pool);
        // Degenerate face should have been removed in Step 1
        for f in &faces {
            assert!(
                f.vertices[0] != f.vertices[1]
                    && f.vertices[1] != f.vertices[2]
                    && f.vertices[0] != f.vertices[2],
                "no degenerate faces should remain"
            );
        }
    }

    // ── Duplicate face removal ───────────────────────────────────────────

    /// Step 2 of patch deduplicates faces with the same vertex set.
    #[test]
    fn patch_removes_duplicate_faces() {
        let (pool, v) = pool_with_positions(&[
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.5, 0.5, 1.0),
        ]);
        // Closed tetrahedron with one face duplicated
        let mut faces = vec![
            FaceData::untagged(v[0], v[2], v[1]),
            FaceData::untagged(v[0], v[1], v[3]),
            FaceData::untagged(v[1], v[2], v[3]),
            FaceData::untagged(v[2], v[0], v[3]),
            FaceData::untagged(v[0], v[1], v[3]), // duplicate of face 1
        ];
        patch_small_boundary_holes(&mut faces, &pool);
        // Should have exactly 4 unique faces
        assert_eq!(faces.len(), 4, "duplicate face must be removed");
    }

    // ── Non-manifold edge repair ─────────────────────────────────────────

    /// Step 3 of patch resolves non-manifold edges by keeping the larger-area
    /// face.
    #[test]
    fn patch_resolves_non_manifold_by_keeping_larger() {
        let (pool, v) = pool_with_positions(&[
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.5, 1.0, 0.0),  // large triangle apex
            p(0.5, 0.01, 0.0), // sliver triangle apex (near the base)
        ]);
        // Two faces share directed edge v0->v1: big triangle and tiny sliver
        let mut faces = vec![
            FaceData::untagged(v[0], v[1], v[2]), // big area
            FaceData::untagged(v[0], v[1], v[3]), // tiny area
        ];
        patch_small_boundary_holes(&mut faces, &pool);
        // The sliver should be removed (smaller area on same directed edge)
        assert_eq!(faces.len(), 1, "non-manifold repair must keep only one face per half-edge");
        // The surviving face should be the big one
        let surviving = &faces[0];
        assert!(
            surviving.vertices.contains(&v[2]),
            "bigger-area face must survive non-manifold resolution"
        );
    }

    // ── Boundary hole patching (triangle hole) ───────────────────────────

    /// A box missing one quad face (= 2 triangles) creates a rectangular
    /// boundary loop. The patcher should fill it.
    #[test]
    fn patch_fills_quad_hole_in_box() {
        // Build an axis-aligned unit cube with one face (2 triangles) missing.
        //
        //   v4──v5
        //   │    │   top (z=1)
        //   v7──v6
        //
        //   v0──v1
        //   │    │   bottom (z=0)
        //   v3──v2
        let pts = [
            p(0.0, 0.0, 0.0), // 0
            p(1.0, 0.0, 0.0), // 1
            p(1.0, 1.0, 0.0), // 2
            p(0.0, 1.0, 0.0), // 3
            p(0.0, 0.0, 1.0), // 4
            p(1.0, 0.0, 1.0), // 5
            p(1.0, 1.0, 1.0), // 6
            p(0.0, 1.0, 1.0), // 7
        ];
        let (pool, v) = pool_with_positions(&pts);

        // 5 faces × 2 triangles = 10 triangles; omit the +Y face (v2,v3,v7,v6)
        let mut faces = vec![
            // bottom (z=0, normal -Z)
            FaceData::untagged(v[0], v[2], v[1]),
            FaceData::untagged(v[0], v[3], v[2]),
            // top (z=1, normal +Z)
            FaceData::untagged(v[4], v[5], v[6]),
            FaceData::untagged(v[4], v[6], v[7]),
            // front (y=0, normal -Y)
            FaceData::untagged(v[0], v[1], v[5]),
            FaceData::untagged(v[0], v[5], v[4]),
            // back: OMITTED — this is the hole
            // left (x=0, normal -X)
            FaceData::untagged(v[0], v[4], v[7]),
            FaceData::untagged(v[0], v[7], v[3]),
            // right (x=1, normal +X)
            FaceData::untagged(v[1], v[2], v[6]),
            FaceData::untagged(v[1], v[6], v[5]),
        ];

        let boundary_before = boundary_half_edges(&faces);
        assert!(
            !boundary_before.is_empty(),
            "box with missing face must have boundary edges"
        );

        patch_small_boundary_holes(&mut faces, &pool);

        let boundary_after = boundary_half_edges(&faces);
        assert!(
            boundary_after.is_empty(),
            "patcher must close the quad hole — {} boundary edges remain",
            boundary_after.len()
        );
        // Should have original 10 + 2 patch triangles = 12
        assert!(
            faces.len() >= 12,
            "patched box must have at least 12 faces, got {}",
            faces.len()
        );
    }

    // ── Sliver face tolerance ────────────────────────────────────────────

    /// Face with area ratio < 1e-8 * max_edge should be removed by Step 1.
    #[test]
    fn patch_removes_extreme_sliver() {
        let (pool, v) = pool_with_positions(&[
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.5, 1e-10, 0.0), // extreme sliver
            p(0.5, 0.5, 0.0),   // normal vertex for a good triangle
        ]);
        let mut faces = vec![
            FaceData::untagged(v[0], v[1], v[2]), // sliver
            FaceData::untagged(v[0], v[1], v[3]), // good
        ];
        patch_small_boundary_holes(&mut faces, &pool);
        // Sliver should be gone
        assert_eq!(faces.len(), 1, "extreme sliver must be removed");
        assert!(
            faces[0].vertices.contains(&v[3]),
            "good triangle must survive"
        );
    }

    // ── Determinism ──────────────────────────────────────────────────────

    /// Running patch twice on the same input must produce the same output.
    #[test]
    fn patch_is_idempotent() {
        let (pool, v) = pool_with_positions(&[
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(1.0, 1.0, 0.0),
        ]);
        let make_faces = || {
            vec![
                FaceData::untagged(v[0], v[1], v[2]),
                FaceData::untagged(v[1], v[3], v[2]),
            ]
        };

        let mut faces1 = make_faces();
        patch_small_boundary_holes(&mut faces1, &pool);
        let snapshot1: Vec<[VertexId; 3]> = faces1.iter().map(|f| f.vertices).collect();

        let mut faces2 = make_faces();
        patch_small_boundary_holes(&mut faces2, &pool);
        let snapshot2: Vec<[VertexId; 3]> = faces2.iter().map(|f| f.vertices).collect();

        assert_eq!(snapshot1, snapshot2, "patch must be deterministic");
    }
}
