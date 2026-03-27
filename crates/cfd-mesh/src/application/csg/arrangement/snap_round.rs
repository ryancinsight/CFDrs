//! T-junction snap-round repair for CSG boolean operations.
//!
//! Resolves T-junctions in the result face set by splitting face edges at
//! boundary vertices that lie on those edges.  Uses exact orientation
//! predicates with tolerance-based fallback for numerical robustness.
//!
//! ## Algorithm
//!
//! A T-junction occurs when boundary vertex V lies on an edge [A,B] of face
//! [A,B,C] but V is not part of that face.  Splitting [A,B,C] into [A,V,C]
//! and [V,B,C] inserts V into the mesh topology, pairing the boundary
//! half-edge.
//!
//! ## Key properties
//!
//! - Only ADDS faces (splits), never removes them
//! - Cannot create degenerate faces (V is verified to be interior to the edge)
//! - Cannot create topology cascades
//!
//! ## References
//!
//! - Attene (2010) "A lightweight approach to repairing digitized polygon meshes"
//! - Shewchuk (1996) exact orientation predicates for on-segment detection

use hashbrown::{HashMap, HashSet};

use super::mesh_ops::{boundary_half_edges, dedup_faces_unordered};
#[cfg(test)]
use crate::application::csg::diagnostics::trace_enabled;
use crate::application::csg::predicates3d::point_on_segment_exact;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Real;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Squared perpendicular distance tolerance for snap-round.
/// A boundary vertex V is considered to lie on face edge [A,B] if:
///   distance(V, line(A,B)) < sqrt(SNAP_TOL_SQ) * |AB|
/// i.e., sin(angle(AV, AB)) < sqrt(SNAP_TOL_SQ) ≈ 0.001
const SNAP_TOL_SQ: Real = 1e-6;
/// Endpoint margin for edge-parameter test `t`.
///
/// Candidate split vertices must be strictly interior to the edge and at least
/// this fraction away from either endpoint.
const SNAP_EDGE_PARAM_EPS: Real = 5e-3;

type FaceEdgeRef = (usize, VertexId, VertexId);

type SplitCandidate = (bool, Real, Real, VertexId, VertexId, VertexId);

// ── Helper functions ─────────────────────────────────────────────────────────

/// Build undirected boundary vertex adjacency from directed boundary edges.
fn build_boundary_adjacency(boundary: &[(VertexId, VertexId)]) -> HashMap<VertexId, Vec<VertexId>> {
    let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
    for &(a, b) in boundary {
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }
    for neigh in adj.values_mut() {
        neigh.sort();
        neigh.dedup();
    }
    adj
}

/// Build an endpoint-to-face-edge index for all current triangle edges.
///
/// # Algorithm
///
/// For each face edge `(a,b)` in face `fi`, insert `(fi,a,b)` into both
/// `index[a]` and `index[b]`.
///
/// # Theorem — Endpoint Index Completeness for Constrained T-Junction Search
///
/// Let `N(v)` be boundary-neighbors of boundary vertex `v`, and let guard
/// `endpoint_constrained(v,a,b)` be true iff `a∈N(v) or b∈N(v)`.
///
/// Any edge `(a,b)` that can pass `endpoint_constrained` must be incident to
/// at least one vertex in `N(v)`, hence appears in
/// `⋃_{u∈N(v)} index[u]`. Therefore searching only this union is equivalent to
/// full-face scanning under the same guard, but with much lower work.
fn build_endpoint_edge_index(faces: &[FaceData]) -> HashMap<VertexId, Vec<FaceEdgeRef>> {
    let mut index: HashMap<VertexId, Vec<FaceEdgeRef>> =
        HashMap::with_capacity(faces.len().saturating_mul(2));
    for (fi, face) in faces.iter().enumerate() {
        for edge_idx in 0..3_usize {
            let a = face.vertices[edge_idx];
            let b = face.vertices[(edge_idx + 1) % 3];
            index.entry(a).or_default().push((fi, a, b));
            index.entry(b).or_default().push((fi, a, b));
        }
    }
    index
}

/// Gather candidate face-edges for boundary vertex `v` using endpoint index.
///
/// Returned tuples are `(face_index, edge_a, edge_b)` with de-duplication by
/// canonical edge per-face key.
fn candidate_face_edges_for_vertex(
    v: VertexId,
    bnd_adj: &HashMap<VertexId, Vec<VertexId>>,
    endpoint_index: &HashMap<VertexId, Vec<FaceEdgeRef>>,
) -> Vec<FaceEdgeRef> {
    let Some(neighbors) = bnd_adj.get(&v) else {
        return Vec::new();
    };
    let mut seen: HashSet<(usize, VertexId, VertexId)> = HashSet::new();
    let mut out: Vec<FaceEdgeRef> = Vec::new();
    for &u in neighbors {
        let Some(edges) = endpoint_index.get(&u) else {
            continue;
        };
        for &(fi, a, b) in edges {
            let (mn, mx) = if a < b { (a, b) } else { (b, a) };
            if seen.insert((fi, mn, mx)) {
                out.push((fi, a, b));
            }
        }
    }
    out
}

#[inline]
fn endpoint_constrained(
    v: VertexId,
    a: VertexId,
    b: VertexId,
    bnd_adj: &HashMap<VertexId, Vec<VertexId>>,
) -> bool {
    bnd_adj
        .get(&v)
        .is_some_and(|neigh| neigh.contains(&a) || neigh.contains(&b))
}

fn update_best_split(
    best_split_for_face: &mut HashMap<usize, SplitCandidate>,
    face_index: usize,
    candidate: SplitCandidate,
) {
    match best_split_for_face.get_mut(&face_index) {
        Some(best) => {
            let best_key = (!best.0, best.1, best.2, best.5.raw(), best.3.raw(), best.4.raw());
            let cand_key = (
                !candidate.0,
                candidate.1,
                candidate.2,
                candidate.5.raw(),
                candidate.3.raw(),
                candidate.4.raw(),
            );
            if cand_key < best_key {
                *best = candidate;
            }
        }
        None => {
            best_split_for_face.insert(face_index, candidate);
        }
    }
}

// ── Main entry point ─────────────────────────────────────────────────────────

/// Resolve T-junctions by splitting face edges at boundary vertices.
///
/// A T-junction occurs when boundary vertex V lies on an edge [A,B] of
/// a face [A,B,C] but V is not part of that face.  This creates an
/// unpaired half-edge at V.  Splitting [A,B,C] into [A,V,C] and [V,B,C]
/// inserts V into the mesh topology, pairing the boundary half-edge.
///
/// Candidate splits are constrained to vertices that are boundary-adjacent to
/// one of the edge endpoints.  On-segment detection uses exact orientation
/// predicates first, then falls back to the legacy distance gate.
///
/// Unlike vertex merging, edge splitting:
/// - Only ADDS faces (splits), never removes them
/// - Cannot create degenerate faces (V is verified to be interior to the edge)
/// - Cannot create topology cascades
///
/// Iterates up to `MAX_SNAP_ITERS` times until no more T-junctions are found
/// or boundary edges are eliminated.
pub(crate) fn snap_round_tjunctions(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    const MAX_SNAP_ITERS: usize = 8;

    for iter_idx in 0..MAX_SNAP_ITERS {
        #[cfg(not(test))]
        let _ = iter_idx;
        let boundary = boundary_half_edges(faces);
        if boundary.is_empty() {
            break;
        }

        // Collect unique boundary vertices.
        let mut bnd_verts: Vec<VertexId> = boundary.iter().flat_map(|&(a, b)| [a, b]).collect();
        bnd_verts.sort();
        bnd_verts.dedup();
        let bnd_set: HashSet<VertexId> = bnd_verts.iter().copied().collect();
        let bnd_adj = build_boundary_adjacency(&boundary);
        let endpoint_index = build_endpoint_edge_index(faces);
        // For each face, keep the best constrained split candidate:
        // exact-hit first, then shortest residual distance, then centrality.
        let mut best_split_for_face: HashMap<usize, SplitCandidate> = HashMap::new();

        for &v in &bnd_verts {
            let pv = pool.position(v);
            let mut found_exact_candidate = false;
            for (fi, a, b) in candidate_face_edges_for_vertex(v, &bnd_adj, &endpoint_index) {
                let Some(face) = faces.get(fi) else {
                    continue;
                };
                // Skip faces that already contain this vertex.
                if face.vertices.contains(&v) {
                    continue;
                }
                if !endpoint_constrained(v, a, b, &bnd_adj) {
                    continue;
                }

                let pa = pool.position(a);
                let pb = pool.position(b);
                let ab = pb - pa;
                let edge_len_sq = ab.norm_squared();
                if edge_len_sq < 1e-30 {
                    continue;
                }

                let mut exact_hit = false;
                let t;
                let dist_metric;
                if let Some(t_exact) = point_on_segment_exact(pa, pb, pv) {
                    if t_exact <= SNAP_EDGE_PARAM_EPS || t_exact >= 1.0 - SNAP_EDGE_PARAM_EPS {
                        continue;
                    }
                    exact_hit = true;
                    found_exact_candidate = true;
                    t = t_exact;
                    dist_metric = 0.0;
                } else {
                    let av = pv - pa;

                    // Projection parameter: t = dot(AV, AB) / |AB|²
                    t = av.dot(&ab) / edge_len_sq;
                    // V must be strictly interior to edge (not at endpoints).
                    if t <= SNAP_EDGE_PARAM_EPS || t >= 1.0 - SNAP_EDGE_PARAM_EPS {
                        continue;
                    }

                    // Perpendicular distance: |cross(AB, AV)|² / |AB|²
                    let cross = ab.cross(&av);
                    let dist_sq = cross.norm_squared() / edge_len_sq;
                    if dist_sq > SNAP_TOL_SQ * edge_len_sq {
                        continue;
                    }
                    dist_metric = dist_sq / edge_len_sq.max(1e-30);
                }

                // Found a T-junction: V lies on edge [A, B] of face fi.
                // Keep the historical seam-ribbon guard for tolerance-only
                // hits, but allow exact constrained hits through.
                if !exact_hit && bnd_set.contains(&a) && bnd_set.contains(&b) {
                    continue;
                }

                let center_bias = (t - 0.5).abs();
                let candidate = (exact_hit, dist_metric, center_bias, a, b, v);
                update_best_split(&mut best_split_for_face, fi, candidate);
            }

            if found_exact_candidate {
                continue;
            }

            for (fi, face) in faces.iter().enumerate() {
                if face.vertices.contains(&v) {
                    continue;
                }
                for edge_idx in 0..3_usize {
                    let a = face.vertices[edge_idx];
                    let b = face.vertices[(edge_idx + 1) % 3];
                    let pa = pool.position(a);
                    let pb = pool.position(b);
                    let Some(t_exact) = point_on_segment_exact(pa, pb, pv) else {
                        continue;
                    };
                    if t_exact <= SNAP_EDGE_PARAM_EPS || t_exact >= 1.0 - SNAP_EDGE_PARAM_EPS {
                        continue;
                    }

                    let center_bias = (t_exact - 0.5).abs();
                    let candidate = (true, 0.0, center_bias, a, b, v);
                    update_best_split(&mut best_split_for_face, fi, candidate);
                }
            }
        }

        let mut splits: Vec<(usize, VertexId, VertexId, VertexId)> = best_split_for_face
            .into_iter()
            .map(|(fi, (_, _, _, a, b, v))| (fi, a, b, v))
            .collect();
        if splits.is_empty() {
            break;
        }

        #[cfg(test)]
        if trace_enabled() {
            eprintln!(
                "[snap-round {}] {} bnd edges, {} bnd verts, {} splits",
                iter_idx,
                boundary.len(),
                bnd_verts.len(),
                splits.len(),
            );
        }

        // Apply splits: replace each face with two sub-faces.
        // Process in reverse order of face index to avoid invalidating indices.
        splits.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        for (fi, a, b, v) in &splits {
            let face = &faces[*fi];
            // Find the third vertex (not a or b).
            let c = face
                .vertices
                .iter()
                .copied()
                .find(|&vid| vid != *a && vid != *b)
                .expect("face must have 3 distinct vertices");
            let region = face.region;

            // Determine winding order: face has edge [a, b] at some position.
            // The original face is [.., a, b, ..] in some winding.
            // We need [a, v, c] and [v, b, c] with correct winding.
            let verts = face.vertices;
            let (f1, f2) = if (verts[0] == *a && verts[1] == *b)
                || (verts[1] == *a && verts[2] == *b)
                || (verts[2] == *a && verts[0] == *b)
            {
                // Edge a→b is in forward winding: face is [a, b, c]
                // Split into [a, v, c] and [v, b, c]
                (
                    FaceData::new(*a, *v, c, region),
                    FaceData::new(*v, *b, c, region),
                )
            } else {
                // Edge b→a is in forward winding: face is [b, a, c] (effectively)
                // This means the original has a after b in winding.
                // Split into [b, v, c] and [v, a, c]
                (
                    FaceData::new(*b, *v, c, region),
                    FaceData::new(*v, *a, c, region),
                )
            };

            // Replace the original face with f1, append f2.
            faces[*fi] = f1;
            faces.push(f2);
        }

        // Remove any degenerate faces created by numerical edge cases.
        faces.retain(|f| {
            f.vertices[0] != f.vertices[1]
                && f.vertices[1] != f.vertices[2]
                && f.vertices[2] != f.vertices[0]
        });
        dedup_faces_unordered(faces);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};

    #[test]
    fn snap_round_splits_exact_tjunction_with_endpoint_constraint() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), n);
        let m = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(2.0, 1.0, 0.0), n);

        let mut faces = vec![FaceData::untagged(a, b, c), FaceData::untagged(m, b, d)];
        snap_round_tjunctions(&mut faces, &pool);

        assert_eq!(
            faces.len(),
            3,
            "one constrained exact split should be applied"
        );
        let with_mc = faces
            .iter()
            .filter(|f| f.vertices.contains(&m) && f.vertices.contains(&c))
            .count();
        assert_eq!(
            with_mc, 2,
            "split should replace [a,b,c] with two triangles using split vertex m"
        );
    }

    #[test]
    fn snap_round_splits_exact_collinear_vertex_without_endpoint_constraint() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), n);
        let m = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(1.0, 1.0, 0.0), n);
        let e = pool.insert_or_weld(Point3r::new(2.0, 1.0, 0.0), n);

        // m lies exactly on [a,b] but is not boundary-adjacent to a or b.
        let mut faces = vec![FaceData::untagged(a, b, c), FaceData::untagged(m, d, e)];
        snap_round_tjunctions(&mut faces, &pool);

        assert_eq!(
            faces.len(),
            3,
            "exact on-edge boundary vertices should trigger a split even without endpoint adjacency"
        );
        assert!(
            faces.iter().filter(|face| face.vertices.contains(&m)).count() >= 2,
            "split faces should contain the exact T-junction vertex"
        );
    }

    #[test]
    fn adversarial_endpoint_index_matches_full_scan_candidates() {
        use super::super::mesh_ops::boundary_half_edges;

        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let ids: Vec<VertexId> = (0..7)
            .map(|i| pool.insert_or_weld(Point3r::new(Real::from(i), 0.0, 0.0), n))
            .collect();
        let [v0, v1, v2, v3, v4, v5, v6] = <[VertexId; 7]>::try_from(ids).expect("7 ids");

        let faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v0, v2, v3),
            FaceData::untagged(v0, v3, v4),
            FaceData::untagged(v0, v4, v1),
            FaceData::untagged(v2, v5, v4),
            FaceData::untagged(v4, v5, v6),
        ];

        let boundary = boundary_half_edges(&faces);
        let bnd_adj = build_boundary_adjacency(&boundary);
        let endpoint_index = build_endpoint_edge_index(&faces);

        let mut boundary_vertices: Vec<VertexId> =
            boundary.iter().flat_map(|&(a, b)| [a, b]).collect();
        boundary_vertices.sort();
        boundary_vertices.dedup();

        for v in boundary_vertices {
            let mut full_scan: HashSet<(usize, VertexId, VertexId)> = HashSet::new();
            for (fi, face) in faces.iter().enumerate() {
                for edge_idx in 0..3_usize {
                    let a = face.vertices[edge_idx];
                    let b = face.vertices[(edge_idx + 1) % 3];
                    if !endpoint_constrained(v, a, b, &bnd_adj) {
                        continue;
                    }
                    let (mn, mx) = if a < b { (a, b) } else { (b, a) };
                    full_scan.insert((fi, mn, mx));
                }
            }

            let indexed: HashSet<(usize, VertexId, VertexId)> =
                candidate_face_edges_for_vertex(v, &bnd_adj, &endpoint_index)
                    .into_iter()
                    .map(|(fi, a, b)| {
                        let (mn, mx) = if a < b { (a, b) } else { (b, a) };
                        (fi, mn, mx)
                    })
                    .collect();

            assert_eq!(
                indexed, full_scan,
                "endpoint-index candidate set must match full-scan set"
            );
        }
    }
}
