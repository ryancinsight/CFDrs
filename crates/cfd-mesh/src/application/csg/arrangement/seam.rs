//! Seam stitching passes for arrangement output repair.
//!
//! Keeps exact/constrained fixes first and only falls back to bounded
//! tolerance-based vertex merges when topology-preserving passes stall.

use hashbrown::HashMap;

use super::mesh_ops::{apply_vertex_merge, boundary_half_edges, merge_root};
use super::snap_round;
#[cfg(test)]
use crate::application::csg::diagnostics::trace_enabled;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Real;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

fn cell_key(p: &nalgebra::Point3<Real>, inv_cell: Real) -> (i64, i64, i64) {
    (
        (p.x * inv_cell + 0.5).floor() as i64,
        (p.y * inv_cell + 0.5).floor() as i64,
        (p.z * inv_cell + 0.5).floor() as i64,
    )
}

/// Build merge pairs by mutual-nearest-neighbor (MNN) matching within
/// `max_dist_sq`, accelerated by a 3-D uniform grid.
///
/// # Algorithm
///
/// 1. Hash every boundary vertex into a cubic cell of edge `h = sqrt(max_dist_sq)`.
/// 2. For each vertex `v`, search only its 27 neighboring cells for the nearest
///    neighbor `nn(v)` within distance bound.
/// 3. Accept pair `(u,v)` only if `nn(u)=v` and `nn(v)=u` (mutual nearest).
///
/// # Theorem — Symmetric Pairing Safety
///
/// MNN yields a matching (vertex degree ≤ 1) because each accepted edge is
/// bidirectional nearest-neighbor; therefore one vertex cannot be paired with
/// two distinct vertices in the same pass. This prevents fan-collapse artifacts
/// common in one-sided greedy nearest merges at V-branch seams. ∎
fn build_mutual_nearest_merge_map(
    bnd_verts: &[VertexId],
    max_dist_sq: Real,
    pool: &VertexPool,
) -> HashMap<VertexId, VertexId> {
    if bnd_verts.len() < 2 || max_dist_sq <= 0.0 {
        return HashMap::new();
    }
    let cell = max_dist_sq.sqrt().max(1e-12);
    let inv_cell = 1.0 / cell;

    let mut grid: HashMap<(i64, i64, i64), Vec<VertexId>> = HashMap::new();
    for &vid in bnd_verts {
        grid.entry(cell_key(pool.position(vid), inv_cell))
            .or_default()
            .push(vid);
    }

    let mut nearest: HashMap<VertexId, VertexId> = HashMap::new();
    for &vi in bnd_verts {
        let pi = pool.position(vi);
        let (cx, cy, cz) = cell_key(pi, inv_cell);
        let mut best_d = max_dist_sq;
        let mut best_v: Option<VertexId> = None;
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    let key = (cx + dx, cy + dy, cz + dz);
                    let Some(cands) = grid.get(&key) else {
                        continue;
                    };
                    for &vj in cands {
                        if vj == vi {
                            continue;
                        }
                        let d = (pool.position(vj) - pi).norm_squared();
                        if d < best_d {
                            best_d = d;
                            best_v = Some(vj);
                        } else if d == best_d && best_v.is_none_or(|best| vj.raw() < best.raw()) {
                            best_v = Some(vj);
                        }
                    }
                }
            }
        }
        if let Some(vj) = best_v {
            nearest.insert(vi, vj);
        }
    }

    let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
    for (&u, &v) in &nearest {
        if u.raw() >= v.raw() {
            continue;
        }
        if nearest.get(&v).is_some_and(|&vv| vv == u) {
            let (keep, discard) = if u < v { (u, v) } else { (v, u) };
            merge_map.insert(discard, keep);
        }
    }
    merge_map
}

/// Grid-accelerated greedy nearest merge map used as last-resort fallback when
/// MNN finds no pairs.
///
/// # Algorithm
///
/// 1. Hash boundary vertices into cubic cells of edge
///    `h = sqrt(max_dist_sq)`.
/// 2. For each unmerged vertex `vi` (in sorted ID order), search only the 27
///    neighboring cells for candidate `vj` with `vj > vi` and unmerged.
/// 3. Pick the nearest such `vj` within threshold and merge `vj -> vi`.
///
/// # Theorem — 27-cell candidate completeness
///
/// Any candidate with distance `< sqrt(max_dist_sq)` must lie in one of the
/// 27 neighboring cells when cell size equals that threshold. Therefore the
/// local-cell query returns exactly the same candidate set as a global scan
/// under the same distance bound. ∎
fn build_greedy_nearest_merge_map(
    bnd_verts: &[VertexId],
    max_dist_sq: Real,
    pool: &VertexPool,
) -> HashMap<VertexId, VertexId> {
    if bnd_verts.len() < 2 || max_dist_sq <= 0.0 {
        return HashMap::new();
    }

    let cell = max_dist_sq.sqrt().max(1e-12);
    let inv_cell = 1.0 / cell;

    let mut grid: HashMap<(i64, i64, i64), Vec<VertexId>> = HashMap::new();
    for &vid in bnd_verts {
        grid.entry(cell_key(pool.position(vid), inv_cell))
            .or_default()
            .push(vid);
    }
    for verts in grid.values_mut() {
        verts.sort_unstable();
    }

    let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
    for &vi in bnd_verts {
        if merge_map.contains_key(&vi) {
            continue;
        }

        let pi = pool.position(vi);
        let mut best_d = max_dist_sq;
        let mut best_j: Option<VertexId> = None;

        let (cx, cy, cz) = cell_key(pi, inv_cell);
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    let key = (cx + dx, cy + dy, cz + dz);
                    let Some(cands) = grid.get(&key) else {
                        continue;
                    };
                    for &vj in cands {
                        if vj.raw() <= vi.raw() || merge_map.contains_key(&vj) {
                            continue;
                        }
                        let d = (pool.position(vj) - pi).norm_squared();
                        if d < best_d
                            || (d == best_d && best_j.is_some_and(|best| vj.raw() < best.raw()))
                        {
                            best_d = d;
                            best_j = Some(vj);
                        }
                    }
                }
            }
        }

        if let Some(vj) = best_j {
            // vj > vi by construction, matching legacy merge direction.
            merge_map.insert(vj, vi);
        }
    }
    merge_map
}

// Phase 5.5 helper: stitch boundary seams from unresolved intersection gaps.
//
// When two surfaces meet at a shallow angle, the CDT co-refinement may produce
// zigzag seam boundaries. This helper now runs exact/constrained repair first
// (T-junction split + CDT loop fill), then falls back to bounded tolerance
// merges only when exact passes make no progress.
pub(crate) fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    // Exact-first prepass before any tolerance merge.
    if !boundary_half_edges(faces).is_empty() {
        snap_round::snap_round_tjunctions(faces, pool);
    }

    // === Pass 1: iterative short-edge collapse (fallback) ===
    for iter_idx in 0..6_usize {
        #[cfg(not(test))]
        let _ = iter_idx;
        let mut boundary_edges = boundary_half_edges(faces);
        if boundary_edges.is_empty() {
            return;
        }

        // Exact/constrained operations first for this iteration.
        let before_exact = boundary_edges.len();
        snap_round::snap_round_tjunctions(faces, pool);
        boundary_edges = boundary_half_edges(faces);
        if boundary_edges.is_empty() {
            return;
        }
        if boundary_edges.len() < before_exact {
            continue;
        }

        // Compute edge lengths.
        let mut edge_info: Vec<(Real, VertexId, VertexId)> = boundary_edges
            .iter()
            .map(|&(vi, vj)| {
                let d = (pool.position(vj) - pool.position(vi)).norm_squared();
                (d, vi, vj)
            })
            .collect();
        edge_info.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let min_len_sq = edge_info.first().map_or(0.0, |e| e.0);
        let max_len_sq = edge_info.last().map_or(0.0, |e| e.0);

        // Need at least 2x spread to identify a bimodal distribution.
        if min_len_sq <= 0.0 || max_len_sq < 2.0 * min_len_sq {
            break;
        }

        // Threshold: geometric mean of min and max.
        //
        // An additional hard cap here proved too aggressive for shallow-angle
        // elbow and branch seams, leaving open zipper gaps that the exact split
        // pass had already localized to a narrow bimodal edge-length band.
        let threshold_sq = (min_len_sq * max_len_sq).sqrt();

        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
        for &(len_sq, vi, vj) in &edge_info {
            if len_sq >= threshold_sq {
                break;
            }
            let root_i = merge_root(&merge_map, vi);
            let root_j = merge_root(&merge_map, vj);
            if root_i == root_j {
                continue;
            }
            let (keep, discard) = if root_i < root_j {
                (root_i, root_j)
            } else {
                (root_j, root_i)
            };
            merge_map.insert(discard, keep);
        }

        if merge_map.is_empty() {
            break;
        }

        #[cfg(test)]
        if trace_enabled() {
            eprintln!(
                "[stitch-p1 {}] {} bnd, {} short (< {:.6}), {} merges",
                iter_idx,
                boundary_edges.len(),
                edge_info.iter().filter(|e| e.0 < threshold_sq).count(),
                threshold_sq.sqrt(),
                merge_map.len(),
            );
        }

        apply_vertex_merge(faces, &merge_map, pool);
    }

    // === Pass 2: bounded nearest-boundary-vertex merge (last resort) ===
    // This pass is entered only when exact passes fail to reduce boundary edges.
    for iter_idx in 0..4_usize {
        #[cfg(not(test))]
        let _ = iter_idx;
        let mut boundary_edges = boundary_half_edges(faces);
        if boundary_edges.is_empty() {
            return;
        }

        let before_exact = boundary_edges.len();
        snap_round::snap_round_tjunctions(faces, pool);
        boundary_edges = boundary_half_edges(faces);
        if boundary_edges.is_empty() {
            return;
        }
        if boundary_edges.len() < before_exact {
            continue;
        }

        let mut bnd_verts: Vec<VertexId> =
            boundary_edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        bnd_verts.sort();
        bnd_verts.dedup();

        // Adaptive tolerance: 50% of average boundary edge length,
        // capped to prevent merging vertices across tube cross-sections.
        let avg_len_sq: Real = boundary_edges
            .iter()
            .map(|&(vi, vj)| (pool.position(vj) - pool.position(vi)).norm_squared())
            .sum::<Real>()
            / boundary_edges.len().max(1) as Real;
        let wide_tol_sq = (avg_len_sq * 0.25).min(0.01); // (0.5 * avg_len)^2, capped at 0.1

        let mut merge_map = build_mutual_nearest_merge_map(&bnd_verts, wide_tol_sq, pool);
        if merge_map.is_empty() {
            merge_map = build_greedy_nearest_merge_map(&bnd_verts, wide_tol_sq, pool);
        }

        if merge_map.is_empty() {
            break;
        }

        #[cfg(test)]
        if trace_enabled() {
            eprintln!(
                "[stitch-p2 {}] {} bnd edges, {} bnd verts, tol={:.6}, {} merges",
                iter_idx,
                boundary_edges.len(),
                bnd_verts.len(),
                wide_tol_sq.sqrt(),
                merge_map.len(),
            );
        }

        apply_vertex_merge(faces, &merge_map, pool);
    }
}

/// Conservative boundary seam stitch: pass-1-only with hard capped threshold.
///
/// Used after patch cleanup exposes boundary edges to avoid collapsing broader
/// geometric features while still closing narrow residual seam gaps.
pub(crate) fn stitch_boundary_seams_conservative(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    // Hard cap: never merge edges longer than 0.02 units.
    // This closes residual seam gaps exposed by patch cleanup without
    // collapsing macroscopic geometry across a tube cross-section.
    const MAX_THRESHOLD_SQ: Real = 4e-4; // 0.02^2

    for iter_idx in 0..4_usize {
        #[cfg(not(test))]
        let _ = iter_idx;
        let mut boundary_edges = boundary_half_edges(faces);
        if boundary_edges.is_empty() {
            return;
        }

        // Exact/constrained operations first; merge fallback only if no progress.
        let before_exact = boundary_edges.len();
        snap_round::snap_round_tjunctions(faces, pool);
        boundary_edges = boundary_half_edges(faces);
        if boundary_edges.is_empty() {
            return;
        }
        if boundary_edges.len() < before_exact {
            continue;
        }

        let mut edge_info: Vec<(Real, VertexId, VertexId)> = boundary_edges
            .iter()
            .map(|&(vi, vj)| {
                let d = (pool.position(vj) - pool.position(vi)).norm_squared();
                (d, vi, vj)
            })
            .collect();
        edge_info.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let min_len_sq = edge_info.first().map_or(0.0, |e| e.0);
        let max_len_sq = edge_info.last().map_or(0.0, |e| e.0);

        if min_len_sq <= 0.0 || max_len_sq < 2.0 * min_len_sq {
            break;
        }

        // Capped threshold: geometric mean but never above MAX_THRESHOLD_SQ.
        let threshold_sq = (min_len_sq * max_len_sq).sqrt().min(MAX_THRESHOLD_SQ);

        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
        for &(len_sq, vi, vj) in &edge_info {
            if len_sq >= threshold_sq {
                break;
            }
            let root_i = merge_root(&merge_map, vi);
            let root_j = merge_root(&merge_map, vj);
            if root_i == root_j {
                continue;
            }
            let (keep, discard) = if root_i < root_j {
                (root_i, root_j)
            } else {
                (root_j, root_i)
            };
            merge_map.insert(discard, keep);
        }

        if merge_map.is_empty() {
            break;
        }

        #[cfg(test)]
        if trace_enabled() {
            eprintln!(
                "[stitch-cons {}] {} bnd, {} short (< {:.6}), {} merges",
                iter_idx,
                boundary_edges.len(),
                edge_info.iter().filter(|e| e.0 < threshold_sq).count(),
                threshold_sq.sqrt(),
                merge_map.len(),
            );
        }

        apply_vertex_merge(faces, &merge_map, pool);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};
    use proptest::prelude::*;

    fn build_greedy_nearest_merge_map_bruteforce_reference(
        bnd_verts: &[VertexId],
        max_dist_sq: Real,
        pool: &VertexPool,
    ) -> HashMap<VertexId, VertexId> {
        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
        for (i, &vi) in bnd_verts.iter().enumerate() {
            if merge_map.contains_key(&vi) {
                continue;
            }
            let pi = pool.position(vi);
            let mut best_d = max_dist_sq;
            let mut best_j: Option<VertexId> = None;
            for &vj in bnd_verts.iter().skip(i + 1) {
                if merge_map.contains_key(&vj) {
                    continue;
                }
                let d = (pool.position(vj) - pi).norm_squared();
                if d < best_d {
                    best_d = d;
                    best_j = Some(vj);
                }
            }
            if let Some(vj) = best_j {
                merge_map.insert(vj, vi);
            }
        }
        merge_map
    }

    #[test]
    fn mnn_pairs_two_disjoint_close_pairs() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.000, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(0.001, 0.0, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(1.000, 0.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(1.001, 0.0, 0.0), n);
        let mut verts = vec![a, b, c, d];
        verts.sort();
        let map = build_mutual_nearest_merge_map(&verts, 0.01, &pool);
        assert_eq!(map.len(), 2, "two disjoint pairs should be matched");
        assert!(map.get(&b).is_some_and(|&k| k == a) || map.get(&a).is_some_and(|&k| k == b));
        assert!(map.get(&d).is_some_and(|&k| k == c) || map.get(&c).is_some_and(|&k| k == d));
    }

    #[test]
    fn adversarial_star_keeps_single_mutual_pair() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let center = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let l1 = pool.insert_or_weld(Point3r::new(0.001, 0.0, 0.0), n);
        let l2 = pool.insert_or_weld(Point3r::new(-0.001, 0.0, 0.0), n);
        let l3 = pool.insert_or_weld(Point3r::new(0.0, 0.001, 0.0), n);
        let mut verts = vec![center, l1, l2, l3];
        verts.sort();
        let map = build_mutual_nearest_merge_map(&verts, 0.01, &pool);
        assert_eq!(
            map.len(),
            1,
            "MNN should avoid collapsing a fan into the center in one pass"
        );
    }

    #[test]
    fn greedy_grid_matches_bruteforce_reference_small_case() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let mut verts = vec![
            pool.insert_or_weld(Point3r::new(0.000, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(0.003, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(0.007, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(0.100, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(0.102, 0.0, 0.0), n),
        ];
        verts.sort();
        verts.dedup();

        let max_dist_sq = 0.01 * 0.01;
        let fast = build_greedy_nearest_merge_map(&verts, max_dist_sq, &pool);
        let brute = build_greedy_nearest_merge_map_bruteforce_reference(&verts, max_dist_sq, &pool);
        assert_eq!(fast, brute);
    }

    proptest! {
        #[test]
        fn greedy_grid_matches_bruteforce_reference_property(
            coords in prop::collection::vec((-50_i16..50_i16, -50_i16..50_i16, -10_i16..10_i16), 4..28),
            max_step in 1_i16..20_i16
        ) {
            let mut pool = VertexPool::default_millifluidic();
            let n = Vector3r::new(0.0, 0.0, 1.0);
            let mut verts = Vec::new();

            for (x, y, z) in coords {
                let p = Point3r::new(f64::from(x) * 0.01, f64::from(y) * 0.01, f64::from(z) * 0.01);
                verts.push(pool.insert_or_weld(p, n));
            }

            verts.sort();
            verts.dedup();

            let max_d = f64::from(max_step) * 0.01;
            let max_dist_sq = max_d * max_d;

            let fast = build_greedy_nearest_merge_map(&verts, max_dist_sq, &pool);
            let brute = build_greedy_nearest_merge_map_bruteforce_reference(&verts, max_dist_sq, &pool);
            prop_assert_eq!(fast, brute);
        }
    }

    // ── Seam stitch tests ────────────────────────────────────────────────

    /// Stitching an empty face set must not panic.
    #[test]
    fn stitch_boundary_seams_empty_faces() {
        let pool = VertexPool::default_millifluidic();
        let mut faces: Vec<FaceData> = Vec::new();
        stitch_boundary_seams(&mut faces, &pool);
        assert!(faces.is_empty());
    }

    /// A closed tetrahedron should have no boundary seams; stitching is a no-op.
    #[test]
    fn stitch_boundary_seams_closed_mesh_noop() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_or_weld(Point3r::new(0.5, 1.0, 0.0), n);
        let v3 = pool.insert_or_weld(Point3r::new(0.5, 0.5, 1.0), n);
        let mut faces = vec![
            FaceData::untagged(v0, v2, v1),
            FaceData::untagged(v0, v1, v3),
            FaceData::untagged(v1, v2, v3),
            FaceData::untagged(v2, v0, v3),
        ];
        let before = faces.len();
        stitch_boundary_seams(&mut faces, &pool);
        assert_eq!(faces.len(), before, "closed mesh must not change under stitch");
    }

    /// Conservative stitch on empty faces must not panic.
    #[test]
    fn stitch_boundary_seams_conservative_empty() {
        let pool = VertexPool::default_millifluidic();
        let mut faces: Vec<FaceData> = Vec::new();
        stitch_boundary_seams_conservative(&mut faces, &pool);
        assert!(faces.is_empty());
    }

    /// Determinism: stitching the same input twice must produce the same output.
    #[test]
    fn stitch_boundary_seams_deterministic() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_or_weld(Point3r::new(0.5, 1.0, 0.0), n);

        let make_faces = || vec![FaceData::untagged(v0, v1, v2)];

        let mut faces1 = make_faces();
        stitch_boundary_seams(&mut faces1, &pool);

        let mut faces2 = make_faces();
        stitch_boundary_seams(&mut faces2, &pool);

        let s1: Vec<_> = faces1.iter().map(|f| f.vertices).collect();
        let s2: Vec<_> = faces2.iter().map(|f| f.vertices).collect();
        assert_eq!(s1, s2, "stitch must be deterministic");
    }
}
