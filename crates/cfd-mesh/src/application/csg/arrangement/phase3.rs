//! Phase 3.5 post-corefine vertex consolidation.
//!
//! Performs cross-mesh-only near-duplicate unification after independent CDTs.

use std::collections::{HashMap, HashSet};

use super::classify::FragRecord;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Real;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Build cross-mesh consolidation map from near-duplicate A/B vertex pairs.
///
/// # Algorithm
///
/// 1. Hash the smaller partition (`A` or `B`) into a uniform 3-D grid.
/// 2. For each probe vertex in the larger partition, search 27 neighboring
///    cells and collect near-duplicate A/B pairs within `tol_sq`.
/// 3. Union these pairs in a dense indexed DSU; unions are directed `B -> A`
///    to preserve canonical A-side roots.
/// 4. Emit `vid -> root_vid` substitutions for non-root vertices.
///
/// # Theorem — Indexed-DSU Equivalence
///
/// Replacing a key-addressed DSU (`HashMap<VertexId, VertexId>`) with an
/// index-addressed DSU over a bijection `VertexId <-> index` preserves the
/// partition relation exactly because `find`/`union` are performed over the
/// same pair sequence under an isomorphic state encoding. ∎
///
/// # Theorem — A-Root Invariant
///
/// Every accepted union is applied as `parent[root_b] = root_a`, where `a` is
/// drawn from partition `A` and `b` from partition `B`. Therefore every mixed
/// connected component has at least one A vertex and, by induction over unions,
/// its representative remains an A-side vertex. ∎
fn build_cross_mesh_merge_map(
    pure_a: &[VertexId],
    pure_b: &[VertexId],
    pool: &VertexPool,
    tol_sq: Real,
) -> HashMap<VertexId, VertexId> {
    if pure_a.is_empty() || pure_b.is_empty() {
        return HashMap::new();
    }

    let tol = tol_sq.sqrt();
    if tol <= 0.0 {
        return HashMap::new();
    }
    let inv_cell = 1.0 / tol;
    let grid_on_a = pure_a.len() <= pure_b.len();
    let (grid_source, probe_source): (&[VertexId], &[VertexId]) = if grid_on_a {
        (pure_a, pure_b)
    } else {
        (pure_b, pure_a)
    };

    let mut grid: HashMap<(i64, i64, i64), Vec<VertexId>> =
        HashMap::with_capacity(grid_source.len());
    for &vid in grid_source {
        let p = pool.position(vid);
        let ix = (p.x * inv_cell).floor() as i64;
        let iy = (p.y * inv_cell).floor() as i64;
        let iz = (p.z * inv_cell).floor() as i64;
        grid.entry((ix, iy, iz)).or_default().push(vid);
    }

    let all_merge_vids: Vec<VertexId> = pure_a.iter().chain(pure_b.iter()).copied().collect();
    let mut index_of: HashMap<VertexId, usize> = HashMap::with_capacity(all_merge_vids.len());
    for (i, &vid) in all_merge_vids.iter().enumerate() {
        index_of.insert(vid, i);
    }

    let mut parent: Vec<usize> = (0..all_merge_vids.len()).collect();

    fn find_root(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            let p = parent[x];
            let gp = parent[p];
            parent[x] = gp;
            x = gp;
        }
        x
    }

    for &probe in probe_source {
        let probe_p = pool.position(probe);
        let ix = (probe_p.x * inv_cell).floor() as i64;
        let iy = (probe_p.y * inv_cell).floor() as i64;
        let iz = (probe_p.z * inv_cell).floor() as i64;
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    if let Some(cands) = grid.get(&(ix + dx, iy + dy, iz + dz)) {
                        for &cand in cands {
                            let (va, vb, pa, pb) = if grid_on_a {
                                // cand from A grid, probe from B set
                                (cand, probe, pool.position(cand), probe_p)
                            } else {
                                // probe from A set, cand from B grid
                                (probe, cand, probe_p, pool.position(cand))
                            };
                            if (pb - pa).norm_squared() < tol_sq {
                                let ia = index_of[&va];
                                let ib = index_of[&vb];
                                let ra = find_root(&mut parent, ia);
                                let rb = find_root(&mut parent, ib);
                                if ra != rb {
                                    // A-root invariant: attach B-root under A-root.
                                    parent[rb] = ra;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
    for &vid in &all_merge_vids {
        let i = index_of[&vid];
        let root_i = find_root(&mut parent, i);
        let root = all_merge_vids[root_i];
        if root != vid {
            merge_map.insert(vid, root);
        }
    }
    merge_map
}

/// Consolidate near-duplicate corefined vertices across meshes A and B.
///
/// Uses a spatial hash and union-find over vertices that are exclusive to A or
/// exclusive to B. Vertices shared by both sources are preserved as seam anchors.
pub(crate) fn consolidate_cross_mesh_vertices(frags: &mut Vec<FragRecord>, pool: &VertexPool) {
    // Tolerance: 2x the corefine Steiner snap tolerance (WELD_TOL = 1e-4).
    const CONSOLIDATE_TOL: Real = 2e-4;
    const CONSOLIDATE_TOL_SQ: Real = CONSOLIDATE_TOL * CONSOLIDATE_TOL;

    let mut vids_a: HashSet<VertexId> = HashSet::new();
    let mut vids_b: HashSet<VertexId> = HashSet::new();
    for fr in &*frags {
        for &v in &fr.face.vertices {
            if fr.from_a {
                vids_a.insert(v);
            } else {
                vids_b.insert(v);
            }
        }
    }

    // Already shared seam IDs must never move.
    let seam_vids: HashSet<VertexId> = vids_a.intersection(&vids_b).copied().collect();
    let pure_a: Vec<VertexId> = vids_a.difference(&seam_vids).copied().collect();
    let pure_b: Vec<VertexId> = vids_b.difference(&seam_vids).copied().collect();

    if pure_a.is_empty() || pure_b.is_empty() {
        return;
    }
    let merge_map = build_cross_mesh_merge_map(&pure_a, &pure_b, pool, CONSOLIDATE_TOL_SQ);
    if merge_map.is_empty() {
        return;
    }

    for fr in frags.iter_mut() {
        for v in &mut fr.face.vertices {
            if let Some(&root) = merge_map.get(v) {
                *v = root;
            }
        }
    }
    frags.retain(|fr| {
        let v = fr.face.vertices;
        v[0] != v[1] && v[1] != v[2] && v[2] != v[0]
    });

    let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(frags.len());
    frags.retain(|fr| {
        let mut key = fr.face.vertices;
        key.sort();
        seen.insert(key)
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};

    #[test]
    fn merge_map_prefers_a_side_root() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        // Keep points outside weld tolerance (1e-4) but within consolidate tol (2e-4).
        let b0 = pool.insert_or_weld(Point3r::new(1.5e-4, 0.0, 0.0), n);

        let map = build_cross_mesh_merge_map(&[a0], &[b0], &pool, 4.0e-8);
        assert_eq!(
            map.get(&b0),
            Some(&a0),
            "B-side near-duplicate should map to A-side root"
        );
        assert!(
            !map.contains_key(&a0),
            "A-side canonical root should not be remapped"
        );
    }

    #[test]
    fn adversarial_bridge_keeps_a_root_in_component() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);

        // b0 is close to both a0 and a1, forcing transitive unions.
        let a0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let a1 = pool.insert_or_weld(Point3r::new(3.0e-4, 0.0, 0.0), n);
        let b0 = pool.insert_or_weld(Point3r::new(1.5e-4, 0.0, 0.0), n);

        let map = build_cross_mesh_merge_map(&[a0, a1], &[b0], &pool, 4.0e-8);

        let root_of_b = map.get(&b0).copied().unwrap_or(b0);
        assert!(
            root_of_b == a0 || root_of_b == a1,
            "component representative must remain on A-side"
        );
    }
}
