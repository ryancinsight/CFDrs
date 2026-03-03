//! Phase 2c/2d coplanar group processing for arrangement CSG.
//!
//! Handles exact 2-D coplanar boolean resolution and seam-vertex injection
//! into adjacent non-coplanar barrel faces.

use std::collections::{HashMap, HashSet};

use super::super::boolean::BooleanOp;
use super::super::coplanar::basis::PlaneBasis;
use super::super::intersect::SnapSegment;
use super::dsu::DisjointSet;
use super::propagate::inject_cap_seam_into_barrels;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Point3r;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

pub(crate) struct CoplanarGroupIndex {
    pub(crate) plane_a: HashMap<usize, Vec<usize>>,
    pub(crate) plane_b: HashMap<usize, Vec<usize>>,
    pub(crate) rep_a: HashMap<usize, usize>,
}

pub(crate) struct CoplanarPhaseResult {
    pub(crate) a_used: HashSet<usize>,
    pub(crate) b_used: HashSet<usize>,
    pub(crate) results: HashMap<usize, Vec<FaceData>>,
}

/// Build coplanar connected components from narrow-phase coplanar face pairs.
///
/// # Algorithm
///
/// Build a bipartite graph over `A` and `B` face indices where each coplanar
/// pair `(a,b)` is an edge. Connected components are computed by disjoint-set
/// union (DSU). Each component becomes one coplanar group with:
/// - `plane_a[group] = {a_i}`
/// - `plane_b[group] = {b_j}`
///
/// # Theorem — Component Correctness
///
/// Let `G=(V,E)` be the bipartite graph from coplanar pairs. DSU unions exactly
/// the endpoints of each edge in `E`, so two vertices share a DSU root iff they
/// are connected by a path in `G`. Therefore DSU roots are exactly the connected
/// components of coplanar pair adjacency, which is precisely the grouping needed
/// for independent per-component 2-D coplanar Boolean processing. ∎
pub(crate) fn build_coplanar_group_index(
    num_faces_a: usize,
    num_faces_b: usize,
    coplanar_pairs: &[(usize, usize)],
) -> CoplanarGroupIndex {
    if coplanar_pairs.is_empty() {
        return CoplanarGroupIndex {
            plane_a: HashMap::new(),
            plane_b: HashMap::new(),
            rep_a: HashMap::new(),
        };
    }

    let mut pairs: Vec<(usize, usize)> = coplanar_pairs
        .iter()
        .copied()
        .filter(|&(a, b)| a < num_faces_a && b < num_faces_b)
        .collect();
    pairs.sort_unstable();
    pairs.dedup();

    if pairs.is_empty() {
        return CoplanarGroupIndex {
            plane_a: HashMap::new(),
            plane_b: HashMap::new(),
            rep_a: HashMap::new(),
        };
    }

    let total = num_faces_a + num_faces_b;
    let mut dsu = DisjointSet::new(total);

    for &(a, b) in &pairs {
        dsu.union(a, num_faces_a + b);
    }

    let mut root_to_group: HashMap<usize, usize> = HashMap::new();
    let mut plane_a: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut plane_b: HashMap<usize, Vec<usize>> = HashMap::new();

    for &(a, b) in &pairs {
        let root = dsu.find(a);
        let group = if let Some(&g) = root_to_group.get(&root) {
            g
        } else {
            let g = root_to_group.len();
            root_to_group.insert(root, g);
            g
        };
        plane_a.entry(group).or_default().push(a);
        plane_b.entry(group).or_default().push(b);
    }

    for v in plane_a.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    for v in plane_b.values_mut() {
        v.sort_unstable();
        v.dedup();
    }

    let mut rep_a: HashMap<usize, usize> = HashMap::with_capacity(plane_a.len());
    for (&group, a_faces) in &plane_a {
        if let Some(&a0) = a_faces.first() {
            rep_a.insert(group, a0);
        }
    }

    CoplanarGroupIndex {
        plane_a,
        plane_b,
        rep_a,
    }
}

/// Resolve all coplanar plane groups (Phase 2c) and propagate seam vertices
/// into barrel faces (Phase 2d).
pub(crate) fn process_coplanar_groups(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
    group_index: &CoplanarGroupIndex,
    segs_a: &mut HashMap<usize, Vec<SnapSegment>>,
    segs_b: &mut HashMap<usize, Vec<SnapSegment>>,
) -> CoplanarPhaseResult {
    let mut coplanar_a_used: HashSet<usize> = HashSet::new();
    let mut coplanar_b_used: HashSet<usize> = HashSet::new();
    let mut coplanar_results: HashMap<usize, Vec<FaceData>> = HashMap::new();

    for (key, a_idxs) in &group_index.plane_a {
        let Some(b_idxs) = group_index.plane_b.get(key) else {
            continue;
        };

        let cap_a: Vec<FaceData> = a_idxs.iter().map(|&i| faces_a[i]).collect();
        let cap_b: Vec<FaceData> = b_idxs.iter().map(|&i| faces_b[i]).collect();

        let basis = match basis_from_group(*key, &group_index.rep_a, faces_a, &cap_a, pool) {
            Some(b) => b,
            None => continue, // fallback to Phase 3/4
        };

        let coplanar_result =
            crate::application::csg::coplanar::boolean_coplanar(op, &cap_a, &cap_b, pool, &basis);

        // Collect original cap vertex IDs for exact seam-vertex detection.
        let mut original_vids: HashSet<VertexId> =
            HashSet::with_capacity((cap_a.len() + cap_b.len()) * 3);
        for f in cap_a.iter().chain(cap_b.iter()) {
            original_vids.insert(f.vertices[0]);
            original_vids.insert(f.vertices[1]);
            original_vids.insert(f.vertices[2]);
        }

        // New seam vertices are exactly those output IDs absent from originals.
        let mut seam_vids: HashSet<VertexId> = HashSet::with_capacity(coplanar_result.len() * 2);
        for face in &coplanar_result {
            for &vid in &face.vertices {
                if !original_vids.contains(&vid) {
                    seam_vids.insert(vid);
                }
            }
        }
        let mut seam_vids: Vec<VertexId> = seam_vids.into_iter().collect();
        seam_vids.sort_unstable();
        let seam_positions: Vec<Point3r> =
            seam_vids.iter().map(|&vid| *pool.position(vid)).collect();

        for &i in a_idxs {
            coplanar_a_used.insert(i);
        }
        for &i in b_idxs {
            coplanar_b_used.insert(i);
        }

        if !seam_positions.is_empty() {
            inject_cap_seam_into_barrels(
                faces_a,
                &coplanar_a_used,
                &basis.origin,
                &basis.normal,
                &seam_positions,
                segs_a,
                pool,
            );
            inject_cap_seam_into_barrels(
                faces_b,
                &coplanar_b_used,
                &basis.origin,
                &basis.normal,
                &seam_positions,
                segs_b,
                pool,
            );
        }

        coplanar_results.insert(*key, coplanar_result);
    }

    CoplanarPhaseResult {
        a_used: coplanar_a_used,
        b_used: coplanar_b_used,
        results: coplanar_results,
    }
}

/// Build a robust plane basis for the coplanar group.
///
/// Modern fast path: derive directly from the group representative triangle.
/// Fallback: full exact flat-plane detection over `cap_a`.
fn basis_from_group(
    group_id: usize,
    rep_a_by_group: &HashMap<usize, usize>,
    faces_a: &[FaceData],
    cap_a: &[FaceData],
    pool: &VertexPool,
) -> Option<PlaneBasis> {
    if let Some(&rep_idx) = rep_a_by_group.get(&group_id) {
        let rep_face = faces_a.get(rep_idx)?;
        let a = *pool.position(rep_face.vertices[0]);
        let b = *pool.position(rep_face.vertices[1]);
        let c = *pool.position(rep_face.vertices[2]);
        if let Some(basis) = PlaneBasis::from_triangle(&a, &b, &c) {
            return Some(basis);
        }
    }
    crate::application::csg::coplanar::detect_flat_plane(cap_a, pool)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn canonical_groups(index: &CoplanarGroupIndex) -> Vec<(Vec<usize>, Vec<usize>)> {
        let mut out: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
        let mut gids: Vec<usize> = index.plane_a.keys().copied().collect();
        gids.sort_unstable();
        for g in gids {
            let mut a = index.plane_a.get(&g).cloned().unwrap_or_default();
            let mut b = index.plane_b.get(&g).cloned().unwrap_or_default();
            a.sort_unstable();
            b.sort_unstable();
            out.push((a, b));
        }
        out.sort_unstable();
        out
    }

    #[test]
    fn adversarial_chain_bridge_merges_single_component() {
        // Chain across bipartite edges:
        // a0-b0-a1-b1-a2-b2
        let pairs = vec![(0, 0), (1, 0), (1, 1), (2, 1), (2, 2)];
        let idx = build_coplanar_group_index(4, 4, &pairs);
        let groups = canonical_groups(&idx);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0, vec![0, 1, 2]);
        assert_eq!(groups[0].1, vec![0, 1, 2]);
    }

    #[test]
    fn adversarial_permuted_duplicates_invariant() {
        let base = vec![(0, 0), (1, 0), (2, 1), (3, 2), (3, 3), (4, 3)];
        let mut noisy = vec![
            (3, 2),
            (0, 0),
            (3, 3),
            (1, 0),
            (3, 2),
            (4, 3),
            (2, 1),
            (0, 0),
        ];

        let a = build_coplanar_group_index(8, 8, &base);
        let b = build_coplanar_group_index(8, 8, &noisy);
        assert_eq!(canonical_groups(&a), canonical_groups(&b));

        noisy.reverse();
        let c = build_coplanar_group_index(8, 8, &noisy);
        assert_eq!(canonical_groups(&a), canonical_groups(&c));
    }
}
