//! Winding orientation consistency.
//!
//! For a closed manifold, every shared edge should be traversed in opposite
//! directions by its two adjacent faces (consistent orientation).

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::FaceId;
use crate::infrastructure::storage::edge_store::EdgeStore;
use crate::infrastructure::storage::face_store::FaceStore;

/// Check if all faces have consistent winding orientation.
///
/// For each manifold edge shared by two faces, the edge should appear as
/// `(a, b)` in one face and `(b, a)` in the other.
pub fn check_orientation(face_store: &FaceStore, edge_store: &EdgeStore) -> MeshResult<()> {
    for edge in edge_store.iter() {
        if edge.faces.len() != 2 {
            continue; // Skip non-manifold / boundary edges
        }

        let f0 = face_store.get(edge.faces[0]);
        let f1 = face_store.get(edge.faces[1]);
        let (ea, eb) = edge.vertices;

        // Find the directed edge in each face
        let dir0 = directed_edge_order(f0.vertices, ea, eb);
        let dir1 = directed_edge_order(f1.vertices, ea, eb);

        // For consistent orientation, the edge must be traversed in
        // opposite directions: one face has (a→b), the other (b→a).
        if let (Some(d0), Some(d1)) = (dir0, dir1) {
            if d0 == d1 {
                return Err(MeshError::InconsistentWinding {
                    face: edge.faces[1],
                });
            }
        }
    }
    Ok(())
}

/// Determine if edge `(a, b)` appears as `a→b` (true) or `b→a` (false) in the face.
fn directed_edge_order(
    verts: [crate::domain::core::index::VertexId; 3],
    a: crate::domain::core::index::VertexId,
    b: crate::domain::core::index::VertexId,
) -> Option<bool> {
    for i in 0..3 {
        let j = (i + 1) % 3;
        if verts[i] == a && verts[j] == b {
            return Some(true);
        }
        if verts[i] == b && verts[j] == a {
            return Some(false);
        }
    }
    None
}

/// Attempt to fix inconsistent winding by flipping faces.
///
/// Uses dual-graph BFS per connected component, flipping faces that have
/// inconsistent orientation relative to their neighbors.
///
/// # Theorem — component-complete consistency repair
///
/// Let `G` be the face-adjacency graph whose vertices are faces and whose edges
/// connect pairs of faces sharing a mesh edge. For each connected component of
/// `G`, BFS propagation from any seed face assigns an orientation parity label
/// to all reachable faces; flipping exactly those neighbors that violate the
/// local opposite-edge-direction constraint yields a consistent orientation on
/// that component. Repeating this for every unvisited seed covers all
/// components, so no disconnected component is skipped. ∎
pub fn fix_orientation(face_store: &mut FaceStore, edge_store: &EdgeStore) -> usize {
    let n_faces = face_store.len();
    if n_faces == 0 {
        return 0;
    }

    // Dense bitset lowers memory and branch overhead versus hash-set
    // membership in hot repair paths.
    let mut visited = vec![false; n_faces];
    let mut flipped = 0usize;
    let mut queue = std::collections::VecDeque::new();

    // Traverse every disconnected component. Seeding only face 0 leaves
    // orientation defects untouched in later components.
    for seed_idx in 0..n_faces {
        if visited[seed_idx] {
            continue;
        }
        visited[seed_idx] = true;
        let seed = FaceId::from_usize(seed_idx);
        queue.push_back(seed);

        while let Some(current) = queue.pop_front() {
            let current_face = *face_store.get(current);
            let current_idx = current.as_usize();

            // For each edge of the current face, check the neighbor
            for (ea, eb) in current_face.edges_canonical() {
                if let Some(eid) = edge_store.find_edge(ea, eb) {
                    let edge = edge_store.get(eid);
                    for &neighbor_fid in &edge.faces {
                        let neighbor_idx = neighbor_fid.as_usize();
                        if neighbor_idx == current_idx || visited[neighbor_idx] {
                            continue;
                        }
                        visited[neighbor_idx] = true;

                        let dir_current = directed_edge_order(current_face.vertices, ea, eb);
                        let neighbor_face = face_store.get(neighbor_fid);
                        let dir_neighbor = directed_edge_order(neighbor_face.vertices, ea, eb);

                        if let (Some(dc), Some(dn)) = (dir_current, dir_neighbor) {
                            if dc == dn {
                                // Same direction -> flip neighbor.
                                face_store.get_mut(neighbor_fid).flip();
                                flipped += 1;
                            }
                        }

                        queue.push_back(neighbor_fid);
                    }
                }
            }
        }
    }

    flipped
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::index::{RegionId, VertexId};
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use crate::infrastructure::storage::face_store::{FaceData, FaceStore};
    use proptest::prelude::*;

    fn tetra_faces(base: u32) -> [FaceData; 4] {
        let v0 = VertexId::new(base);
        let v1 = VertexId::new(base + 1);
        let v2 = VertexId::new(base + 2);
        let v3 = VertexId::new(base + 3);
        [
            FaceData::new(v0, v2, v1, RegionId::INVALID),
            FaceData::new(v0, v1, v3, RegionId::INVALID),
            FaceData::new(v1, v2, v3, RegionId::INVALID),
            FaceData::new(v2, v0, v3, RegionId::INVALID),
        ]
    }

    fn make_two_component_store() -> FaceStore {
        let mut store = FaceStore::new();
        for face in tetra_faces(0) {
            let _ = store.push(face);
        }
        for face in tetra_faces(4) {
            let _ = store.push(face);
        }
        store
    }

    #[test]
    fn fix_orientation_repairs_disconnected_component() {
        let mut faces = make_two_component_store();
        // Corrupt only the second connected component.
        faces.get_mut(FaceId::from_usize(5)).flip();
        let edges = EdgeStore::from_face_store(&faces);
        assert!(check_orientation(&faces, &edges).is_err());

        let flipped = fix_orientation(&mut faces, &edges);
        assert!(flipped > 0);
        let edges_after = EdgeStore::from_face_store(&faces);
        assert!(check_orientation(&faces, &edges_after).is_ok());
    }

    proptest! {
        #[test]
        fn fix_orientation_restores_consistency_under_random_flips(mask in prop::array::uniform8(any::<bool>())) {
            let mut faces = make_two_component_store();
            for (idx, flip) in mask.into_iter().enumerate() {
                if flip {
                    faces.get_mut(FaceId::from_usize(idx)).flip();
                }
            }

            let edges = EdgeStore::from_face_store(&faces);
            let _ = fix_orientation(&mut faces, &edges);
            let edges_after = EdgeStore::from_face_store(&faces);
            prop_assert!(check_orientation(&faces, &edges_after).is_ok());
        }
    }
}
