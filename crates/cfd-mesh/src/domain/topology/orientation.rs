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
/// Uses BFS from an arbitrary seed face, flipping faces that have inconsistent
/// orientation relative to their neighbors.
pub fn fix_orientation(face_store: &mut FaceStore, edge_store: &EdgeStore) -> usize {
    use hashbrown::HashSet;

    let n_faces = face_store.len();
    if n_faces == 0 {
        return 0;
    }

    let mut visited = HashSet::with_capacity(n_faces);
    let mut flipped = 0usize;
    let mut queue = std::collections::VecDeque::new();

    // Seed with face 0
    let seed = FaceId::from_usize(0);
    queue.push_back(seed);
    visited.insert(seed);

    while let Some(current) = queue.pop_front() {
        let current_face = *face_store.get(current);

        // For each edge of the current face, check the neighbor
        for (ea, eb) in current_face.edges_canonical() {
            if let Some(eid) = edge_store.find_edge(ea, eb) {
                let edge = edge_store.get(eid);
                for &neighbor_fid in &edge.faces {
                    if neighbor_fid == current || !visited.insert(neighbor_fid) {
                        continue;
                    }

                    let dir_current = directed_edge_order(current_face.vertices, ea, eb);
                    let neighbor_face = face_store.get(neighbor_fid);
                    let dir_neighbor = directed_edge_order(neighbor_face.vertices, ea, eb);

                    if let (Some(dc), Some(dn)) = (dir_current, dir_neighbor) {
                        if dc == dn {
                            // Same direction → flip neighbor
                            face_store.get_mut(neighbor_fid).flip();
                            flipped += 1;
                        }
                    }

                    queue.push_back(neighbor_fid);
                }
            }
        }
    }

    flipped
}
