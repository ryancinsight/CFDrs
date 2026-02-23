//! Connected component analysis.

use hashbrown::HashSet;

use crate::core::index::FaceId;
use crate::topology::AdjacencyGraph;
use crate::storage::face_store::FaceStore;

/// Identify connected components using BFS on face adjacency.
///
/// Returns a list of components, each being a set of face IDs.
pub fn connected_components(
    face_store: &FaceStore,
    adjacency: &AdjacencyGraph,
) -> Vec<Vec<FaceId>> {
    let total_faces = face_store.len();
    let mut visited = HashSet::with_capacity(total_faces);
    let mut components = Vec::new();

    for (fid, _) in face_store.iter_enumerated() {
        if visited.contains(&fid) {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(fid);
        visited.insert(fid);

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for &neighbor in adjacency.face_neighbors(current) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    components
}
