//! Connected component analysis via BFS on face adjacency.
//!
//! # Algorithm
//!
//! Repeated BFS seeds from the lowest-index unvisited face, expanding via
//! the face-adjacency graph until all reachable faces are collected into one
//! component.  The outer loop guarantees every disconnected component is found.
//!
//! # Theorem — BFS Component Completeness
//!
//! **Statement.**  Let *G = (V, E)* be an undirected graph.  BFS from any
//! seed vertex visits exactly the connected component containing that seed.
//! Repeating with a fresh seed from the unvisited set until empty partitions
//! *V* into its connected components.
//!
//! **Proof sketch.**  BFS explores all vertices reachable from the seed via
//! the queue-expansion invariant: every vertex dequeued has all its neighbours
//! enqueued (or already visited).  By induction on the BFS wavefront, this
//! reaches the full connected component.  Since each vertex is visited at most
//! once (`visited` guard), and the outer loop picks a new seed from the
//! complement of all visited vertices, every component is discovered exactly
//! once.  ∎
//!
//! # Complexity
//!
//! Time **O(V + E)**, space **O(V)** where *V* = number of faces and
//! *E* = number of adjacency edges.  The dense `Vec<bool>` visited set
//! replaces the previous `HashSet<FaceId>` to eliminate per-element hashing
//! overhead and reduce memory from ~56 bytes/entry (hash-set bucket) to
//! 1 byte/entry (bool vec) — a ~56× memory reduction for the visited set.

use crate::domain::core::index::FaceId;
use crate::domain::topology::AdjacencyGraph;
use crate::infrastructure::storage::face_store::FaceStore;

/// Identify connected components using BFS on face adjacency.
///
/// Returns a list of components, each being a vector of face IDs.
///
/// # Performance
///
/// Uses a dense `Vec<bool>` for the visited set instead of `HashSet<FaceId>`,
/// giving cache-friendly sequential access and eliminating hash overhead.
/// The `VecDeque` is allocated once and reused across components via `clear()`.
#[must_use]
pub fn connected_components(
    face_store: &FaceStore,
    adjacency: &AdjacencyGraph,
) -> Vec<Vec<FaceId>> {
    let total_faces = face_store.len();
    let mut visited = vec![false; total_faces];
    let mut components = Vec::new();
    let mut queue = std::collections::VecDeque::new();

    for (fid, _) in face_store.iter_enumerated() {
        let idx = fid.as_usize();
        if visited[idx] {
            continue;
        }

        let mut component = Vec::new();
        queue.clear();
        queue.push_back(fid);
        visited[idx] = true;

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for &neighbor in adjacency.face_neighbors(current) {
                let ni = neighbor.as_usize();
                if !visited[ni] {
                    visited[ni] = true;
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    components
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::index::RegionId;
    use crate::domain::topology::AdjacencyGraph;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use crate::infrastructure::storage::face_store::{FaceData, FaceStore};

    fn tetra_faces(base: u32) -> [FaceData; 4] {
        let v0 = crate::domain::core::index::VertexId::new(base);
        let v1 = crate::domain::core::index::VertexId::new(base + 1);
        let v2 = crate::domain::core::index::VertexId::new(base + 2);
        let v3 = crate::domain::core::index::VertexId::new(base + 3);
        [
            FaceData::new(v0, v2, v1, RegionId::INVALID),
            FaceData::new(v0, v1, v3, RegionId::INVALID),
            FaceData::new(v1, v2, v3, RegionId::INVALID),
            FaceData::new(v2, v0, v3, RegionId::INVALID),
        ]
    }

    /// Single tetrahedron — must produce exactly 1 connected component
    /// containing all 4 faces.
    ///
    /// # Theorem — BFS Single Component (trivial case)
    ///
    /// A tetrahedron has 4 mutually-edge-adjacent faces.  BFS from any
    /// face reaches all others within 2 hops (diameter of K₄ = 1).
    /// Therefore `connected_components` returns exactly `[{f0,f1,f2,f3}]`.  ∎
    #[test]
    fn single_tetrahedron_one_component() {
        let mut store = FaceStore::new();
        for face in tetra_faces(0) {
            store.push(face);
        }
        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);
        assert_eq!(comps.len(), 1);
        assert_eq!(comps[0].len(), 4);
    }

    /// Two disjoint tetrahedra — must produce exactly 2 connected components.
    ///
    /// # Theorem — BFS Disconnected Graph
    ///
    /// Two tetrahedra on disjoint vertex sets share no edges so the face
    /// adjacency graph has two connected components.  BFS from a face in
    /// tetrahedron A can never reach a face in tetrahedron B.  The outer
    /// loop of `connected_components` seeds a second BFS for B.  ∎
    #[test]
    fn two_disjoint_tetrahedra_two_components() {
        let mut store = FaceStore::new();
        for face in tetra_faces(0) {
            store.push(face);
        }
        for face in tetra_faces(10) {
            store.push(face);
        }
        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);
        assert_eq!(comps.len(), 2);
        assert_eq!(comps[0].len(), 4);
        assert_eq!(comps[1].len(), 4);
    }

    /// Three disjoint tetrahedra — BFS must find exactly 3 components.
    ///
    /// # Theorem — k-Component Discovery
    ///
    /// For k disjoint connected subgraphs, the outer loop of
    /// `connected_components` seeds exactly k BFS traversals.  Each
    /// traversal exhausts one subgraph (by the BFS component completeness
    /// theorem), and the outer loop terminates after all faces are visited.  ∎
    #[test]
    fn three_disjoint_tetrahedra_three_components() {
        let mut store = FaceStore::new();
        for face in tetra_faces(0) {
            store.push(face);
        }
        for face in tetra_faces(10) {
            store.push(face);
        }
        for face in tetra_faces(20) {
            store.push(face);
        }
        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);
        assert_eq!(comps.len(), 3);
        for comp in &comps {
            assert_eq!(comp.len(), 4);
        }
    }

    /// Empty face store — must return zero components (not panic).
    ///
    /// # Theorem — Empty Graph
    ///
    /// A graph with V = 0 has zero connected components by convention.
    /// The outer loop of `connected_components` iterates zero times.  ∎
    #[test]
    fn empty_mesh_zero_components() {
        let store = FaceStore::new();
        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);
        assert_eq!(comps.len(), 0);
    }

    /// Two tetrahedra sharing one face (vertex-connected) — must still
    /// produce 1 component because they share an edge through the common
    /// face.
    ///
    /// # Theorem — Shared Face Connectivity
    ///
    /// If two polyhedra share a face, the shared face's edges connect
    /// faces from both polyhedra in the adjacency graph, yielding a
    /// single connected component.  ∎
    #[test]
    fn shared_face_one_component() {
        let v0 = crate::domain::core::index::VertexId::new(0);
        let v1 = crate::domain::core::index::VertexId::new(1);
        let v2 = crate::domain::core::index::VertexId::new(2);
        let v3 = crate::domain::core::index::VertexId::new(3);
        let v4 = crate::domain::core::index::VertexId::new(4);

        let mut store = FaceStore::new();
        // Tetra 1: v0, v1, v2, v3
        store.push(FaceData::new(v0, v2, v1, RegionId::INVALID));
        store.push(FaceData::new(v0, v1, v3, RegionId::INVALID));
        store.push(FaceData::new(v1, v2, v3, RegionId::INVALID));
        store.push(FaceData::new(v2, v0, v3, RegionId::INVALID));
        // Tetra 2: shares face (v0, v1, v2), apex v4
        store.push(FaceData::new(v0, v1, v4, RegionId::INVALID));
        store.push(FaceData::new(v1, v2, v4, RegionId::INVALID));
        store.push(FaceData::new(v2, v0, v4, RegionId::INVALID));

        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);
        assert_eq!(comps.len(), 1);
        assert_eq!(comps[0].len(), 7);
    }

    /// All faces in each component must reference valid face IDs.
    #[test]
    fn component_face_ids_are_valid() {
        let mut store = FaceStore::new();
        for face in tetra_faces(0) {
            store.push(face);
        }
        for face in tetra_faces(10) {
            store.push(face);
        }
        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);
        for comp in &comps {
            for &fid in comp {
                assert!(fid.as_usize() < store.len());
            }
        }
    }

    /// The union of all components must cover every face exactly once.
    ///
    /// # Theorem — BFS Partition Property
    ///
    /// The connected components returned by `connected_components` form
    /// a partition of the face set: every face appears in exactly one
    /// component, and every component is non-empty.
    ///
    /// **Proof sketch.** Each face is visited at most once (the `visited`
    /// guard prevents re-insertion).  The outer loop seeds BFS from every
    /// unvisited face, so every face is eventually processed.  Each face
    /// is added to exactly one component vector.  ∎
    #[test]
    fn components_partition_all_faces() {
        let mut store = FaceStore::new();
        for face in tetra_faces(0) {
            store.push(face);
        }
        for face in tetra_faces(10) {
            store.push(face);
        }
        for face in tetra_faces(20) {
            store.push(face);
        }
        let edges = EdgeStore::from_face_store(&store);
        let adj = AdjacencyGraph::build(&store, &edges);
        let comps = connected_components(&store, &adj);

        let total: usize = comps.iter().map(|c| c.len()).sum();
        assert_eq!(total, store.len());

        // No face appears in more than one component.
        let mut all_fids: Vec<FaceId> = comps.into_iter().flatten().collect();
        all_fids.sort_unstable();
        all_fids.dedup();
        assert_eq!(all_fids.len(), store.len());
    }
}
