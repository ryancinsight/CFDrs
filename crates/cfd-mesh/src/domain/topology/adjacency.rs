//! Vertex and face adjacency graph.
//!
//! Built from the `EdgeStore` — O(V+E) construction, O(1) lookups.

use hashbrown::HashMap;

use crate::domain::core::index::{FaceId, VertexId};
use crate::infrastructure::storage::edge_store::EdgeStore;
use crate::infrastructure::storage::face_store::FaceStore;

/// Pre-built adjacency graph for vertex-vertex and vertex-face queries.
pub struct AdjacencyGraph {
    /// vertex → list of adjacent vertices (1-ring neighborhood).
    vertex_neighbors: HashMap<VertexId, Vec<VertexId>>,
    /// vertex → list of incident faces.
    vertex_faces: HashMap<VertexId, Vec<FaceId>>,
    /// face → list of adjacent faces (sharing an edge).
    face_neighbors: HashMap<FaceId, Vec<FaceId>>,
}

impl AdjacencyGraph {
    /// Build the adjacency graph from edge and face stores.
    pub fn build(face_store: &FaceStore, edge_store: &EdgeStore) -> Self {
        let mut vertex_neighbors: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
        let mut vertex_faces: HashMap<VertexId, Vec<FaceId>> = HashMap::new();
        let mut face_neighbors: HashMap<FaceId, Vec<FaceId>> = HashMap::new();

        // Build vertex → faces from face store
        for (fid, face) in face_store.iter_enumerated() {
            for &vid in &face.vertices {
                vertex_faces.entry(vid).or_default().push(fid);
            }
        }

        // Build vertex neighbors and face neighbors from edge store
        for edge in edge_store.iter() {
            let (a, b) = edge.vertices;
            vertex_neighbors.entry(a).or_default().push(b);
            vertex_neighbors.entry(b).or_default().push(a);

            // Two faces sharing an edge are neighbors
            for i in 0..edge.faces.len() {
                for j in (i + 1)..edge.faces.len() {
                    let fi = edge.faces[i];
                    let fj = edge.faces[j];
                    face_neighbors.entry(fi).or_default().push(fj);
                    face_neighbors.entry(fj).or_default().push(fi);
                }
            }
        }

        // Deduplicate
        for v in vertex_neighbors.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        for v in face_neighbors.values_mut() {
            v.sort_unstable();
            v.dedup();
        }

        Self {
            vertex_neighbors,
            vertex_faces,
            face_neighbors,
        }
    }

    /// Get the 1-ring vertex neighborhood.
    pub fn vertex_neighbors(&self, v: VertexId) -> &[VertexId] {
        self.vertex_neighbors
            .get(&v)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Get faces incident to a vertex.
    pub fn vertex_faces(&self, v: VertexId) -> &[FaceId] {
        self.vertex_faces
            .get(&v)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Get faces neighboring a given face (sharing an edge).
    pub fn face_neighbors(&self, f: FaceId) -> &[FaceId] {
        self.face_neighbors
            .get(&f)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Vertex valence (number of adjacent vertices).
    pub fn vertex_valence(&self, v: VertexId) -> usize {
        self.vertex_neighbors(v).len()
    }
}
