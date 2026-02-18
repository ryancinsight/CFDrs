//! Edge storage with half-edge connectivity.
//!
//! Each edge is stored as a canonical `(min_vertex, max_vertex)` pair with
//! references to adjacent faces. This replaces csgrs's on-demand adjacency
//! rebuilding with a persistent, incrementally-maintained structure.

use hashbrown::HashMap;

use crate::core::index::{VertexId, FaceId, EdgeId};
use crate::storage::face_store::FaceData;

/// Data stored per edge.
#[derive(Clone, Debug)]
pub struct EdgeData {
    /// The two endpoint vertex IDs (canonical: v0 < v1).
    pub vertices: (VertexId, VertexId),
    /// Faces sharing this edge (0 = boundary, 1 = boundary, 2 = manifold, >2 = non-manifold).
    pub faces: Vec<FaceId>,
}

impl EdgeData {
    /// Is this a boundary edge (shared by exactly 1 face)?
    #[inline]
    pub fn is_boundary(&self) -> bool {
        self.faces.len() == 1
    }

    /// Is this a manifold interior edge (shared by exactly 2 faces)?
    #[inline]
    pub fn is_manifold(&self) -> bool {
        self.faces.len() == 2
    }

    /// Is this a non-manifold edge (shared by >2 faces)?
    #[inline]
    pub fn is_non_manifold(&self) -> bool {
        self.faces.len() > 2
    }

    /// Valence: number of adjacent faces.
    #[inline]
    pub fn valence(&self) -> usize {
        self.faces.len()
    }
}

/// Storage for edges, built from faces.
///
/// Edges are identified by their canonical vertex pair `(min, max)`.
#[derive(Clone)]
pub struct EdgeStore {
    /// Edge data indexed by `EdgeId`.
    edges: Vec<EdgeData>,
    /// Lookup: canonical vertex pair → edge ID.
    edge_map: HashMap<(VertexId, VertexId), EdgeId>,
}

impl EdgeStore {
    /// Create an empty edge store.
    pub fn new() -> Self {
        Self {
            edges: Vec::new(),
            edge_map: HashMap::new(),
        }
    }

    /// Build the edge store from a slice of faces.
    ///
    /// This scans all face edges and constructs the edge adjacency in O(F)
    /// where F = number of faces.
    pub fn from_faces(faces: &[(FaceId, &FaceData)]) -> Self {
        let mut store = Self::new();
        store.edges.reserve(faces.len() * 2); // heuristic: E ≈ 1.5F for manifold

        for &(face_id, face) in faces {
            for (a, b) in face.edges_canonical() {
                store.register_edge(a, b, face_id);
            }
        }

        store
    }

    /// Build from a face store directly.
    pub fn from_face_store(face_store: &crate::storage::face_store::FaceStore) -> Self {
        let pairs: Vec<_> = face_store
            .iter_enumerated()
            .collect();
        Self::from_faces(&pairs)
    }

    /// Register an edge between `a` and `b` as belonging to `face`.
    fn register_edge(&mut self, a: VertexId, b: VertexId, face: FaceId) {
        let key = if a.0 <= b.0 { (a, b) } else { (b, a) };

        if let Some(&edge_id) = self.edge_map.get(&key) {
            self.edges[edge_id.as_usize()].faces.push(face);
        } else {
            let edge_id = EdgeId::from_usize(self.edges.len());
            self.edges.push(EdgeData {
                vertices: key,
                faces: vec![face],
            });
            self.edge_map.insert(key, edge_id);
        }
    }

    /// Number of edges.
    #[inline]
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// Is the store empty?
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Get edge data by ID.
    #[inline]
    pub fn get(&self, id: EdgeId) -> &EdgeData {
        &self.edges[id.as_usize()]
    }

    /// Look up an edge by its canonical vertex pair.
    pub fn find_edge(&self, a: VertexId, b: VertexId) -> Option<EdgeId> {
        let key = if a.0 <= b.0 { (a, b) } else { (b, a) };
        self.edge_map.get(&key).copied()
    }

    /// Iterate over all edges.
    pub fn iter(&self) -> impl Iterator<Item = &EdgeData> {
        self.edges.iter()
    }

    /// Iterate with IDs.
    pub fn iter_enumerated(&self) -> impl Iterator<Item = (EdgeId, &EdgeData)> {
        self.edges
            .iter()
            .enumerate()
            .map(|(i, e)| (EdgeId::from_usize(i), e))
    }

    /// All boundary edges (valence == 1).
    pub fn boundary_edges(&self) -> Vec<EdgeId> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_boundary())
            .map(|(i, _)| EdgeId::from_usize(i))
            .collect()
    }

    /// All non-manifold edges (valence > 2).
    pub fn non_manifold_edges(&self) -> Vec<EdgeId> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_non_manifold())
            .map(|(i, _)| EdgeId::from_usize(i))
            .collect()
    }

    /// Count boundary edges.
    pub fn boundary_edge_count(&self) -> usize {
        self.edges.iter().filter(|e| e.is_boundary()).count()
    }

    /// Clear all edges.
    pub fn clear(&mut self) {
        self.edges.clear();
        self.edge_map.clear();
    }
}

impl Default for EdgeStore {
    fn default() -> Self {
        Self::new()
    }
}
