//! Edge storage with half-edge connectivity.
//!
//! Each edge is stored as a canonical `(min_vertex, max_vertex)` pair with
//! references to adjacent faces. This replaces csgrs's on-demand adjacency
//! rebuilding with a persistent, incrementally-maintained structure.

use hashbrown::HashMap;

use crate::domain::core::index::{EdgeId, FaceId, VertexId};
use crate::infrastructure::storage::face_store::FaceData;

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
    #[must_use]
    pub fn is_boundary(&self) -> bool {
        self.faces.len() == 1
    }

    /// Is this a manifold interior edge (shared by exactly 2 faces)?
    #[inline]
    #[must_use]
    pub fn is_manifold(&self) -> bool {
        self.faces.len() == 2
    }

    /// Is this a non-manifold edge (shared by >2 faces)?
    #[inline]
    #[must_use]
    pub fn is_non_manifold(&self) -> bool {
        self.faces.len() > 2
    }

    /// Valence: number of adjacent faces.
    #[inline]
    #[must_use]
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
    #[must_use]
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
    ///
    /// Capacity hint: for a closed manifold triangle mesh, E = 3F/2 by
    /// the Euler relation.  Pre-allocating both the `edges` vec and the
    /// `edge_map` hash avoids incremental rehashing during construction.
    #[must_use]
    pub fn from_faces(faces: &[(FaceId, &FaceData)]) -> Self {
        let cap = faces.len().saturating_mul(3) / 2;
        let mut store = Self {
            edges: Vec::with_capacity(cap),
            edge_map: HashMap::with_capacity(cap),
        };

        for &(face_id, face) in faces {
            for (a, b) in face.edges_canonical() {
                store.register_edge(a, b, face_id);
            }
        }

        store
    }

    /// Build from a face store directly.
    #[must_use]
    pub fn from_face_store(
        face_store: &crate::infrastructure::storage::face_store::FaceStore,
    ) -> Self {
        let pairs: Vec<_> = face_store.iter_enumerated().collect();
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
    #[must_use]
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// Is the store empty?
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Get edge data by ID.
    #[inline]
    #[must_use]
    pub fn get(&self, id: EdgeId) -> &EdgeData {
        &self.edges[id.as_usize()]
    }

    /// Look up an edge by its canonical vertex pair.
    #[must_use]
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
    #[must_use]
    pub fn boundary_edges(&self) -> Vec<EdgeId> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_boundary())
            .map(|(i, _)| EdgeId::from_usize(i))
            .collect()
    }

    /// All non-manifold edges (valence > 2).
    #[must_use]
    pub fn non_manifold_edges(&self) -> Vec<EdgeId> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_non_manifold())
            .map(|(i, _)| EdgeId::from_usize(i))
            .collect()
    }

    /// Count boundary edges.
    #[must_use]
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::index::RegionId;
    use crate::infrastructure::storage::face_store::{FaceData, FaceStore};

    fn vid(n: u32) -> VertexId {
        VertexId::new(n)
    }

    fn tetra_store() -> FaceStore {
        let mut s = FaceStore::new();
        s.push(FaceData::new(vid(0), vid(2), vid(1), RegionId::INVALID));
        s.push(FaceData::new(vid(0), vid(1), vid(3), RegionId::INVALID));
        s.push(FaceData::new(vid(1), vid(2), vid(3), RegionId::INVALID));
        s.push(FaceData::new(vid(2), vid(0), vid(3), RegionId::INVALID));
        s
    }

    /// A tetrahedron has exactly 6 edges: $\binom{4}{2} = 6$.
    ///
    /// # Theorem — Tetrahedron Edge Count
    ///
    /// **Statement.** A tetrahedron on 4 vertices has exactly 6 edges.
    ///
    /// **Proof.** Each pair of distinct vertices is connected by an edge.
    /// There are $\binom{4}{2} = 6$ such pairs.  ∎
    #[test]
    fn tet_edge_count_is_six() {
        let fs = tetra_store();
        let es = EdgeStore::from_face_store(&fs);
        assert_eq!(es.len(), 6);
    }

    /// Every edge of a closed tetrahedron is manifold (shared by exactly 2 faces).
    ///
    /// # Theorem — Tetrahedron Manifold Edges
    ///
    /// **Statement.** In a closed tetrahedron, each of the 6 edges is shared
    /// by exactly 2 of the 4 triangular faces.
    ///
    /// **Proof.** Each edge `{a, b}` is the intersection of the two faces
    /// containing both `a` and `b`.  Since $\binom{4-2}{1} + 1 = 2$ faces
    /// contain any given edge (choose 1 of the remaining 2 vertices for
    /// each face), every edge has valence 2.  ∎
    #[test]
    fn tet_all_edges_manifold() {
        let fs = tetra_store();
        let es = EdgeStore::from_face_store(&fs);
        for edge in es.iter() {
            assert!(
                edge.is_manifold(),
                "edge {:?} has valence {} (expected 2)",
                edge.vertices,
                edge.valence()
            );
        }
        assert_eq!(es.boundary_edge_count(), 0);
    }

    /// `find_edge` is canonical — order of arguments does not matter.
    ///
    /// # Theorem — Canonical Edge Lookup
    ///
    /// **Statement.** `find_edge(a, b)` and `find_edge(b, a)` return the
    /// same `EdgeId`.
    ///
    /// **Proof.** Both calls compute the canonical key
    /// `(min(a, b), max(a, b))` and look it up in the same `edge_map`.  ∎
    #[test]
    fn find_edge_canonical_order() {
        let fs = tetra_store();
        let es = EdgeStore::from_face_store(&fs);
        for edge in es.iter() {
            let (a, b) = edge.vertices;
            let id_ab = es.find_edge(a, b);
            let id_ba = es.find_edge(b, a);
            assert_eq!(id_ab, id_ba, "canonical order violated for {:?}", (a, b));
            assert!(id_ab.is_some());
        }
    }

    /// A single triangle has 3 boundary edges (valence 1).
    #[test]
    fn single_triangle_boundary_edges() {
        let mut fs = FaceStore::new();
        fs.push(FaceData::new(vid(0), vid(1), vid(2), RegionId::INVALID));
        let es = EdgeStore::from_face_store(&fs);
        assert_eq!(es.len(), 3);
        assert_eq!(es.boundary_edge_count(), 3);
        for edge in es.iter() {
            assert!(edge.is_boundary());
        }
    }

    /// Non-manifold edge detection: 3 faces sharing one edge.
    ///
    /// # Theorem — Non-Manifold Detection
    ///
    /// **Statement.** An edge shared by $k > 2$ faces is classified as
    /// non-manifold (valence $k$).
    ///
    /// **Proof.** `register_edge` pushes each face_id into the edge's
    /// face list.  After processing all faces, `edge.faces.len() == k`,
    /// and `is_non_manifold()` returns `k > 2`.  ∎
    #[test]
    fn non_manifold_edge_detected() {
        let mut fs = FaceStore::new();
        // 3 triangles sharing edge (v0, v1)
        fs.push(FaceData::new(vid(0), vid(1), vid(2), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(1), vid(3), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(1), vid(4), RegionId::INVALID));
        let es = EdgeStore::from_face_store(&fs);
        let nm = es.non_manifold_edges();
        assert_eq!(nm.len(), 1, "exactly one non-manifold edge expected");
        let e = es.get(nm[0]);
        assert_eq!(e.valence(), 3);
        assert!(e.is_non_manifold());
    }

    /// Empty face store produces empty edge store.
    #[test]
    fn empty_face_store_empty_edges() {
        let fs = FaceStore::new();
        let es = EdgeStore::from_face_store(&fs);
        assert!(es.is_empty());
        assert_eq!(es.len(), 0);
        assert_eq!(es.boundary_edge_count(), 0);
    }

    /// Cube (12 triangles, 8 vertices) has 18 edges — all manifold.
    ///
    /// # Theorem — Cube Edge Count
    ///
    /// **Statement.** A triangulated cube with 8 vertices, 12 triangles,
    /// and 6 diagonal edges has E = 18.  By Euler: V - E + F = 2 →
    /// 8 - E + 12 = 2 → E = 18.
    ///
    /// **Proof.** Direct application of the Euler formula for a closed
    /// genus-0 surface.  ∎
    #[test]
    fn cube_edge_count_and_manifold() {
        let mut fs = FaceStore::new();
        // Front z=0: v0(0,0,0) v1(1,0,0) v2(1,1,0) v3(0,1,0)
        // Back z=1: v4(0,0,1) v5(1,0,1) v6(1,1,1) v7(0,1,1)
        fs.push(FaceData::new(vid(0), vid(1), vid(2), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(2), vid(3), RegionId::INVALID));
        fs.push(FaceData::new(vid(4), vid(6), vid(5), RegionId::INVALID));
        fs.push(FaceData::new(vid(4), vid(7), vid(6), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(5), vid(1), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(4), vid(5), RegionId::INVALID));
        fs.push(FaceData::new(vid(3), vid(2), vid(6), RegionId::INVALID));
        fs.push(FaceData::new(vid(3), vid(6), vid(7), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(3), vid(7), RegionId::INVALID));
        fs.push(FaceData::new(vid(0), vid(7), vid(4), RegionId::INVALID));
        fs.push(FaceData::new(vid(1), vid(5), vid(6), RegionId::INVALID));
        fs.push(FaceData::new(vid(1), vid(6), vid(2), RegionId::INVALID));
        let es = EdgeStore::from_face_store(&fs);
        assert_eq!(es.len(), 18);
        assert_eq!(es.boundary_edge_count(), 0);
        for edge in es.iter() {
            assert!(
                edge.is_manifold(),
                "cube edge {:?} has valence {} (expected 2)",
                edge.vertices,
                edge.valence()
            );
        }
    }

    /// `clear()` resets both the edge vec and the edge map.
    #[test]
    fn clear_resets_store() {
        let fs = tetra_store();
        let mut es = EdgeStore::from_face_store(&fs);
        assert!(!es.is_empty());
        es.clear();
        assert!(es.is_empty());
        assert_eq!(es.len(), 0);
    }
}
