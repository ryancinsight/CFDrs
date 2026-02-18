//! Triangular face storage.
//!
//! Each face stores three `VertexId` indices and an optional `RegionId`.
//! This is the indexed-mesh counterpart to csgrs's `Vec<Polygon<S>>` —
//! instead of each polygon owning its vertices, faces reference shared
//! vertices in the `VertexPool`.

use crate::core::index::{VertexId, FaceId, RegionId};

/// A triangular face referencing three vertices by index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceData {
    /// The three vertex indices in counter-clockwise winding order.
    pub vertices: [VertexId; 3],
    /// Optional region tag (channel ID, junction ID, boundary, etc.).
    pub region: RegionId,
}

impl FaceData {
    /// Create a new face.
    pub const fn new(v0: VertexId, v1: VertexId, v2: VertexId, region: RegionId) -> Self {
        Self {
            vertices: [v0, v1, v2],
            region,
        }
    }

    /// Create a face with no region tag.
    pub const fn untagged(v0: VertexId, v1: VertexId, v2: VertexId) -> Self {
        Self {
            vertices: [v0, v1, v2],
            region: RegionId::INVALID,
        }
    }

    /// Flip the winding order (swap v1 ↔ v2).
    pub fn flip(&mut self) {
        self.vertices.swap(1, 2);
    }

    /// Return a flipped copy.
    pub const fn flipped(&self) -> Self {
        Self {
            vertices: [self.vertices[0], self.vertices[2], self.vertices[1]],
            region: self.region,
        }
    }

    /// The three edges as ordered vertex-ID pairs `(smaller, larger)`.
    ///
    /// Canonical ordering ensures edge identity regardless of face winding.
    pub fn edges_canonical(&self) -> [(VertexId, VertexId); 3] {
        let [a, b, c] = self.vertices;
        [
            canonical_edge(a, b),
            canonical_edge(b, c),
            canonical_edge(c, a),
        ]
    }

    /// The three directed edges (preserving winding).
    pub fn edges_directed(&self) -> [(VertexId, VertexId); 3] {
        let [a, b, c] = self.vertices;
        [(a, b), (b, c), (c, a)]
    }

    /// Check if this face contains a specific vertex.
    pub fn contains_vertex(&self, v: VertexId) -> bool {
        self.vertices.contains(&v)
    }

    /// The vertex opposite to a given edge.
    pub fn opposite_vertex(&self, edge: (VertexId, VertexId)) -> Option<VertexId> {
        for &v in &self.vertices {
            if v != edge.0 && v != edge.1 {
                return Some(v);
            }
        }
        None
    }
}

/// Canonical edge: always `(min, max)`.
#[inline]
pub fn canonical_edge(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a.0 <= b.0 { (a, b) } else { (b, a) }
}

/// Storage for triangular faces.
pub struct FaceStore {
    faces: Vec<FaceData>,
}

impl FaceStore {
    /// Create an empty face store.
    pub fn new() -> Self {
        Self { faces: Vec::new() }
    }

    /// Create with capacity.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            faces: Vec::with_capacity(cap),
        }
    }

    /// Number of faces.
    #[inline]
    pub fn len(&self) -> usize {
        self.faces.len()
    }

    /// Is the store empty?
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.faces.is_empty()
    }

    /// Add a face, returning its ID.
    pub fn push(&mut self, face: FaceData) -> FaceId {
        let id = FaceId::from_usize(self.faces.len());
        self.faces.push(face);
        id
    }

    /// Remove the last face.
    pub fn pop(&mut self) -> Option<FaceData> {
        self.faces.pop()
    }

    /// Add a triangle from three vertex IDs (untagged).
    pub fn add_triangle(&mut self, v0: VertexId, v1: VertexId, v2: VertexId) -> FaceId {
        self.push(FaceData::untagged(v0, v1, v2))
    }

    /// Add a triangle with a region tag.
    pub fn add_triangle_with_region(
        &mut self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        region: RegionId,
    ) -> FaceId {
        self.push(FaceData::new(v0, v1, v2, region))
    }

    /// Get face data by ID.
    #[inline]
    pub fn get(&self, id: FaceId) -> &FaceData {
        &self.faces[id.as_usize()]
    }

    /// Get face data mutably by ID.
    #[inline]
    pub fn get_mut(&mut self, id: FaceId) -> &mut FaceData {
        &mut self.faces[id.as_usize()]
    }

    /// Iterate all faces.
    pub fn iter(&self) -> impl Iterator<Item = &FaceData> {
        self.faces.iter()
    }

    /// Mutable iterate.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut FaceData> {
        self.faces.iter_mut()
    }

    /// Iterate with IDs.
    pub fn iter_enumerated(&self) -> impl Iterator<Item = (FaceId, &FaceData)> {
        self.faces
            .iter()
            .enumerate()
            .map(|(i, f)| (FaceId::from_usize(i), f))
    }

    /// Mutable iterate with IDs.
    pub fn iter_mut_enumerated(&mut self) -> impl Iterator<Item = (FaceId, &mut FaceData)> {
        self.faces
            .iter_mut()
            .enumerate()
            .map(|(i, f)| (FaceId::from_usize(i), f))
    }

    /// Get faces belonging to a specific region.
    pub fn faces_in_region(&self, region: RegionId) -> Vec<FaceId> {
        self.faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.region == region)
            .map(|(i, _)| FaceId::from_usize(i))
            .collect()
    }

    /// Access the underlying slice.
    pub fn as_slice(&self) -> &[FaceData] {
        &self.faces
    }

    /// Clear all faces.
    pub fn clear(&mut self) {
        self.faces.clear();
    }
}

impl Default for FaceStore {
    fn default() -> Self {
        Self::new()
    }
}
