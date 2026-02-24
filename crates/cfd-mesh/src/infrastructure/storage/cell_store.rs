//! Volumetric cell storage for [`crate::domain::mesh::IndexedMesh`].
//!
//! Each cell references its bounding faces (by [`FaceId`]) and directly lists
//! its vertex indices (by [`VertexId`]) for O(1) element connectivity lookups.

use crate::domain::core::index::{CellId, FaceId, VertexId};
use crate::domain::topology::ElementType;

/// A volumetric mesh cell referencing faces and vertices by typed index.
#[derive(Clone, Debug)]
pub struct CellData {
    /// Face indices forming this cell's boundary.
    pub faces: Vec<FaceId>,
    /// Element type (tetrahedron, hexahedron, etc.).
    pub element_type: ElementType,
    /// Vertex indices forming this cell (flat, for direct connectivity access).
    pub vertex_ids: Vec<VertexId>,
}

impl CellData {
    /// Create a tetrahedral cell from four triangular face IDs.
    pub fn tetrahedron(f0: FaceId, f1: FaceId, f2: FaceId, f3: FaceId) -> Self {
        Self {
            faces: vec![f0, f1, f2, f3],
            element_type: ElementType::Tetrahedron,
            vertex_ids: Vec::new(),
        }
    }

    /// Create a tetrahedral cell with explicit vertex IDs.
    pub fn tetrahedron_with_vertices(
        f0: FaceId,
        f1: FaceId,
        f2: FaceId,
        f3: FaceId,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        v3: VertexId,
    ) -> Self {
        Self {
            faces: vec![f0, f1, f2, f3],
            element_type: ElementType::Tetrahedron,
            vertex_ids: vec![v0, v1, v2, v3],
        }
    }

    /// Create a hexahedral cell from six face IDs.
    pub fn hexahedron(faces: [FaceId; 6]) -> Self {
        Self {
            faces: faces.to_vec(),
            element_type: ElementType::Hexahedron,
            vertex_ids: Vec::new(),
        }
    }
}

/// Storage for volumetric mesh cells.
pub struct CellStore {
    cells: Vec<CellData>,
}

impl CellStore {
    /// Create an empty cell store.
    pub fn new() -> Self {
        Self { cells: Vec::new() }
    }

    /// Number of cells.
    #[inline]
    pub fn len(&self) -> usize {
        self.cells.len()
    }

    /// Is the store empty?
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// Add a cell, returning its ID.
    pub fn push(&mut self, cell: CellData) -> CellId {
        let id = CellId::from_usize(self.cells.len());
        self.cells.push(cell);
        id
    }

    /// Get cell data by ID.
    #[inline]
    pub fn get(&self, id: CellId) -> &CellData {
        &self.cells[id.as_usize()]
    }

    /// Iterate all cells.
    pub fn iter(&self) -> impl Iterator<Item = &CellData> {
        self.cells.iter()
    }

    /// Iterate with IDs.
    pub fn iter_enumerated(&self) -> impl Iterator<Item = (CellId, &CellData)> {
        self.cells
            .iter()
            .enumerate()
            .map(|(i, c)| (CellId::from_usize(i), c))
    }

    /// Access the underlying slice.
    pub fn as_slice(&self) -> &[CellData] {
        &self.cells
    }

    /// Clear all cells.
    pub fn clear(&mut self) {
        self.cells.clear();
    }
}

impl Default for CellStore {
    fn default() -> Self {
        Self::new()
    }
}
