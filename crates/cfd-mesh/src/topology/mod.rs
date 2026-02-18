//! Topological queries on the indexed mesh.
//!
//! Built from the indexed storage layer. No redundant adjacency rebuilding —
//! the `EdgeStore` provides persistent adjacency.

pub mod adjacency;
pub mod connectivity;
pub mod manifold;
pub mod orientation;

pub use adjacency::AdjacencyGraph;

// ── Mesh-element types needed by cfd-3d ──────────────────────────────────────

use nalgebra::Point3;

/// Element type of a mesh cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementType {
    Tetrahedron,
    Hexahedron,
    Triangle,
    Quadrilateral,
    Wedge,
    Pyramid,
}

/// A mesh vertex carrying a 3-D position.
#[derive(Debug, Clone)]
pub struct Vertex<T: nalgebra::Scalar + Copy> {
    /// 3-D position in space.
    pub position: Point3<T>,
}

impl<T: nalgebra::Scalar + Copy + nalgebra::RealField> Vertex<T> {
    /// Create a vertex at `position`.
    pub fn new(position: Point3<T>) -> Self {
        Self { position }
    }
}

/// A mesh face (polygon) referencing vertex indices.
#[derive(Debug, Clone)]
pub struct Face {
    /// Ordered list of vertex indices that define this face.
    pub vertices: Vec<usize>,
}

impl Face {
    /// Triangular face.
    pub fn triangle(v0: usize, v1: usize, v2: usize) -> Self {
        Self { vertices: vec![v0, v1, v2] }
    }
    /// Quadrilateral face.
    pub fn quad(v0: usize, v1: usize, v2: usize, v3: usize) -> Self {
        Self { vertices: vec![v0, v1, v2, v3] }
    }
}

/// A volumetric mesh cell referencing face indices.
#[derive(Debug, Clone)]
pub struct Cell {
    /// Indices into the mesh face list.
    pub faces: Vec<usize>,
    /// Element type.
    pub element_type: ElementType,
    /// Indices of all vertices forming this cell (flat, for convenience).
    pub vertex_ids: Vec<usize>,
}

impl Cell {
    /// Tetrahedral cell — four triangular faces.
    pub fn tetrahedron(f0: usize, f1: usize, f2: usize, f3: usize) -> Self {
        Self {
            faces: vec![f0, f1, f2, f3],
            element_type: ElementType::Tetrahedron,
            vertex_ids: Vec::new(),
        }
    }
    /// Hexahedral cell — six quadrilateral faces.
    pub fn hexahedron(
        f0: usize, f1: usize, f2: usize,
        f3: usize, f4: usize, f5: usize,
    ) -> Self {
        Self {
            faces: vec![f0, f1, f2, f3, f4, f5],
            element_type: ElementType::Hexahedron,
            vertex_ids: Vec::new(),
        }
    }
}
