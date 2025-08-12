//! Mesh data structures and operations.

use nalgebra::{Point3, RealField};

/// Vertex in a mesh
#[derive(Debug, Clone)]
pub struct Vertex<T: RealField> {
    /// Position
    pub position: Point3<T>,
    /// Vertex ID
    pub id: usize,
}

/// Edge in a mesh
#[derive(Debug, Clone)]
pub struct Edge {
    /// Vertex indices
    pub vertices: [usize; 2],
    /// Edge ID
    pub id: usize,
}

/// Face in a mesh
#[derive(Debug, Clone)]
pub struct Face {
    /// Vertex indices
    pub vertices: Vec<usize>,
    /// Face ID
    pub id: usize,
}

/// Cell in a mesh
#[derive(Debug, Clone)]
pub struct Cell {
    /// Face indices
    pub faces: Vec<usize>,
    /// Cell ID
    pub id: usize,
}

/// Mesh topology
#[derive(Debug, Clone)]
pub struct MeshTopology {
    /// Number of vertices
    pub num_vertices: usize,
    /// Number of edges
    pub num_edges: usize,
    /// Number of faces
    pub num_faces: usize,
    /// Number of cells
    pub num_cells: usize,
}

/// Generic mesh structure
#[derive(Debug, Clone)]
pub struct Mesh<T: RealField> {
    /// Vertices
    pub vertices: Vec<Vertex<T>>,
    /// Edges
    pub edges: Vec<Edge>,
    /// Faces
    pub faces: Vec<Face>,
    /// Cells
    pub cells: Vec<Cell>,
    /// Topology information
    pub topology: MeshTopology,
}

impl<T: RealField> Mesh<T> {
    /// Create a new empty mesh
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            cells: Vec::new(),
            topology: MeshTopology {
                num_vertices: 0,
                num_edges: 0,
                num_faces: 0,
                num_cells: 0,
            },
        }
    }
    
    /// Update mesh topology counts
    pub fn update_topology(&mut self) {
        self.topology = MeshTopology {
            num_vertices: self.vertices.len(),
            num_edges: self.edges.len(),
            num_faces: self.faces.len(),
            num_cells: self.cells.len(),
        };
    }
}

impl<T: RealField> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}