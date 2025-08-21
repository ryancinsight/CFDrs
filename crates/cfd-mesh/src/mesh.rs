//! Mesh data structures and operations.

use nalgebra::{Point3, RealField};
use std::collections::HashSet;

/// Mesh element representation
#[derive(Debug, Clone)]
pub struct Element<T: RealField + Copy> {
    /// Indices of vertices that form this element
    pub vertices: Vec<usize>,
    /// Type of the element
    pub element_type: ElementType,
    /// Phantom data for type parameter
    _phantom: std::marker::PhantomData<T>,
}

/// Element type classification for mesh cells
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ElementType {
    /// 1D line element
    Line,
    /// 2D triangle element
    Triangle,
    /// 2D quadrilateral element
    Quadrilateral,
    /// 3D tetrahedron element
    Tetrahedron,
    /// 3D hexahedron (cube) element
    Hexahedron,
    /// 3D pentahedron (triangular prism) element
    Pentahedron,
    /// 3D pyramid element
    Pyramid,
    /// General polyhedron element
    Polyhedron,
}

impl ElementType {
    /// Get the expected number of vertices for this element type
    pub fn expected_vertex_count(&self) -> usize {
        match self {
            ElementType::Line => 2,
            ElementType::Triangle => 3,
            ElementType::Quadrilateral => 4,
            ElementType::Tetrahedron => 4,
            ElementType::Hexahedron => 8,
            ElementType::Pentahedron => 6,
            ElementType::Pyramid => 5,
            ElementType::Polyhedron => 0, // Variable
        }
    }

    /// Check if this is a 3D element type
    pub fn is_3d(&self) -> bool {
        matches!(self, 
            ElementType::Tetrahedron | 
            ElementType::Hexahedron | 
            ElementType::Pentahedron | 
            ElementType::Pyramid | 
            ElementType::Polyhedron
        )
    }

    /// Check if this is a 2D element type
    pub fn is_2d(&self) -> bool {
        matches!(self, ElementType::Triangle | ElementType::Quadrilateral)
    }

    /// Infer element type from vertex count (fallback when type is unknown)
    pub fn infer_from_vertex_count(count: usize) -> Self {
        match count {
            2 => ElementType::Line,
            3 => ElementType::Triangle,
            4 => ElementType::Tetrahedron, // Default to 3D for 4 vertices
            5 => ElementType::Pyramid,
            6 => ElementType::Pentahedron,
            8 => ElementType::Hexahedron,
            _ => ElementType::Polyhedron,
        }
    }
}

/// Vertex in a mesh
#[derive(Debug, Clone)]
pub struct Vertex<T: RealField + Copy> {
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

/// Cell in a mesh with explicit element type
#[derive(Debug, Clone)]
pub struct Cell {
    /// Face indices
    pub faces: Vec<usize>,
    /// Cell ID
    pub id: usize,
    /// Explicit element type
    pub element_type: ElementType,
}

impl Cell {
    /// Create a new cell with explicit element type
    pub fn new(id: usize, faces: Vec<usize>, element_type: ElementType) -> Self {
        Self {
            faces,
            id,
            element_type,
        }
    }

    /// Get unique vertices for this cell from the mesh
    /// 
    /// This method correctly handles the vertex collection by deduplicating
    /// vertices from all faces, fixing the issue where flattening face vertices
    /// resulted in many duplicates.
    pub fn unique_vertices<'a, T: RealField + Copy>(&self, mesh: &'a Mesh<T>) -> Vec<&'a Point3<T>> {
        let mut vertex_indices = HashSet::new();
        
        // Collect unique vertex indices from all faces
        for &face_idx in &self.faces {
            if let Some(face) = mesh.faces.get(face_idx) {
                for &vertex_idx in &face.vertices {
                    vertex_indices.insert(vertex_idx);
                }
            }
        }
        
        // Convert indices to vertex positions
        vertex_indices
            .into_iter()
            .filter_map(|idx| mesh.vertices.get(idx).map(|v| &v.position))
            .collect()
    }

    /// Get unique vertex indices for this cell
    pub fn unique_vertex_indices<T: RealField + Copy>(&self, mesh: &Mesh<T>) -> Vec<usize> {
        let mut vertex_indices = HashSet::new();
        
        for &face_idx in &self.faces {
            if let Some(face) = mesh.faces.get(face_idx) {
                for &vertex_idx in &face.vertices {
                    vertex_indices.insert(vertex_idx);
                }
            }
        }
        
        vertex_indices.into_iter().collect()
    }

    /// Validate that the cell has the expected number of vertices for its type
    pub fn validate_vertex_count<T: RealField + Copy>(&self, mesh: &Mesh<T>) -> bool {
        let actual_count = self.unique_vertex_indices(mesh).len();
        let expected_count = self.element_type.expected_vertex_count();
        
        if expected_count == 0 {
            // Polyhedron can have any number of vertices
            actual_count >= 4
        } else {
            actual_count == expected_count
        }
    }
}

/// Mesh topology
#[derive(Debug, Clone, Default)]
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
pub struct Mesh<T: RealField + Copy> {
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

impl<T: RealField + Copy> Mesh<T> {
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

    /// Get mesh elements (cells as elements)
    pub fn elements(&self) -> Vec<Element<T>> {
        self.cells.iter().map(|cell| {
            // Collect all unique vertices from the cell's faces
            let mut vertex_indices = HashSet::new();
            for &face_idx in &cell.faces {
                if let Some(face) = self.faces.get(face_idx) {
                    for &vertex_idx in &face.vertices {
                        vertex_indices.insert(vertex_idx);
                    }
                }
            }
            
            Element {
                vertices: vertex_indices.into_iter().collect(),
                element_type: cell.element_type,
                _phantom: std::marker::PhantomData,
            }
        }).collect()
    }
    
    /// Validate mesh consistency
    pub fn validate(&self) -> Result<(), String> {
        // Check that all cells have valid element types and vertex counts
        for cell in &self.cells {
            if !cell.validate_vertex_count(self) {
                return Err(format!(
                    "Cell {} has incorrect vertex count for element type {:?}",
                    cell.id.clone(), cell.element_type
                ));
            }
        }
        
        Ok(())
    }
}

impl<T: RealField + Copy> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}