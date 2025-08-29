//! Core mesh data structures and operations
//! Following SOLID principles

use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// For now, import ElementType from cfd-core
// TODO: Move ElementType to cfd-mesh in a future refactoring
pub use cfd_core::domains::mesh_operations::ElementType;

/// Vertex in 3D space
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vertex<T: RealField + Copy> {
    /// Position in 3D space
    pub position: Point3<T>,
}

impl<T: RealField + Copy> Vertex<T> {
    /// Create a new vertex at the given position
    pub fn new(position: Point3<T>) -> Self {
        Self { position }
    }

    /// Create a new vertex from coordinates
    pub fn from_coords(x: T, y: T, z: T) -> Self {
        Self {
            position: Point3::new(x, y, z),
        }
    }

    /// Distance to another vertex
    pub fn distance_to(&self, other: &Self) -> T {
        (self.position - other.position).norm()
    }
}

/// Edge connecting two vertices
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Edge {
    /// Start vertex index
    pub start: usize,
    /// End vertex index
    pub end: usize,
}

impl Edge {
    /// Create a new edge with canonical vertex ordering
    #[must_use]
    pub fn new(v1: usize, v2: usize) -> Self {
        // Always order vertices to ensure canonical representation
        if v1 < v2 {
            Self { start: v1, end: v2 }
        } else {
            Self { start: v2, end: v1 }
        }
    }

    /// Check if edge contains vertex
    #[must_use]
    pub fn contains(&self, vertex_id: usize) -> bool {
        self.start == vertex_id || self.end == vertex_id
    }
}

/// Face defined by vertices
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Face {
    /// Vertex indices defining the face
    pub vertices: Vec<usize>,
}

impl Face {
    /// Create a triangular face
    #[must_use]
    pub fn triangle(v0: usize, v1: usize, v2: usize) -> Self {
        Self {
            vertices: vec![v0, v1, v2],
        }
    }

    /// Create a quadrilateral face
    #[must_use]
    pub fn quad(v0: usize, v1: usize, v2: usize, v3: usize) -> Self {
        Self {
            vertices: vec![v0, v1, v2, v3],
        }
    }

    /// Number of vertices
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
}

/// Cell (3D element)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    /// Cell type
    pub element_type: ElementType,
    /// Vertex indices
    pub vertices: Vec<usize>,
}

impl Cell {
    /// Create a tetrahedral cell
    #[must_use]
    pub fn tetrahedron(vertices: Vec<usize>) -> Self {
        assert_eq!(vertices.len(), 4, "Tetrahedron requires 4 vertices");
        Self {
            element_type: ElementType::Tetrahedron,
            vertices,
        }
    }

    /// Create a hexahedral cell
    #[must_use]
    pub fn hexahedron(vertices: Vec<usize>) -> Self {
        assert_eq!(vertices.len(), 8, "Hexahedron requires 8 vertices");
        Self {
            element_type: ElementType::Hexahedron,
            vertices,
        }
    }

    /// Create a prismatic cell
    #[must_use]
    pub fn prism(vertices: Vec<usize>) -> Self {
        assert_eq!(vertices.len(), 6, "Prism requires 6 vertices");
        Self {
            element_type: ElementType::Prism,
            vertices,
        }
    }

    /// Create a pyramid cell
    #[must_use]
    pub fn pyramid(vertices: Vec<usize>) -> Self {
        assert_eq!(vertices.len(), 5, "Pyramid requires 5 vertices");
        Self {
            element_type: ElementType::Pyramid,
            vertices,
        }
    }
}

/// Mesh topology information
#[derive(Debug, Clone, Default)]
pub struct MeshTopology {
    /// Vertex to cell connectivity
    pub vertex_cells: HashMap<usize, Vec<usize>>,
    /// Edge to face connectivity
    pub edge_faces: HashMap<Edge, Vec<usize>>,
    /// Face neighbors
    pub face_neighbors: HashMap<usize, Vec<usize>>,
}

/// Main mesh structure
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
    #[must_use]
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            cells: Vec::new(),
            topology: MeshTopology::default(),
        }
    }

    /// Add a vertex
    pub fn add_vertex(&mut self, vertex: Vertex<T>) -> usize {
        let id = self.vertices.len();
        self.vertices.push(vertex);
        id
    }

    /// Add an edge
    pub fn add_edge(&mut self, edge: Edge) -> usize {
        let id = self.edges.len();
        self.edges.push(edge);
        id
    }

    /// Add a face
    pub fn add_face(&mut self, face: Face) -> usize {
        let id = self.faces.len();
        self.faces.push(face);
        id
    }
    
    /// Add a cell
    pub fn add_cell(&mut self, cell: Cell) -> usize {
        let id = self.cells.len();
        self.cells.push(cell);
        id
    }

    /// Get vertex count
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Get edge count
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Get face count
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Get cell count
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Get vertices of an element
    pub fn get_element_vertices(&self, element: &Cell) -> Vec<Point3<T>> {
        element
            .vertices
            .iter()
            .filter_map(|&idx| self.vertices.get(idx))
            .map(|v| v.position)
            .collect()
    }

    /// Get faces of an element
    pub fn get_element_faces(&self, element: &Cell) -> Vec<Face> {
        match element.element_type {
            ElementType::Tetrahedron | ElementType::Tetrahedron10 => {
                // Tetrahedron has 4 triangular faces
                let v = &element.vertices;
                vec![
                    Face::triangle(v[0], v[1], v[2]), // Face 0
                    Face::triangle(v[0], v[1], v[3]), // Face 1
                    Face::triangle(v[0], v[2], v[3]), // Face 2
                    Face::triangle(v[1], v[2], v[3]), // Face 3
                ]
            }
            ElementType::Hexahedron | ElementType::Hexahedron20 => {
                // Hexahedron has 6 quadrilateral faces
                let v = &element.vertices;
                vec![
                    Face::quad(v[0], v[1], v[2], v[3]), // Bottom face
                    Face::quad(v[4], v[5], v[6], v[7]), // Top face
                    Face::quad(v[0], v[1], v[5], v[4]), // Front face
                    Face::quad(v[2], v[3], v[7], v[6]), // Back face
                    Face::quad(v[0], v[3], v[7], v[4]), // Left face
                    Face::quad(v[1], v[2], v[6], v[5]), // Right face
                ]
            }
            ElementType::Prism => {
                // Prism has 2 triangular faces and 3 quadrilateral faces
                let v = &element.vertices;
                vec![
                    Face::triangle(v[0], v[1], v[2]),   // Bottom triangular face
                    Face::triangle(v[3], v[4], v[5]),   // Top triangular face
                    Face::quad(v[0], v[1], v[4], v[3]), // Side face 1
                    Face::quad(v[1], v[2], v[5], v[4]), // Side face 2
                    Face::quad(v[2], v[0], v[3], v[5]), // Side face 3
                ]
            }
            ElementType::Pyramid => {
                // Pyramid has 1 quadrilateral base and 4 triangular faces
                let v = &element.vertices;
                vec![
                    Face::quad(v[0], v[1], v[2], v[3]), // Base
                    Face::triangle(v[0], v[1], v[4]),   // Side face 1
                    Face::triangle(v[1], v[2], v[4]),   // Side face 2
                    Face::triangle(v[2], v[3], v[4]),   // Side face 3
                    Face::triangle(v[3], v[0], v[4]),   // Side face 4
                ]
            }
            // 2D and 1D elements don't have 3D faces
            _ => Vec::new(),
        }
    }

    /// Build topology information
    pub fn build_topology(&mut self) {
        self.topology = MeshTopology::default();

        // Build vertex to cell connectivity
        for (cell_id, cell) in self.cells.iter().enumerate() {
            for &vertex_id in &cell.vertices {
                self.topology
                    .vertex_cells
                    .entry(vertex_id)
                    .or_default()
                    .push(cell_id);
            }
        }

        // Build edge to face connectivity
        for (face_id, face) in self.faces.iter().enumerate() {
            let n = face.vertices.len();
            for i in 0..n {
                let edge = Edge::new(face.vertices[i], face.vertices[(i + 1) % n]);
                self.topology
                    .edge_faces
                    .entry(edge)
                    .or_default()
                    .push(face_id);
            }
        }

        // Build face neighbors (faces that share an edge)
        // This is crucial for finite volume methods
        for (edge, face_ids) in &self.topology.edge_faces {
            if face_ids.len() == 2 {
                // Interior edge shared by exactly two faces
                let face1 = face_ids[0];
                let face2 = face_ids[1];

                self.topology
                    .face_neighbors
                    .entry(face1)
                    .or_default()
                    .push(face2);

                self.topology
                    .face_neighbors
                    .entry(face2)
                    .or_default()
                    .push(face1);
            }
        }
    }
}

impl<T: RealField + Copy> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vertex_creation() {
        let v = Vertex::from_coords(1.0, 2.0, 3.0);
        assert_relative_eq!(v.position.x, 1.0);
        assert_relative_eq!(v.position.y, 2.0);
        assert_relative_eq!(v.position.z, 3.0);
    }

    #[test]
    fn test_vertex_distance() {
        let v1 = Vertex::from_coords(0.0, 0.0, 0.0);
        let v2 = Vertex::from_coords(3.0, 4.0, 0.0);
        assert_relative_eq!(v1.distance_to(&v2), 5.0);
    }

    #[test]
    fn test_edge_creation() {
        // Test canonical ordering
        let edge1 = Edge::new(0, 1);
        assert_eq!(edge1.start, 0);
        assert_eq!(edge1.end, 1);

        // Test reverse ordering produces same edge
        let edge2 = Edge::new(1, 0);
        assert_eq!(edge2.start, 0);
        assert_eq!(edge2.end, 1);
        assert_eq!(edge1, edge2);

        assert!(edge1.contains(0));
        assert!(edge1.contains(1));
        assert!(!edge1.contains(2));
    }

    #[test]
    fn test_face_creation() {
        let tri = Face::triangle(0, 1, 2);
        assert_eq!(tri.vertex_count(), 3);

        let quad = Face::quad(0, 1, 2, 3);
        assert_eq!(quad.vertex_count(), 4);
    }

    #[test]
    fn test_mesh_operations() {
        let mut mesh = Mesh::<f64>::new();

        // Add vertices
        let v0 = mesh.add_vertex(Vertex::new(0, 0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex(Vertex::new(1, 1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex(Vertex::new(2, 0.0, 1.0, 0.0));

        assert_eq!(mesh.vertex_count(), 3);

        // Add edge
        mesh.add_edge(Edge::new(v0, v1));
        assert_eq!(mesh.edge_count(), 1);

        // Add face
        mesh.add_face(Face::triangle(v0, v1, v2));
        assert_eq!(mesh.face_count(), 1);
    }

    #[test]
    fn test_mesh_topology() {
        let mut mesh = Mesh::<f64>::new();

        // Create a simple tetrahedron
        mesh.add_vertex(Vertex::new(0, 0.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::new(1, 1.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::new(2, 0.0, 1.0, 0.0));
        mesh.add_vertex(Vertex::new(3, 0.0, 0.0, 1.0));

        // Add faces
        mesh.add_face(Face::triangle(0, 1, 2));
        mesh.add_face(Face::triangle(0, 1, 3));
        mesh.add_face(Face::triangle(0, 2, 3));
        mesh.add_face(Face::triangle(1, 2, 3));

        // Add cell
        mesh.cells.push(Cell {
            element_type: ElementType::Tetrahedron,
            vertices: vec![0, 1, 2, 3],
        });

        // Build topology
        mesh.build_topology();

        // Check vertex-cell connectivity
        assert_eq!(mesh.topology.vertex_cells[&0].len(), 1);
        assert_eq!(mesh.topology.vertex_cells[&1].len(), 1);
    }

    #[test]
    fn test_element_types() {
        let tet = ElementType::Tetrahedron;
        let hex = ElementType::Hexahedron;
        assert_ne!(tet, hex);
    }
}
