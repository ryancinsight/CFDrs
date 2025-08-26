//! Core mesh data structures and operations
//! Following SOLID principles

use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Import ElementType from cfd-core as the single source of truth
pub use cfd_core::domains::mesh_operations::ElementType;

/// Vertex in 3D space
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vertex<T: RealField + Copy> {
    /// Position in 3D space
    pub position: Point3<T>,
    /// Vertex ID
    pub id: usize,
}

impl<T: RealField + Copy> Vertex<T> {
    /// Create a new vertex
    pub fn new(id: usize, x: T, y: T, z: T) -> Self {
        Self {
            position: Point3::new(x, y, z),
            id,
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
    /// Create a new edge
    #[must_use]
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
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
    /// Face ID
    pub id: usize,
    /// Vertex indices
    pub vertices: Vec<usize>,
}

impl Face {
    /// Create a triangular face
    #[must_use]
    pub fn triangle(v0: usize, v1: usize, v2: usize) -> Self {
        Self {
            id: 0, // ID should be set when adding to mesh
            vertices: vec![v0, v1, v2],
        }
    }

    /// Create a quadrilateral face
    #[must_use]
    pub fn quad(v0: usize, v1: usize, v2: usize, v3: usize) -> Self {
        Self {
            id: 0, // ID should be set when adding to mesh
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
    pub fn get_element_faces(&self, _element: &Cell) -> Vec<&Face> {
        // For now, return empty vector
        // Full implementation would compute faces from element vertices
        Vec::new()
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
        let v = Vertex::new(0, 1.0, 2.0, 3.0);
        assert_eq!(v.id, 0);
        assert_relative_eq!(v.position.x, 1.0);
        assert_relative_eq!(v.position.y, 2.0);
        assert_relative_eq!(v.position.z, 3.0);
    }

    #[test]
    fn test_vertex_distance() {
        let v1 = Vertex::new(0, 0.0, 0.0, 0.0);
        let v2 = Vertex::new(1, 3.0, 4.0, 0.0);
        assert_relative_eq!(v1.distance_to(&v2), 5.0);
    }

    #[test]
    fn test_edge_creation() {
        let edge = Edge::new(0, 1);
        assert_eq!(edge.start, 0);
        assert_eq!(edge.end, 1);
        assert!(edge.contains(0));
        assert!(edge.contains(1));
        assert!(!edge.contains(2));
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
