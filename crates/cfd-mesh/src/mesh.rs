//! Core mesh data structure with topology and connectivity

use crate::topology::{Cell, Edge, ElementType, Face, Vertex};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Main mesh data structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh<T: RealField + Copy> {
    /// Vertices in the mesh
    vertices: Vec<Vertex<T>>,
    /// Edges connecting vertices
    edges: Vec<Edge>,
    /// Faces formed by edges
    faces: Vec<Face>,
    /// Cells formed by faces
    cells: Vec<Cell>,
    /// Boundary markers for faces
    boundary_markers: HashMap<usize, String>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Create an empty mesh
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            cells: Vec::new(),
            boundary_markers: HashMap::new(),
        }
    }

    /// Add a vertex to the mesh
    pub fn add_vertex(&mut self, vertex: Vertex<T>) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(vertex);
        idx
    }

    /// Add an edge to the mesh
    pub fn add_edge(&mut self, edge: Edge) -> usize {
        let idx = self.edges.len();
        self.edges.push(edge);
        idx
    }

    /// Add a face to the mesh
    pub fn add_face(&mut self, face: Face) -> usize {
        let idx = self.faces.len();
        self.faces.push(face);
        idx
    }

    /// Add a cell to the mesh
    pub fn add_cell(&mut self, cell: Cell) -> usize {
        let idx = self.cells.len();
        self.cells.push(cell);
        idx
    }

    /// Mark a face as boundary with a label
    pub fn mark_boundary(&mut self, face_idx: usize, label: String) {
        self.boundary_markers.insert(face_idx, label);
    }

    /// Get vertex by index
    pub fn vertex(&self, idx: usize) -> Option<&Vertex<T>> {
        self.vertices.get(idx)
    }

    /// Get edge by index
    pub fn edge(&self, idx: usize) -> Option<&Edge> {
        self.edges.get(idx)
    }

    /// Get face by index
    pub fn face(&self, idx: usize) -> Option<&Face> {
        self.faces.get(idx)
    }

    /// Get cell by index
    pub fn cell(&self, idx: usize) -> Option<&Cell> {
        self.cells.get(idx)
    }

    /// Number of vertices
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of edges
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Number of faces
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Number of cells
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Get all boundary face indices
    pub fn boundary_faces(&self) -> Vec<usize> {
        self.boundary_markers.keys().copied().collect()
    }

    /// Get boundary label for a face
    pub fn boundary_label(&self, face_idx: usize) -> Option<&str> {
        self.boundary_markers.get(&face_idx).map(String::as_str)
    }

    /// Get all cells
    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    /// Get all vertices  
    pub fn vertices(&self) -> &[Vertex<T>] {
        &self.vertices
    }

    /// Get all edges
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Get all faces
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    /// Get faces of a cell as an iterator
    pub fn element_faces<'a>(&'a self, cell: &'a Cell) -> impl Iterator<Item = &'a Face> + 'a {
        cell.faces
            .iter()
            .filter_map(move |&face_idx| self.face(face_idx))
    }

    /// Get faces of a cell (allocating version for compatibility)
    #[deprecated(note = "Use element_faces() iterator for zero-copy access")]
    pub fn get_element_faces(&self, cell: &Cell) -> Vec<&Face> {
        self.element_faces(cell).collect()
    }

    /// Get vertices of a cell as an iterator (zero-copy)
    pub fn element_vertices<'a>(
        &'a self,
        cell: &'a Cell,
    ) -> impl Iterator<Item = &'a Vertex<T>> + 'a {
        use std::collections::HashSet;

        // Collect unique indices first (necessary for deduplication)
        let mut vertex_indices = HashSet::new();
        for &face_idx in &cell.faces {
            if let Some(face) = self.face(face_idx) {
                vertex_indices.extend(&face.vertices);
            }
        }

        // Return iterator over vertices
        vertex_indices
            .into_iter()
            .filter_map(move |idx| self.vertex(idx))
    }

    /// Get vertices of a cell (allocating version for compatibility)
    #[deprecated(note = "Use element_vertices() iterator for zero-copy access")]
    pub fn get_element_vertices(&self, cell: &Cell) -> Vec<&Vertex<T>> {
        self.element_vertices(cell).collect()
    }

    /// Check mesh validity
    pub fn validate(&self) -> Result<(), String> {
        // Check vertex indices in edges
        for edge in &self.edges {
            if edge.start >= self.vertices.len() || edge.end >= self.vertices.len() {
                return Err(format!("Edge references invalid vertex: {:?}", edge));
            }
        }

        // Check vertex indices in faces
        for face in &self.faces {
            for &v in &face.vertices {
                if v >= self.vertices.len() {
                    return Err(format!("Face references invalid vertex: {}", v));
                }
            }
        }

        // Check face indices in cells
        for cell in &self.cells {
            for &f in &cell.faces {
                if f >= self.faces.len() {
                    return Err(format!("Cell references invalid face: {}", f));
                }
            }
        }

        Ok(())
    }

    /// Compute mesh statistics
    pub fn statistics(&self) -> MeshStatistics {
        let mut stats = MeshStatistics::default();
        stats.vertex_count = self.vertices.len();
        stats.edge_count = self.edges.len();
        stats.face_count = self.faces.len();
        stats.cell_count = self.cells.len();
        stats.boundary_face_count = self.boundary_markers.len();

        // Count element types
        for cell in &self.cells {
            match cell.element_type {
                ElementType::Tetrahedron => stats.tetrahedra += 1,
                ElementType::Hexahedron => stats.hexahedra += 1,
                ElementType::Pyramid => stats.pyramids += 1,
                ElementType::Prism => stats.prisms += 1,
                _ => {}
            }
        }

        stats
    }
}

impl<T: RealField + Copy> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Mesh statistics
#[derive(Debug, Default, Clone)]
pub struct MeshStatistics {
    /// Number of vertices
    pub vertex_count: usize,
    /// Number of edges
    pub edge_count: usize,
    /// Number of faces
    pub face_count: usize,
    /// Number of cells
    pub cell_count: usize,
    /// Number of boundary faces
    pub boundary_face_count: usize,
    /// Number of tetrahedral cells
    pub tetrahedra: usize,
    /// Number of hexahedral cells
    pub hexahedra: usize,
    /// Number of pyramid cells
    pub pyramids: usize,
    /// Number of prism cells
    pub prisms: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mesh_operations() {
        let mut mesh = Mesh::<f64>::new();

        // Add vertices
        let v0 = mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));

        assert_eq!(mesh.vertex_count(), 3);

        // Add edge
        mesh.add_edge(Edge::new(v0, v1));
        assert_eq!(mesh.edge_count(), 1);

        // Add face
        mesh.add_face(Face::triangle(v0, v1, v2));
        assert_eq!(mesh.face_count(), 1);

        // Validate mesh
        assert!(mesh.validate().is_ok());
    }

    #[test]
    fn test_mesh_topology() {
        let mut mesh = Mesh::<f64>::new();

        // Create a simple tetrahedron
        mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 1.0));

        // Add faces
        mesh.add_face(Face::triangle(0, 1, 2));
        mesh.add_face(Face::triangle(0, 1, 3));
        mesh.add_face(Face::triangle(0, 2, 3));
        mesh.add_face(Face::triangle(1, 2, 3));

        // Add cell
        mesh.add_cell(Cell::tetrahedron(0, 1, 2, 3));

        assert_eq!(mesh.vertex_count(), 4);
        assert_eq!(mesh.face_count(), 4);
        assert_eq!(mesh.cell_count(), 1);

        // Check statistics
        let stats = mesh.statistics();
        assert_eq!(stats.tetrahedra, 1);
        assert_eq!(stats.hexahedra, 0);
    }

    #[test]
    fn test_boundary_marking() {
        let mut mesh = Mesh::<f64>::new();

        // Add a simple face
        mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));
        let face_idx = mesh.add_face(Face::triangle(0, 1, 2));

        // Mark as boundary
        mesh.mark_boundary(face_idx, "inlet".to_string());

        assert_eq!(mesh.boundary_faces().len(), 1);
        assert_eq!(mesh.boundary_label(face_idx), Some("inlet"));
    }
}
