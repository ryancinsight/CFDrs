//! Face representation for 2D elements in 3D meshes

use serde::{Deserialize, Serialize};

/// Face defined by vertex indices
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Face {
    /// Ordered list of vertex indices
    pub vertices: Vec<usize>,
}

impl Face {
    /// Create a triangular face
    pub fn triangle(v0: usize, v1: usize, v2: usize) -> Self {
        Self {
            vertices: vec![v0, v1, v2],
        }
    }

    /// Create a quadrilateral face
    pub fn quad(v0: usize, v1: usize, v2: usize, v3: usize) -> Self {
        Self {
            vertices: vec![v0, v1, v2, v3],
        }
    }

    /// Create a face from vertex indices
    pub fn from_vertices(vertices: Vec<usize>) -> Self {
        Self { vertices }
    }

    /// Number of vertices in the face
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Check if face contains a vertex
    pub fn contains(&self, vertex: usize) -> bool {
        self.vertices.contains(&vertex)
    }

    /// Get edges of the face
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let n = self.vertices.len();
        (0..n)
            .map(|i| (self.vertices[i], self.vertices[(i + 1) % n]))
            .collect()
    }

    /// Check if face is a triangle
    pub fn is_triangle(&self) -> bool {
        self.vertices.len() == 3
    }

    /// Check if face is a quadrilateral
    pub fn is_quad(&self) -> bool {
        self.vertices.len() == 4
    }
}
