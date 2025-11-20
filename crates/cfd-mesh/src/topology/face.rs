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

    /// Create a face from vertex indices
    #[must_use]
    pub fn from_vertices(vertices: Vec<usize>) -> Self {
        Self { vertices }
    }

    /// Number of vertices in the face
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Check if face contains a vertex
    #[must_use]
    pub fn contains(&self, vertex: usize) -> bool {
        self.vertices.contains(&vertex)
    }

    /// Get edges of the face
    #[must_use]
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let n = self.vertices.len();
        (0..n)
            .map(|i| (self.vertices[i], self.vertices[(i + 1) % n]))
            .collect()
    }

    /// Check if face is a triangle
    #[must_use]
    pub fn is_triangle(&self) -> bool {
        self.vertices.len() == 3
    }

    /// Check if face is a quadrilateral
    #[must_use]
    pub fn is_quad(&self) -> bool {
        self.vertices.len() == 4
    }
}
