//! Cell representation for 3D volumetric elements

use super::ElementType;
use serde::{Deserialize, Serialize};

/// Cell defined by faces
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    /// Indices of faces forming the cell
    pub faces: Vec<usize>,
    /// Type of element
    pub element_type: ElementType,
}

impl Cell {
    /// Create a tetrahedral cell
    pub fn tetrahedron(f0: usize, f1: usize, f2: usize, f3: usize) -> Self {
        Self {
            faces: vec![f0, f1, f2, f3],
            element_type: ElementType::Tetrahedron,
        }
    }

    /// Create a hexahedral cell
    pub fn hexahedron(faces: Vec<usize>) -> Self {
        assert_eq!(faces.len(), 6, "Hexahedron must have 6 faces");
        Self {
            faces,
            element_type: ElementType::Hexahedron,
        }
    }

    /// Create a pyramidal cell
    pub fn pyramid(faces: Vec<usize>) -> Self {
        assert_eq!(faces.len(), 5, "Pyramid must have 5 faces");
        Self {
            faces,
            element_type: ElementType::Pyramid,
        }
    }

    /// Create a prism cell
    pub fn prism(faces: Vec<usize>) -> Self {
        assert_eq!(faces.len(), 5, "Prism must have 5 faces");
        Self {
            faces,
            element_type: ElementType::Prism,
        }
    }

    /// Number of faces in the cell
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Check if cell contains a face
    pub fn contains_face(&self, face: usize) -> bool {
        self.faces.contains(&face)
    }
}
