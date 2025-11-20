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
    /// Global ID for distributed meshes
    pub global_id: Option<usize>,
    /// Partition ID (rank) that owns this cell
    pub partition_id: Option<usize>,
}

impl Cell {
    /// Create a tetrahedral cell
    #[must_use]
    pub fn tetrahedron(f0: usize, f1: usize, f2: usize, f3: usize) -> Self {
        Self {
            faces: vec![f0, f1, f2, f3],
            element_type: ElementType::Tetrahedron,
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a hexahedral cell
    #[must_use]
    pub fn hexahedron(faces: Vec<usize>) -> Self {
        assert_eq!(faces.len(), 6, "Hexahedron must have 6 faces");
        Self {
            faces,
            element_type: ElementType::Hexahedron,
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a pyramidal cell
    #[must_use]
    pub fn pyramid(faces: Vec<usize>) -> Self {
        assert_eq!(faces.len(), 5, "Pyramid must have 5 faces");
        Self {
            faces,
            element_type: ElementType::Pyramid,
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a prism cell
    #[must_use]
    pub fn prism(faces: Vec<usize>) -> Self {
        assert_eq!(faces.len(), 5, "Prism must have 5 faces");
        Self {
            faces,
            element_type: ElementType::Prism,
            global_id: None,
            partition_id: None,
        }
    }

    /// Set distributed mesh properties
    pub fn with_distributed_info(mut self, global_id: usize, partition_id: usize) -> Self {
        self.global_id = Some(global_id);
        self.partition_id = Some(partition_id);
        self
    }

    /// Number of faces in the cell
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Check if cell contains a face
    #[must_use]
    pub fn contains_face(&self, face: usize) -> bool {
        self.faces.contains(&face)
    }
}
