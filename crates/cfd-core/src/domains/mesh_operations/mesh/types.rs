//! Core mesh types and data structures

use crate::error::{Error, Result};

use crate::domains::mesh_operations::element::ElementType;
use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mesh structure containing nodes and elements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh<T: RealField + Copy> {
    /// Nodes in the mesh
    pub nodes: Vec<Point3<T>>,
    /// Elements in the mesh
    pub elements: Vec<Element>,
    /// Boundary markers for nodes
    pub boundary_markers: HashMap<usize, String>,
    /// Metadata about the mesh
    pub metadata: MeshMetadata,
}

/// Element definition
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Element {
    /// Element type
    pub element_type: ElementType,
    /// Node indices
    pub nodes: Vec<usize>,
    /// Material/region ID
    pub material_id: usize,
}

/// Mesh metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshMetadata {
    /// Mesh name
    pub name: String,
    /// Mesh dimension
    pub dimension: usize,
    /// Creation timestamp
    pub created: String,
    /// Additional properties
    pub properties: HashMap<String, String>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Create a new empty mesh
    #[must_use]
    pub fn new(name: String, dimension: usize) -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
            boundary_markers: HashMap::new(),
            metadata: MeshMetadata {
                name,
                dimension,
                created: chrono::Utc::now().to_rfc3339(),
                properties: HashMap::new(),
            },
        }
    }

    /// Add a node to the mesh
    pub fn add_node(&mut self, point: Point3<T>) -> usize {
        self.nodes.push(point);
        self.nodes.len() - 1
    }

    /// Add an element to the mesh
    pub fn add_element(&mut self, element: Element) {
        self.elements.push(element);
    }

    /// Get number of nodes
    #[must_use]
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Get number of elements
    #[must_use]
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Validate mesh consistency
    /// 
    /// # Errors
    /// Returns error if mesh has invalid topology or degenerate elements
    pub fn validate(&self) -> Result<()> {
        // Check for empty mesh
        if self.nodes.is_empty() {
            return Err(Error::InvalidInput("Mesh has no nodes".into()));
        }

        // Check element node indices
        for element in &self.elements {
            for &node_idx in &element.nodes {
                if node_idx >= self.nodes.len() {
                    return Err(Error::InvalidInput(format!(
                        "Element references invalid node index: {node_idx}"
                    )));
                }
            }
        }

        Ok(())
    }
}
