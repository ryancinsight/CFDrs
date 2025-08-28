//! Mesh data structures and operations

use super::element::ElementType;
use crate::Result;
use nalgebra::{Point3, RealField, Vector3};
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

/// Mesh statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshStatistics {
    /// Number of nodes
    pub num_nodes: usize,
    /// Number of elements
    pub num_elements: usize,
    /// Number of boundary nodes
    pub num_boundary_nodes: usize,
    /// Element type distribution
    pub element_distribution: HashMap<ElementType, usize>,
}

// Constants for mesh operations
const DIMENSION_TOLERANCE: f64 = 1e-10;
const MIN_ELEMENT_VOLUME: f64 = 1e-12;
const MAX_ASPECT_RATIO: f64 = 100.0;

impl<T: RealField + Copy> Mesh<T> {
    /// Create a new mesh
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
    pub fn add_node(&mut self, node: Point3<T>) -> usize {
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Add an element to the mesh
    pub fn add_element(&mut self, element: Element) -> Result<usize> {
        // Validate element nodes exist
        for &node_idx in &element.nodes {
            if node_idx >= self.nodes.len() {
                return Err(crate::error::Error::InvalidInput(format!(
                    "Node index {} out of bounds",
                    node_idx
                )));
            }
        }

        // Validate element has correct number of nodes
        if element.nodes.len() != element.element_type.num_nodes() {
            return Err(crate::error::Error::InvalidInput(format!(
                "Element type {:?} requires {} nodes, got {}",
                element.element_type,
                element.element_type.num_nodes(),
                element.nodes.len()
            )));
        }

        self.elements.push(element);
        Ok(self.elements.len() - 1)
    }

    /// Mark a node as boundary
    pub fn mark_boundary(&mut self, node_idx: usize, marker: String) -> Result<()> {
        if node_idx >= self.nodes.len() {
            return Err(crate::error::Error::InvalidInput(format!(
                "Node index {} out of bounds",
                node_idx
            )));
        }
        self.boundary_markers.insert(node_idx, marker);
        Ok(())
    }

    /// Get mesh statistics
    pub fn statistics(&self) -> MeshStatistics {
        let mut element_distribution = HashMap::new();
        for element in &self.elements {
            *element_distribution
                .entry(element.element_type)
                .or_insert(0) += 1;
        }

        MeshStatistics {
            num_nodes: self.nodes.len(),
            num_elements: self.elements.len(),
            num_boundary_nodes: self.boundary_markers.len(),
            element_distribution,
        }
    }

    /// Get element center
    pub fn element_center(&self, element_idx: usize) -> Result<Point3<T>> {
        if element_idx >= self.elements.len() {
            return Err(crate::error::Error::InvalidInput(format!(
                "Element index {} out of bounds",
                element_idx
            )));
        }

        let element = &self.elements[element_idx];
        let mut center = Vector3::zeros();
        let num_nodes = T::from_usize(element.nodes.len()).unwrap_or_else(T::one);

        for &node_idx in &element.nodes {
            center += self.nodes[node_idx].coords;
        }

        Ok(Point3::from(center / num_nodes))
    }

    /// Get element volume/area/length
    pub fn element_measure(&self, element_idx: usize) -> Result<T> {
        if element_idx >= self.elements.len() {
            return Err(crate::error::Error::InvalidInput(format!(
                "Element index {} out of bounds",
                element_idx
            )));
        }

        let element = &self.elements[element_idx];

        match element.element_type {
            ElementType::Line | ElementType::Line3 => {
                // Line length
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                Ok((p1 - p0).norm())
            }
            ElementType::Triangle => {
                // Triangle area using cross product
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                Ok(v1.cross(&v2).norm() / T::from_f64(2.0).unwrap_or_else(T::one))
            }
            ElementType::Tetrahedron => {
                // Tetrahedron volume using scalar triple product
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let p3 = &self.nodes[element.nodes[3]];
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let v3 = p3 - p0;
                Ok(v1.dot(&v2.cross(&v3)).abs() / T::from_f64(6.0).unwrap_or_else(T::one))
            }
            _ => {
                // For other element types, return a placeholder
                // This should be implemented properly for each element type
                Ok(T::one())
            }
        }
    }

    /// Validate mesh consistency
    pub fn validate(&self) -> Result<()> {
        // Check for empty mesh
        if self.nodes.is_empty() {
            return Err(crate::error::Error::InvalidInput(
                "Mesh has no nodes".to_string(),
            ));
        }

        // Check element validity
        for (idx, element) in self.elements.iter().enumerate() {
            // Check node indices
            for &node_idx in &element.nodes {
                if node_idx >= self.nodes.len() {
                    return Err(crate::error::Error::InvalidInput(format!(
                        "Element {} has invalid node index {}",
                        idx, node_idx
                    )));
                }
            }

            // Check element has correct number of nodes
            if element.nodes.len() != element.element_type.num_nodes() {
                return Err(crate::error::Error::InvalidInput(format!(
                    "Element {} has wrong number of nodes",
                    idx
                )));
            }

            // Check for degenerate elements
            let measure = self.element_measure(idx)?;
            if measure < T::from_f64(MIN_ELEMENT_VOLUME).unwrap_or_else(T::zero) {
                return Err(crate::error::Error::InvalidInput(format!(
                    "Element {} is degenerate (measure = {:?})",
                    idx, measure
                )));
            }
        }

        // Check boundary markers
        for &node_idx in self.boundary_markers.keys() {
            if node_idx >= self.nodes.len() {
                return Err(crate::error::Error::InvalidInput(format!(
                    "Boundary marker for non-existent node {}",
                    node_idx
                )));
            }
        }

        Ok(())
    }
}
