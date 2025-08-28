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
    /// Helper function to calculate tetrahedron volume
    fn tet_volume(p0: &Point3<T>, p1: &Point3<T>, p2: &Point3<T>, p3: &Point3<T>) -> T {
        let v1 = p1 - p0;
        let v2 = p2 - p0;
        let v3 = p3 - p0;
        v1.dot(&v2.cross(&v3)).abs() / T::from_f64(6.0).unwrap_or_else(T::one)
    }

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
            ElementType::Quadrilateral => {
                // Quadrilateral area using shoelace formula
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let p3 = &self.nodes[element.nodes[3]];
                // Split into two triangles and sum areas
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let area1 = v1.cross(&v2).norm() / T::from_f64(2.0).unwrap_or_else(T::one);
                let v3 = p3 - p0;
                let area2 = v2.cross(&v3).norm() / T::from_f64(2.0).unwrap_or_else(T::one);
                Ok(area1 + area2)
            }
            ElementType::Hexahedron => {
                // Hexahedron volume using decomposition into tetrahedra
                let nodes: Vec<_> = element.nodes.iter().map(|&i| &self.nodes[i]).collect();
                // Decompose into 6 tetrahedra and sum volumes
                let mut volume = T::zero();
                // Tetrahedron 1: 0-1-3-4
                volume += Self::tet_volume(nodes[0], nodes[1], nodes[3], nodes[4]);
                // Tetrahedron 2: 1-2-3-6
                volume += Self::tet_volume(nodes[1], nodes[2], nodes[3], nodes[6]);
                // Tetrahedron 3: 1-3-4-6
                volume += Self::tet_volume(nodes[1], nodes[3], nodes[4], nodes[6]);
                // Tetrahedron 4: 3-4-6-7
                volume += Self::tet_volume(nodes[3], nodes[4], nodes[6], nodes[7]);
                // Tetrahedron 5: 1-4-5-6
                volume += Self::tet_volume(nodes[1], nodes[4], nodes[5], nodes[6]);
                // Tetrahedron 6: 4-5-6-7
                volume += Self::tet_volume(nodes[4], nodes[5], nodes[6], nodes[7]);
                Ok(volume)
            }
            ElementType::Pyramid => {
                // Pyramid volume = (1/3) * base_area * height
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let p3 = &self.nodes[element.nodes[3]];
                let apex = &self.nodes[element.nodes[4]];

                // Calculate base area (quadrilateral)
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let area1 = v1.cross(&v2).norm() / T::from_f64(2.0).unwrap_or_else(T::one);
                let v3 = p3 - p0;
                let area2 = v2.cross(&v3).norm() / T::from_f64(2.0).unwrap_or_else(T::one);
                let base_area = area1 + area2;

                // Calculate height (perpendicular distance from apex to base)
                let base_normal = v1.cross(&v2).normalize();
                let height_vec = apex - p0;
                let height = height_vec.dot(&base_normal).abs();

                Ok(base_area * height / T::from_f64(3.0).unwrap_or_else(T::one))
            }
            ElementType::Prism => {
                // Prism/Wedge volume = base_area * height
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let p3 = &self.nodes[element.nodes[3]];

                // Calculate triangular base area
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let base_area = v1.cross(&v2).norm() / T::from_f64(2.0).unwrap_or_else(T::one);

                // Calculate height (distance between triangular faces)
                let height_vec = p3 - p0;
                let base_normal = v1.cross(&v2).normalize();
                let height = height_vec.dot(&base_normal).abs();

                Ok(base_area * height)
            }
            ElementType::Line3 => {
                // Use first two nodes for length
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                Ok((p1 - p0).norm())
            }
            ElementType::Triangle6 => {
                // Use first three nodes for area
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                Ok(v1.cross(&v2).norm() / T::from_f64(2.0).unwrap_or_else(T::one))
            }
            ElementType::Quadrilateral9 => {
                // Use first four nodes for area
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let p3 = &self.nodes[element.nodes[3]];
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let area1 = v1.cross(&v2).norm() / T::from_f64(2.0).unwrap_or_else(T::one);
                let v3 = p3 - p0;
                let area2 = v2.cross(&v3).norm() / T::from_f64(2.0).unwrap_or_else(T::one);
                Ok(area1 + area2)
            }
            ElementType::Tetrahedron10 => {
                // Use first four nodes for volume
                let p0 = &self.nodes[element.nodes[0]];
                let p1 = &self.nodes[element.nodes[1]];
                let p2 = &self.nodes[element.nodes[2]];
                let p3 = &self.nodes[element.nodes[3]];
                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let v3 = p3 - p0;
                Ok(v1.dot(&v2.cross(&v3)).abs() / T::from_f64(6.0).unwrap_or_else(T::one))
            }
            ElementType::Hexahedron20 => {
                // Use first eight nodes for volume
                let nodes: Vec<_> = element.nodes[..8].iter().map(|&i| &self.nodes[i]).collect();
                let mut volume = T::zero();
                volume += Self::tet_volume(nodes[0], nodes[1], nodes[3], nodes[4]);
                volume += Self::tet_volume(nodes[1], nodes[2], nodes[3], nodes[6]);
                volume += Self::tet_volume(nodes[1], nodes[3], nodes[4], nodes[6]);
                volume += Self::tet_volume(nodes[3], nodes[4], nodes[6], nodes[7]);
                volume += Self::tet_volume(nodes[1], nodes[4], nodes[5], nodes[6]);
                volume += Self::tet_volume(nodes[4], nodes[5], nodes[6], nodes[7]);
                Ok(volume)
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
