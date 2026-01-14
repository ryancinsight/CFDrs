//! Mesh representation and operations for computational domains
//!
//! This module provides the data structures and traits for managing
//! discretized domains (meshes), including nodes, elements, and connectivity.

use crate::error::{Error, Result};
use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Element types supported by the mesh system
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElementType {
    // 1D Elements
    /// Line element (2 nodes)
    Line,
    /// Quadratic line element (3 nodes)
    Line3,

    // 2D Elements
    /// Triangle element (3 nodes)
    Triangle,
    /// Quadratic triangle element (6 nodes)
    Triangle6,
    /// Quadrilateral element (4 nodes)
    Quadrilateral,
    /// Quadratic quadrilateral element (9 nodes)
    Quadrilateral9,

    // 3D Elements
    /// Tetrahedral element (4 nodes)
    Tetrahedron,
    /// Quadratic tetrahedral element (10 nodes)
    Tetrahedron10,
    /// Hexahedral element (8 nodes)
    Hexahedron,
    /// Quadratic hexahedral element (20 nodes)
    Hexahedron20,
    /// Pyramid element (5 nodes)
    Pyramid,
    /// Prism/Wedge element (6 nodes)
    Prism,
}

impl ElementType {
    /// Get the number of nodes for this element type
    #[must_use]
    pub fn num_nodes(&self) -> usize {
        match self {
            Self::Line => 2,
            Self::Line3 | Self::Triangle => 3,
            Self::Triangle6 | Self::Prism => 6,
            Self::Quadrilateral | Self::Tetrahedron => 4,
            Self::Quadrilateral9 => 9,
            Self::Tetrahedron10 => 10,
            Self::Hexahedron => 8,
            Self::Hexahedron20 => 20,
            Self::Pyramid => 5,
        }
    }

    /// Get the dimension of this element type
    #[must_use]
    pub fn dimension(&self) -> usize {
        match self {
            Self::Line | Self::Line3 => 1,
            Self::Triangle | Self::Triangle6 | Self::Quadrilateral | Self::Quadrilateral9 => 2,
            Self::Tetrahedron
            | Self::Tetrahedron10
            | Self::Hexahedron
            | Self::Hexahedron20
            | Self::Pyramid
            | Self::Prism => 3,
        }
    }
}

/// Element definition in the mesh
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Element {
    /// Element type
    pub element_type: ElementType,
    /// Indices of nodes forming this element
    pub nodes: Vec<usize>,
    /// Material or region ID associated with this element
    pub material_id: usize,
}

/// Core mesh data structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh<T: RealField + Copy> {
    /// List of nodal coordinates
    pub nodes: Vec<Point3<T>>,
    /// List of elements
    pub elements: Vec<Element>,
    /// Boundary markers mapping node indices to boundary names
    pub boundary_markers: HashMap<usize, String>,
    /// Mesh metadata
    pub metadata: MeshMetadata,
}

/// Metadata for mesh identification and tracking
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshMetadata {
    /// Descriptive name
    pub name: String,
    /// Spatial dimension
    pub dimension: usize,
    /// ISO 8601 creation timestamp
    pub created: String,
    /// Arbitrary key-value properties
    pub properties: HashMap<String, String>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Create a new empty mesh
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

    /// Add a node and return its index
    pub fn add_node(&mut self, point: Point3<T>) -> usize {
        self.nodes.push(point);
        self.nodes.len() - 1
    }

    /// Add an element to the mesh
    pub fn add_element(
        &mut self,
        element_type: ElementType,
        nodes: Vec<usize>,
        material_id: usize,
    ) {
        self.elements.push(Element {
            element_type,
            nodes,
            material_id,
        });
    }

    /// Validate mesh topology and consistency
    pub fn validate(&self) -> Result<()> {
        if self.nodes.is_empty() {
            return Err(Error::InvalidInput("Mesh has no nodes".into()));
        }

        for (i, element) in self.elements.iter().enumerate() {
            if element.nodes.len() != element.element_type.num_nodes() {
                return Err(Error::InvalidInput(format!(
                    "Element {i} has incorrect number of nodes"
                )));
            }
            for &node_idx in &element.nodes {
                if node_idx >= self.nodes.len() {
                    return Err(Error::InvalidInput(format!(
                        "Element {i} references non-existent node {node_idx}"
                    )));
                }
            }
        }
        Ok(())
    }
}

/// Trait for mesh generation strategies
pub trait MeshGeneration<T: RealField + Copy>: Send + Sync {
    /// Generate a structured grid
    fn generate_structured(&self, nx: usize, ny: usize, nz: usize) -> Result<Mesh<T>>;

    /// Generate an unstructured mesh with target element size
    fn generate_unstructured(&self, target_size: T) -> Result<Mesh<T>>;
}

pub mod connectivity;
pub mod operations;
pub mod quality;
pub mod refinement;
pub mod service;
pub mod statistics;

pub use connectivity::Connectivity;
pub use operations::MeshOperations;
pub use quality::{
    MeshQuality, MeshQualityService, MetricStatistics, QualityAssessment, QualityLevel,
    QualityReport, QualityStatistics,
};
pub use refinement::{MeshRefinement, RefinementCriteria};
pub use service::MeshOperationsService;
pub use statistics::MeshStatistics;
