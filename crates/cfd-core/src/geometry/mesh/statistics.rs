//! Mesh statistics and analysis
//!
//! This module provides tools for gathering quantitative information about
//! the mesh, such as element counts and material distributions.

use super::{ElementType, Mesh};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Statistical summary of a mesh
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshStatistics {
    /// Total number of nodes
    pub num_nodes: usize,
    /// Total number of elements
    pub num_elements: usize,
    /// Distribution of elements by type
    pub element_types: HashMap<ElementType, usize>,
    /// Number of nodes with boundary markers
    pub num_boundary_nodes: usize,
    /// Distribution of elements by material ID
    pub material_distribution: HashMap<usize, usize>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Compute a statistical summary of the mesh
    #[must_use]
    pub fn compute_statistics(&self) -> MeshStatistics {
        let mut element_types = HashMap::new();
        let mut material_distribution = HashMap::new();

        for element in &self.elements {
            *element_types.entry(element.element_type).or_insert(0) += 1;
            *material_distribution
                .entry(element.material_id)
                .or_insert(0) += 1;
        }

        MeshStatistics {
            num_nodes: self.nodes.len(),
            num_elements: self.elements.len(),
            element_types,
            num_boundary_nodes: self.boundary_markers.len(),
            material_distribution,
        }
    }
}
