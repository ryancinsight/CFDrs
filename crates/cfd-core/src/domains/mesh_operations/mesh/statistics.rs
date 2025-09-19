//! Mesh statistics and analysis

use super::types::Mesh;
use crate::domains::mesh_operations::element::ElementType;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mesh statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshStatistics {
    /// Number of nodes
    pub num_nodes: usize,
    /// Number of elements
    pub num_elements: usize,
    /// Element type distribution
    pub element_types: HashMap<ElementType, usize>,
    /// Boundary nodes count
    pub num_boundary_nodes: usize,
    /// Material distribution
    pub material_distribution: HashMap<usize, usize>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Compute mesh statistics
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
