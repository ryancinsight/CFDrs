//! Mesh refinement and adaptation strategies
//!
//! This module provides tools for dynamically adjusting the mesh resolution
//! based on geometric or solution-based criteria.

use super::Mesh;
use crate::error::Result;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Trait for mesh refinement and adaptation strategies
pub trait MeshRefinement<T: RealField + Copy>: Send + Sync {
    /// Subdivide elements to increase resolution based on criteria
    fn refine(&self, mesh: &mut Mesh<T>, criteria: &RefinementCriteria<T>) -> Result<()>;

    /// Merge elements to decrease resolution based on criteria
    fn coarsen(&self, mesh: &mut Mesh<T>, criteria: &RefinementCriteria<T>) -> Result<()>;

    /// Adapt the mesh resolution based on solution gradients or error estimates
    fn adapt(&self, mesh: &mut Mesh<T>, solution: &[T]) -> Result<()>;
}

/// Criteria for controlling mesh refinement and coarsening
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefinementCriteria<T: RealField + Copy> {
    /// Maximum allowable element size (length scale)
    pub max_size: T,
    /// Minimum allowable element size (length scale)
    pub min_size: T,
    /// Optional target number of elements for the entire mesh
    pub target_elements: Option<usize>,
    /// Optional error threshold for solution-based adaptation
    pub error_threshold: Option<T>,
    /// Optional gradient threshold for refinement in high-gradient regions
    pub gradient_threshold: Option<T>,
}
