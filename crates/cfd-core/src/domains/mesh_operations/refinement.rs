//! Mesh refinement operations and criteria

use super::mesh::Mesh;
use crate::Result;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Mesh refinement trait
pub trait MeshRefinement<T: RealField + Copy>: Send + Sync {
    /// Refine mesh based on criteria
    fn refine(&self, mesh: &mut Mesh<T>, criteria: &RefinementCriteria<T>) -> Result<()>;

    /// Coarsen mesh based on criteria
    fn coarsen(&self, mesh: &mut Mesh<T>, criteria: &RefinementCriteria<T>) -> Result<()>;

    /// Adapt mesh based on solution
    fn adapt(&self, mesh: &mut Mesh<T>, solution: &[T]) -> Result<()>;
}

/// Refinement criteria
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefinementCriteria<T: RealField + Copy> {
    /// Maximum element size
    pub max_size: T,
    /// Minimum element size
    pub min_size: T,
    /// Target number of elements
    pub target_elements: Option<usize>,
    /// Error threshold for adaptation
    pub error_threshold: Option<T>,
    /// Gradient threshold for refinement
    pub gradient_threshold: Option<T>,
}
