//! Mesh quality assessment and metrics

use super::mesh::Mesh;
use crate::error::Result;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mesh quality assessment trait
pub trait MeshQuality<T: RealField + Copy>: Send + Sync {
    /// Compute aspect ratio for an element
    fn aspect_ratio(&self, mesh: &Mesh<T>, element_idx: usize) -> Result<T>;

    /// Compute skewness for an element
    fn skewness(&self, mesh: &Mesh<T>, element_idx: usize) -> Result<T>;

    /// Generate quality report for the mesh
    fn quality_report(&self, mesh: &Mesh<T>) -> Result<QualityReport<T>>;
}

/// Quality report for a mesh
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityReport<T: RealField + Copy> {
    /// Overall quality score (0-1, 1 being perfect)
    pub overall_score: T,
    /// Statistics for various quality metrics
    pub statistics: QualityStatistics<T>,
    /// Elements with quality issues
    pub problematic_elements: Vec<usize>,
    /// Recommendations for improvement
    pub recommendations: Vec<String>,
}

/// Quality statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityStatistics<T: RealField + Copy> {
    /// Minimum aspect ratio
    pub min_aspect_ratio: T,
    /// Maximum aspect ratio
    pub max_aspect_ratio: T,
    /// Average aspect ratio
    pub avg_aspect_ratio: T,
    /// Minimum skewness
    pub min_skewness: T,
    /// Maximum skewness
    pub max_skewness: T,
    /// Average skewness
    pub avg_skewness: T,
    /// Distribution of quality scores
    pub quality_distribution: HashMap<String, usize>,
}
