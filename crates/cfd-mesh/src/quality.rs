//! Mesh quality analysis and metrics.
//!
//! This module provides tools for analyzing mesh quality using various metrics
//! such as aspect ratio, skewness, and orthogonality.

use crate::mesh::{Mesh, Cell};
use nalgebra::{RealField, Point3};
use serde::{Deserialize, Serialize};

/// Mesh quality metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics<T: RealField> {
    /// Aspect ratio statistics
    pub aspect_ratio: QualityStatistics<T>,
    /// Skewness statistics
    pub skewness: QualityStatistics<T>,
    /// Orthogonality statistics
    pub orthogonality: QualityStatistics<T>,
    /// Volume statistics
    pub volume: QualityStatistics<T>,
    /// Number of cells analyzed
    pub num_cells: usize,
}

/// Statistical summary of a quality metric
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityStatistics<T: RealField> {
    /// Minimum value
    pub min: T,
    /// Maximum value
    pub max: T,
    /// Mean value
    pub mean: T,
    /// Standard deviation
    pub std_dev: T,
    /// Number of samples
    pub count: usize,
}

/// Mesh quality analyzer
#[derive(Debug)]
pub struct QualityAnalyzer<T: RealField> {
    /// Minimum acceptable aspect ratio
    pub min_aspect_ratio: T,
    /// Maximum acceptable skewness
    pub max_skewness: T,
    /// Minimum acceptable orthogonality
    pub min_orthogonality: T,
}

impl<T: RealField> Default for QualityAnalyzer<T>
where
    T: From<f64>,
{
    fn default() -> Self {
        Self {
            min_aspect_ratio: T::from(0.1),
            max_skewness: T::from(0.8),
            min_orthogonality: T::from(0.2),
        }
    }
}

impl<T: RealField + std::iter::Sum + Copy> QualityAnalyzer<T>
where
    T: From<f64> + From<usize>,
{
    /// Create a new quality analyzer with custom thresholds
    pub fn new(min_aspect_ratio: T, max_skewness: T, min_orthogonality: T) -> Self {
        Self {
            min_aspect_ratio,
            max_skewness,
            min_orthogonality,
        }
    }

    /// Analyze mesh quality using iterator combinators for efficiency
    pub fn analyze(&self, mesh: &Mesh<T>) -> QualityMetrics<T> {
        // Use iterator patterns for zero-copy quality analysis
        let quality_data: Vec<_> = mesh.cells
            .iter()
            .filter_map(|cell| self.analyze_cell(cell, mesh))
            .collect();

        if quality_data.is_empty() {
            return self.empty_metrics();
        }

        // Extract individual metrics using iterator combinators
        let aspect_ratios: Vec<T> = quality_data.iter().map(|(ar, _, _, _)| *ar).collect();
        let skewness_values: Vec<T> = quality_data.iter().map(|(_, sk, _, _)| *sk).collect();
        let orthogonality_values: Vec<T> = quality_data.iter().map(|(_, _, orth, _)| *orth).collect();
        let volumes: Vec<T> = quality_data.iter().map(|(_, _, _, vol)| *vol).collect();

        QualityMetrics {
            aspect_ratio: self.compute_statistics(&aspect_ratios),
            skewness: self.compute_statistics(&skewness_values),
            orthogonality: self.compute_statistics(&orthogonality_values),
            volume: self.compute_statistics(&volumes),
            num_cells: quality_data.len(),
        }
    }

    /// Analyze individual cell quality
    fn analyze_cell(&self, cell: &Cell, mesh: &Mesh<T>) -> Option<(T, T, T, T)> {
        // Get vertices from faces for this cell
        let vertices: Vec<_> = cell.faces
            .iter()
            .filter_map(|&face_idx| mesh.faces.get(face_idx))
            .flat_map(|face| &face.vertices)
            .filter_map(|&vertex_idx| mesh.vertices.get(vertex_idx))
            .collect();

        if vertices.len() < 3 {
            return None;
        }

        // Extract positions from vertices
        let positions: Vec<_> = vertices.iter().map(|v| &v.position).collect();

        // Calculate basic metrics
        let aspect_ratio = self.calculate_aspect_ratio(&positions);
        let skewness = self.calculate_skewness(&positions);
        let orthogonality = self.calculate_orthogonality(&positions);
        let volume = self.calculate_volume(&positions);

        Some((aspect_ratio, skewness, orthogonality, volume))
    }

    /// Calculate aspect ratio (simplified)
    fn calculate_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() < 2 {
            return T::one();
        }

        // Find min and max distances using iterator patterns
        let distances: Vec<T> = vertices
            .windows(2)
            .map(|pair| (pair[1] - pair[0]).norm())
            .collect();

        let min_dist = distances.iter().cloned().fold(T::infinity(), T::min);
        let max_dist = distances.iter().cloned().fold(T::neg_infinity(), T::max);

        if *max_dist > T::zero() {
            min_dist.clone() / max_dist.clone()
        } else {
            T::one()
        }
    }

    /// Calculate skewness (simplified)
    fn calculate_skewness(&self, _vertices: &[&Point3<T>]) -> T {
        // Simplified skewness calculation
        // In practice, this would be based on angle deviations
        T::from(0.1) // Placeholder
    }

    /// Calculate orthogonality (simplified)
    fn calculate_orthogonality(&self, _vertices: &[&Point3<T>]) -> T {
        // Simplified orthogonality calculation
        // In practice, this would be based on face normal angles
        T::from(0.9) // Placeholder
    }

    /// Calculate cell volume (simplified)
    fn calculate_volume(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() < 4 {
            return T::zero();
        }

        // Simplified volume calculation for tetrahedron
        let v1 = vertices[1] - vertices[0];
        let v2 = vertices[2] - vertices[0];
        let v3 = vertices[3] - vertices[0];

        let six = T::from(6.0);
        v1.cross(&v2).dot(&v3).abs() / six
    }

    /// Compute statistics for a set of values using iterator patterns
    fn compute_statistics(&self, values: &[T]) -> QualityStatistics<T> {
        if values.is_empty() {
            return QualityStatistics {
                min: T::zero(),
                max: T::zero(),
                mean: T::zero(),
                std_dev: T::zero(),
                count: 0,
            };
        }

        // Use iterator combinators for efficient statistics computation
        let min = values.iter().cloned().fold(T::infinity(), T::min);
        let max = values.iter().cloned().fold(T::neg_infinity(), T::max);
        let sum: T = values.iter().cloned().sum();
        let count = values.len();
        let mean = sum / T::from(count);

        // Calculate standard deviation
        let variance: T = values
            .iter()
            .map(|x| {
                let diff = x.clone() - mean.clone();
                diff.clone() * diff
            })
            .sum::<T>() / T::from(count);

        let std_dev = variance.sqrt();

        QualityStatistics {
            min,
            max,
            mean,
            std_dev,
            count,
        }
    }

    /// Create empty metrics for edge cases
    fn empty_metrics(&self) -> QualityMetrics<T> {
        let empty_stats = QualityStatistics {
            min: T::zero(),
            max: T::zero(),
            mean: T::zero(),
            std_dev: T::zero(),
            count: 0,
        };

        QualityMetrics {
            aspect_ratio: empty_stats.clone(),
            skewness: empty_stats.clone(),
            orthogonality: empty_stats.clone(),
            volume: empty_stats,
            num_cells: 0,
        }
    }
}
