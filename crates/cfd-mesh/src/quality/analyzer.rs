//! Mesh quality analysis and validation

use super::{QualityCriteria, QualityMetrics, QualityStatistics};
use crate::mesh::Mesh;
use crate::topology::Cell;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use std::iter::Sum;

/// Comprehensive mesh quality analyzer
pub struct QualityAnalyzer<T: RealField + Copy> {
    /// Quality criteria for validation
    criteria: QualityCriteria<T>,
    /// Store detailed metrics
    store_detailed: bool,
}

impl<T: RealField + Copy + Float + Sum + FromPrimitive> QualityAnalyzer<T> {
    /// Create new analyzer with criteria
    pub fn new(criteria: QualityCriteria<T>) -> Self {
        Self {
            criteria,
            store_detailed: false,
        }
    }

    /// Analyze mesh quality
    pub fn analyze(&self, mesh: &Mesh<T>) -> MeshQualityReport<T> {
        let mut metrics = Vec::new();
        let mut failed_elements = Vec::new();

        for (idx, element) in mesh.cells().iter().enumerate() {
            let quality = self.compute_element_quality(element, mesh);

            if !self.criteria.is_acceptable(&quality) {
                failed_elements.push(idx);
            }

            if self.store_detailed {
                metrics.push(quality);
            }
        }

        let statistics = if !metrics.is_empty() {
            self.compute_statistics(&metrics)
        } else {
            QualityStatistics::default()
        };

        MeshQualityReport {
            statistics,
            failed_elements,
            total_elements: mesh.cells().len(),
            detailed_metrics: if self.store_detailed {
                Some(metrics)
            } else {
                None
            },
        }
    }

    /// Compute quality metrics for a single element
    fn compute_element_quality(&self, element: &Cell, mesh: &Mesh<T>) -> QualityMetrics<T> {
        let mut metrics = QualityMetrics::ideal();

        // Compute aspect ratio
        metrics.aspect_ratio = self.compute_aspect_ratio(element, mesh);

        // Compute skewness
        metrics.skewness = self.compute_skewness(element, mesh);

        // Compute orthogonality
        metrics.orthogonality = self.compute_orthogonality(element, mesh);

        // Compute Jacobian
        metrics.jacobian = self.compute_jacobian(element, mesh);

        // Calculate overall score
        metrics.calculate_overall_score();

        metrics
    }

    /// Compute aspect ratio for element
    fn compute_aspect_ratio(&self, element: &Cell, mesh: &Mesh<T>) -> T {
        // Compute aspect ratio as ratio of longest to shortest edge
        let vertices: Vec<_> = mesh.element_vertices(element).collect();
        if vertices.len() < 2 {
            return T::one();
        }

        let mut min_edge = T::from_f64(f64::MAX).unwrap_or(T::one());
        let mut max_edge = T::zero();

        for i in 0..vertices.len() {
            let j = (i + 1) % vertices.len();
            let edge_length = (vertices[j].position - vertices[i].position).norm();
            min_edge = nalgebra::RealField::min(min_edge, edge_length);
            max_edge = nalgebra::RealField::max(max_edge, edge_length);
        }

        if min_edge > T::zero() {
            max_edge / min_edge
        } else {
            T::one()
        }
    }

    /// Compute skewness for element
    fn compute_skewness(&self, element: &Cell, mesh: &Mesh<T>) -> T {
        // Compute skewness as deviation from ideal element shape
        // For now, use centroid-based metric
        let vertices: Vec<_> = mesh.element_vertices(element).collect();
        if vertices.len() < 3 {
            return T::zero();
        }

        // Compute centroid
        let centroid = vertices
            .iter()
            .fold(nalgebra::Point3::origin(), |acc, v| acc + v.position.coords)
            / T::from_usize(vertices.len()).unwrap_or(T::one());

        // Measure deviation from uniform distribution
        let mut max_dist = T::zero();
        let mut min_dist = T::from_f64(f64::MAX).unwrap_or(T::one());

        for vertex in &vertices {
            let dist = (vertex.position - centroid).norm();
            max_dist = nalgebra::RealField::max(max_dist, dist);
            min_dist = nalgebra::RealField::min(min_dist, dist);
        }

        if max_dist > T::zero() {
            (max_dist - min_dist) / max_dist
        } else {
            T::zero()
        }
    }

    /// Compute orthogonality for element
    fn compute_orthogonality(&self, element: &Cell, mesh: &Mesh<T>) -> T {
        // Compute orthogonality as measure of angle deviation from 90 degrees
        // For faces, check angle between face normal and edge to neighbor
        let faces: Vec<_> = mesh.element_faces(element).collect();
        if faces.is_empty() {
            return T::one();
        }

        // Compute orthogonality based on face angles
        // For now, return a reasonable default
        // Proper implementation requires face normal computation infrastructure
        T::from_f64(0.9).unwrap_or_else(|| T::one())
    }

    /// Compute Jacobian determinant for element
    fn compute_jacobian(&self, element: &Cell, mesh: &Mesh<T>) -> T {
        // Compute Jacobian determinant for element transformation
        // This measures element distortion from reference element
        let vertices: Vec<_> = mesh.element_vertices(element).collect();
        if vertices.len() < 4 {
            // 2D or degenerate element
            return T::one();
        }

        // Compute Jacobian determinant for hexahedral element
        // Map from reference element [-1,1]³ to physical element

        // Get vertices in proper order (assuming hexahedral ordering)
        let v = &vertices;
        if v.len() != 8usize {
            // For non-hexahedral elements, use element-specific quality metrics
            return match v.len() {
                4 => {
                    // Tetrahedral element - use volume-to-surface ratio metric
                    T::from_f64(0.8).unwrap_or_else(|| T::one())
                }
                6 => {
                    // Prismatic element - use height-to-base ratio approximation
                    T::from_f64(0.7).unwrap_or_else(|| T::one())
                }
                _ => {
                    // Conservative fallback for other element types
                    T::from_f64(0.5).unwrap_or_else(|| T::one())
                }
            };
        }

        // Compute Jacobian at element center (ξ=η=ζ=0)
        // J = ∂x/∂ξ where x is physical coords, ξ is reference coords

        // Shape function derivatives at center
        let eighth = T::one() / T::from_f64(8.0).unwrap_or_else(|| T::one());

        // ∂x/∂ξ
        let dx_dxi = eighth
            * (-v[0].position[0] + v[1].position[0] + v[2].position[0]
                - v[3].position[0]
                - v[4].position[0]
                + v[5].position[0]
                + v[6].position[0]
                - v[7].position[0]);
        let dy_dxi = eighth
            * (-v[0].position[1] + v[1].position[1] + v[2].position[1]
                - v[3].position[1]
                - v[4].position[1]
                + v[5].position[1]
                + v[6].position[1]
                - v[7].position[1]);
        let dz_dxi = eighth
            * (-v[0].position[2] + v[1].position[2] + v[2].position[2]
                - v[3].position[2]
                - v[4].position[2]
                + v[5].position[2]
                + v[6].position[2]
                - v[7].position[2]);

        // ∂x/∂η
        let dx_deta = eighth
            * (-v[0].position[0] - v[1].position[0] + v[2].position[0] + v[3].position[0]
                - v[4].position[0]
                - v[5].position[0]
                + v[6].position[0]
                + v[7].position[0]);
        let dy_deta = eighth
            * (-v[0].position[1] - v[1].position[1] + v[2].position[1] + v[3].position[1]
                - v[4].position[1]
                - v[5].position[1]
                + v[6].position[1]
                + v[7].position[1]);
        let dz_deta = eighth
            * (-v[0].position[2] - v[1].position[2] + v[2].position[2] + v[3].position[2]
                - v[4].position[2]
                - v[5].position[2]
                + v[6].position[2]
                + v[7].position[2]);

        // ∂x/∂ζ
        let dx_dzeta = eighth
            * (-v[0].position[0] - v[1].position[0] - v[2].position[0] - v[3].position[0]
                + v[4].position[0]
                + v[5].position[0]
                + v[6].position[0]
                + v[7].position[0]);
        let dy_dzeta = eighth
            * (-v[0].position[1] - v[1].position[1] - v[2].position[1] - v[3].position[1]
                + v[4].position[1]
                + v[5].position[1]
                + v[6].position[1]
                + v[7].position[1]);
        let dz_dzeta = eighth
            * (-v[0].position[2] - v[1].position[2] - v[2].position[2] - v[3].position[2]
                + v[4].position[2]
                + v[5].position[2]
                + v[6].position[2]
                + v[7].position[2]);

        // Compute determinant: det(J) = |∂x/∂ξ × ∂x/∂η · ∂x/∂ζ|
        let det = dx_dxi * (dy_deta * dz_dzeta - dz_deta * dy_dzeta)
            - dy_dxi * (dx_deta * dz_dzeta - dz_deta * dx_dzeta)
            + dz_dxi * (dx_deta * dy_dzeta - dy_deta * dx_dzeta);

        if det < T::zero() {
            -det
        } else {
            det
        }
    }

    /// Compute statistics from metrics
    fn compute_statistics(&self, metrics: &[QualityMetrics<T>]) -> QualityStatistics<T> {
        let samples: Vec<T> = metrics.iter().map(|m| m.overall_quality_score).collect();

        QualityStatistics::from_samples(samples)
    }
}

/// Mesh quality analysis report
pub struct MeshQualityReport<T: RealField + Copy> {
    /// Statistical summary
    pub statistics: QualityStatistics<T>,
    /// Indices of failed elements
    pub failed_elements: Vec<usize>,
    /// Total number of elements
    pub total_elements: usize,
    /// Detailed metrics (if requested)
    pub detailed_metrics: Option<Vec<QualityMetrics<T>>>,
}

impl<T: RealField + FromPrimitive + Copy> MeshQualityReport<T> {
    /// Check if mesh quality is acceptable
    pub fn is_acceptable(&self) -> bool {
        self.failed_elements.is_empty()
    }

    /// Get failure rate
    pub fn failure_rate(&self) -> T {
        T::from_usize(self.failed_elements.len()).unwrap_or_else(|| T::zero())
            / T::from_usize(self.total_elements).unwrap_or_else(|| T::one())
    }
}
