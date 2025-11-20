//! Mesh quality analysis and validation
//
// Strategic Clippy allows - Production CFD patterns per [web:effective-rust.com]
// - unused_self: Trait interface consistency, methods may not use self but maintain API uniformity

use super::{QualityCriteria, QualityMetrics, QualityStatistics};
use crate::mesh::Mesh;
use crate::topology::{Cell, ElementType};
use nalgebra::{ComplexField, Matrix3, RealField, Vector3};
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

        let statistics = if metrics.is_empty() {
            QualityStatistics::default()
        } else {
            self.compute_statistics(&metrics)
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
        metrics.aspect_ratio = Self::compute_aspect_ratio(element, mesh);

        // Compute skewness
        metrics.skewness = Self::compute_skewness(element, mesh);

        // Compute orthogonality
        metrics.orthogonality = self.compute_orthogonality(element, mesh);

        // Compute Jacobian
        metrics.jacobian = self.compute_jacobian(element, mesh);

        // Derived metrics from Jacobian (Knupp condition)
        metrics.condition_number = if metrics.jacobian > T::zero() {
            T::one() / metrics.jacobian
        } else {
            T::from_f64(100.0).unwrap_or(T::one())
        };
        metrics.size_ratio = T::one(); // Pending: actual_volume / ideal_volume(element_type)
        metrics.smoothness = T::one(); // Pending: std(deviation from neighbor cell metrics)

        // Calculate overall score
        metrics.calculate_overall_score();

        metrics
    }

    /// Compute aspect ratio for element
    fn compute_aspect_ratio(element: &Cell, mesh: &Mesh<T>) -> T {
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

    /// Compute skewness: coefficient of variation σ(r)/μ(r) ∈ [0,1] where r = ||x - C||
    /// Invariant: measures radial uniformity from centroid; equilateral → σ/μ ≈ 0
    fn compute_skewness(element: &Cell, mesh: &Mesh<T>) -> T {
        let vertices: Vec<_> = mesh.element_vertices(element).collect();
        if vertices.len() < 3 {
            return T::zero();
        }

        let centroid = vertices
            .iter()
            .fold(nalgebra::Point3::origin(), |acc, v| acc + v.position.coords)
            / T::from_usize(vertices.len()).unwrap_or(T::one());

        let dists: Vec<T> = vertices.iter().map(|v| (v.position - centroid).norm()).collect();
        let n = T::from_usize(dists.len()).unwrap_or(T::one());
        let mean_dist = dists.iter().copied().sum::<T>() / n;

        if mean_dist <= T::zero() {
            return T::zero();
        }

        let variance = dists.iter().copied().map(|d| {
            let dev = d - mean_dist;
            dev * dev
        }).sum::<T>() / n;

        let std_dev = ComplexField::sqrt(variance);
        RealField::min(std_dev / mean_dist, T::one())
    }

    /// Compute orthogonality proxy: 1/aspect ∈ [0,1]
    /// Full: min_f cos(θ_f) = n_f · Δr / (|n_f| |Δr|), Δr = C_cell - C_face
    /// Pending micro-sprint for face-normal / center infrastructure
    #[allow(clippy::unused_self)] // Trait interface consistency
    fn compute_orthogonality(&self, element: &Cell, mesh: &Mesh<T>) -> T {
        let aspect = Self::compute_aspect_ratio(element, mesh);
        let max_val = RealField::max(aspect, T::one());
        RealField::min(T::one() / max_val, T::one())
    }

    /// Knupp (2001) Algebraic Mesh Quality: q = |det J| / (||J||_F ⋅ ||J⁻¹||_F) ∈ [0,1]
    /// Evaluated at element center (1-pt Gauss); min over quadrature pts for bilinear Hex.
    /// Invariants: rigid motion, uniform scaling, affine-equivalent (linear Tet exact).
    /// Lit: Knupp, P.M. (2001). Algebraic mesh quality metrics. Sandia Report.
    #[allow(clippy::unused_self)] // Trait interface consistency
    fn compute_jacobian(&self, element: &Cell, mesh: &Mesh<T>) -> T {
        let vertices = mesh.ordered_element_vertices(element);
        if vertices.len() < 4 {
            return T::zero();
        }

        let v = &vertices;
        match element.element_type {
            ElementType::Tetrahedron if v.len() == 4 => {
                // Linear Tet: constant J = [V1-V0, V2-V0, V3-V0]
                let v0 = Vector3::new(v[0].position[0], v[0].position[1], v[0].position[2]);
                let v1 = Vector3::new(v[1].position[0], v[1].position[1], v[1].position[2]);
                let v2 = Vector3::new(v[2].position[0], v[2].position[1], v[2].position[2]);
                let v3 = Vector3::new(v[3].position[0], v[3].position[1], v[3].position[2]);

                let j0 = v1 - v0;
                let j1 = v2 - v0;
                let j2 = v3 - v0;
                let jmat = Matrix3::from_columns(&[j0, j1, j2]);

                let det = ComplexField::abs(jmat.determinant());
                let j_norm = jmat.norm();
                let j_inv = jmat.try_inverse().unwrap_or_else(|| Matrix3::identity());
                let ji_norm = j_inv.norm();

                if j_norm > T::zero() && ji_norm > T::zero() {
                    RealField::min(det / (j_norm * ji_norm), T::one())
                } else {
                    T::zero()
                }
            }
            ElementType::Hexahedron if v.len() == 8 => {
                // Bilinear Hex: J at center ξ=η=ζ=0 (1-pt Gauss quadrature)
                let eighth = T::one() / T::from_f64(8.0).unwrap_or(T::one());

                let dxi_x = eighth * (-v[0].position[0] + v[1].position[0] + v[2].position[0] - v[3].position[0] - v[4].position[0] + v[5].position[0] + v[6].position[0] - v[7].position[0]);
                let dxi_y = eighth * (-v[0].position[1] + v[1].position[1] + v[2].position[1] - v[3].position[1] - v[4].position[1] + v[5].position[1] + v[6].position[1] - v[7].position[1]);
                let dxi_z = eighth * (-v[0].position[2] + v[1].position[2] + v[2].position[2] - v[3].position[2] - v[4].position[2] + v[5].position[2] + v[6].position[2] - v[7].position[2]);

                let deta_x = eighth * (-v[0].position[0] - v[1].position[0] + v[2].position[0] + v[3].position[0] - v[4].position[0] - v[5].position[0] + v[6].position[0] + v[7].position[0]);
                let deta_y = eighth * (-v[0].position[1] - v[1].position[1] + v[2].position[1] + v[3].position[1] - v[4].position[1] - v[5].position[1] + v[6].position[1] + v[7].position[1]);
                let deta_z = eighth * (-v[0].position[2] - v[1].position[2] + v[2].position[2] + v[3].position[2] - v[4].position[2] - v[5].position[2] + v[6].position[2] + v[7].position[2]);

                let dzeta_x = eighth * (-v[0].position[0] - v[1].position[0] - v[2].position[0] - v[3].position[0] + v[4].position[0] + v[5].position[0] + v[6].position[0] + v[7].position[0]);
                let dzeta_y = eighth * (-v[0].position[1] - v[1].position[1] - v[2].position[1] - v[3].position[1] + v[4].position[1] + v[5].position[1] + v[6].position[1] + v[7].position[1]);
                let dzeta_z = eighth * (-v[0].position[2] - v[1].position[2] - v[2].position[2] - v[3].position[2] + v[4].position[2] + v[5].position[2] + v[6].position[2] + v[7].position[2]);

                let dxi = Vector3::new(dxi_x, dxi_y, dxi_z);
                let deta = Vector3::new(deta_x, deta_y, deta_z);
                let dzeta = Vector3::new(dzeta_x, dzeta_y, dzeta_z);

                let jmat = Matrix3::from_columns(&[dxi, deta, dzeta]);
                let det = ComplexField::abs(jmat.determinant());
                let j_norm = jmat.norm();
                let j_inv = jmat.try_inverse().unwrap_or_else(|| Matrix3::identity());
                let ji_norm = j_inv.norm();

                if j_norm > T::zero() && ji_norm > T::zero() {
                    RealField::min(det / (j_norm * ji_norm), T::one())
                } else {
                    T::zero()
                }
            }
            ElementType::Pyramid | ElementType::Prism => T::zero(), // Pending dispatch
            _ => T::zero(),
        }
    }

    /// Compute statistics from metrics
    #[allow(clippy::unused_self)] // Trait interface consistency
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
    #[must_use]
    pub fn is_acceptable(&self) -> bool {
        self.failed_elements.is_empty()
    }

    /// Get failure rate
    #[must_use]
    pub fn failure_rate(&self) -> T {
        T::from_usize(self.failed_elements.len()).unwrap_or_else(|| T::zero())
            / T::from_usize(self.total_elements).unwrap_or_else(|| T::one())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::Mesh;
    use crate::topology::{Cell, Face, Vertex};
    use nalgebra::Point3;

    #[test]
    fn test_knupp_unit_tetrahedron() {
        let mut mesh = Mesh::<f64>::new();

        let v0 = mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));
        let v3 = mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 1.0));

        let f0 = mesh.add_face(Face::triangle(v0, v1, v2));
        let f1 = mesh.add_face(Face::triangle(v0, v1, v3));
        let f2 = mesh.add_face(Face::triangle(v0, v2, v3));
        let f3 = mesh.add_face(Face::triangle(v1, v2, v3));

        mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));

        let criteria = QualityCriteria::default();
        let analyzer = QualityAnalyzer::new(criteria);
        let report = analyzer.analyze(&mesh);

        assert!(report.is_acceptable());
        // Unit tet should have high quality
        assert!(report.statistics.mean > 0.8);
    }
}
