//! Mesh quality analysis and metrics.
//!
//! This module provides tools for analyzing mesh quality using various metrics
//! such as aspect ratio, skewness, and orthogonality.

use crate::mesh::{Mesh, Cell};
use nalgebra::{ComplexField, RealField, Point3};
use num_traits::Float;
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

impl<T: RealField + std::iter::Sum + Copy + Float> QualityAnalyzer<T>
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

    /// Analyze mesh quality using advanced iterator patterns for zero-copy efficiency
    pub fn analyze(&self, mesh: &Mesh<T>) -> QualityMetrics<T> {
        // Zero-copy quality analysis using iterator chains
        let quality_data: Vec<_> = mesh.cells
            .iter()
            .filter_map(|cell| self.analyze_cell(cell, mesh))
            .collect();

        if quality_data.is_empty() {
            return self.empty_metrics();
        }

        let num_cells = quality_data.len();

        // Single-pass metric extraction using iterator patterns
        let (aspect_ratios, skewness_values, orthogonality_values, volumes): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
            quality_data.into_iter()
                .fold(
                    (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                    |(mut ar, mut sk, mut orth, mut vol), (a, s, o, v)| {
                        ar.push(a);
                        sk.push(s);
                        orth.push(o);
                        vol.push(v);
                        (ar, sk, orth, vol)
                    }
                );

        QualityMetrics {
            aspect_ratio: self.compute_statistics(&aspect_ratios),
            skewness: self.compute_statistics(&skewness_values),
            orthogonality: self.compute_statistics(&orthogonality_values),
            volume: self.compute_statistics(&volumes),
            num_cells,
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

    /// Calculate aspect ratio - ratio of minimum to maximum edge length
    /// For high-quality elements, aspect ratio should be close to 1.0
    /// Reference: Knupp, P. "Algebraic Mesh Quality Metrics" (2001)
    fn calculate_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() < 2 {
            return T::one();
        }

        // Compute all edge lengths based on element type
        let edge_lengths = match vertices.len() {
            3 => {
                // Triangle: 3 edges
                vec![
                    (vertices[1] - vertices[0]).norm(),
                    (vertices[2] - vertices[1]).norm(),
                    (vertices[0] - vertices[2]).norm(),
                ]
            }
            4 => {
                // Tetrahedron or Quadrilateral
                if self.is_planar(vertices) {
                    // Quadrilateral: 4 edges
                    vec![
                        (vertices[1] - vertices[0]).norm(),
                        (vertices[2] - vertices[1]).norm(),
                        (vertices[3] - vertices[2]).norm(),
                        (vertices[0] - vertices[3]).norm(),
                    ]
                } else {
                    // Tetrahedron: 6 edges
                    vec![
                        (vertices[1] - vertices[0]).norm(),
                        (vertices[2] - vertices[0]).norm(),
                        (vertices[3] - vertices[0]).norm(),
                        (vertices[2] - vertices[1]).norm(),
                        (vertices[3] - vertices[1]).norm(),
                        (vertices[3] - vertices[2]).norm(),
                    ]
                }
            }
            8 => {
                // Hexahedron: 12 edges
                vec![
                    // Bottom face edges
                    (vertices[1] - vertices[0]).norm(),
                    (vertices[2] - vertices[1]).norm(),
                    (vertices[3] - vertices[2]).norm(),
                    (vertices[0] - vertices[3]).norm(),
                    // Top face edges
                    (vertices[5] - vertices[4]).norm(),
                    (vertices[6] - vertices[5]).norm(),
                    (vertices[7] - vertices[6]).norm(),
                    (vertices[4] - vertices[7]).norm(),
                    // Vertical edges
                    (vertices[4] - vertices[0]).norm(),
                    (vertices[5] - vertices[1]).norm(),
                    (vertices[6] - vertices[2]).norm(),
                    (vertices[7] - vertices[3]).norm(),
                ]
            }
            _ => {
                // General polyhedron: compute all pairwise edges
                vertices
                    .windows(2)
                    .map(|pair| (pair[1] - pair[0]).norm())
                    .collect()
            }
        };

        // Find min and max edge lengths using iterator combinators
        let min_length = edge_lengths.iter()
            .cloned()
            .fold(T::infinity(), RealField::min);
        let max_length = edge_lengths.iter()
            .cloned()
            .fold(T::neg_infinity(), RealField::max);

        if max_length > T::zero() && min_length > T::zero() {
            min_length / max_length
        } else {
            T::one()
        }
    }
    
    /// Check if vertices are coplanar within tolerance
    fn is_planar(&self, vertices: &[&Point3<T>]) -> bool {
        if vertices.len() < 4 {
            return true;
        }
        
        // Compute normal from first three vertices
        let v1 = vertices[1] - vertices[0];
        let v2 = vertices[2] - vertices[0];
        let normal = v1.cross(&v2);
        
        if normal.norm() < T::from_f64(1e-10).unwrap_or_else(T::zero) {
            return true; // Degenerate case
        }
        
        // Check if fourth vertex is in the same plane
        let v3 = vertices[3] - vertices[0];
        let dot_product = v3.dot(&normal);
        let abs_dot = if dot_product >= T::zero() { dot_product.clone() } else { -dot_product.clone() };
        let distance = abs_dot / normal.norm();
        
        distance < T::from_f64(1e-6).unwrap_or_else(T::zero)
    }

    /// Calculate skewness based on angle deviations
    fn calculate_skewness(&self, vertices: &[&Point3<T>]) -> T {
        // Skewness measures deviation from ideal element shape
        // For a tetrahedron, we check angle deviations from ideal angles
        if vertices.len() < 4 {
            return T::one(); // Maximum skewness for degenerate element
        }
        
        // Calculate all edge vectors
        let edges = vec![
            vertices[1] - vertices[0],
            vertices[2] - vertices[0],
            vertices[3] - vertices[0],
            vertices[2] - vertices[1],
            vertices[3] - vertices[1],
            vertices[3] - vertices[2],
        ];
        
        // Calculate angles between edges
        let mut max_angle_deviation = T::zero();
        let reference_angle = T::from_f64(60.0_f64.to_radians()).unwrap_or(T::one());
        
        for i in 0..edges.len() {
            for j in i+1..edges.len() {
                let dot = edges[i].dot(&edges[j]);
                let mag_i = edges[i].norm();
                let mag_j = edges[j].norm();
                
                if mag_i > T::zero() && mag_j > T::zero() {
                    let cos_angle = dot / (mag_i * mag_j);
                    // Clamp to valid range for acos
                    let cos_f64: f64 = cos_angle.to_subset().unwrap_or(0.0);
                    let cos_clamped = cos_f64.max(-1.0).min(1.0);
                    let angle = T::from_f64(cos_clamped.acos()).unwrap_or(T::zero());
                    let deviation = ComplexField::abs(angle - reference_angle.clone()) / reference_angle.clone();
                    if deviation > max_angle_deviation {
                        max_angle_deviation = deviation;
                    }
                }
            }
        }
        
        // Skewness ranges from 0 (perfect) to 1 (highly skewed)
        if max_angle_deviation > T::one() {
            T::one()
        } else {
            max_angle_deviation
        }
    }

    /// Calculate orthogonality based on face normal angles
    fn calculate_orthogonality(&self, vertices: &[&Point3<T>]) -> T {
        // Orthogonality measures how perpendicular adjacent faces are
        // For a tetrahedron, we check angles between face normals
        if vertices.len() < 4 {
            return T::zero(); // No orthogonality for degenerate element
        }
        
        // Calculate face normals
        let face_normals = vec![
            // Face 0-1-2
            (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0])),
            // Face 0-1-3
            (vertices[1] - vertices[0]).cross(&(vertices[3] - vertices[0])),
            // Face 0-2-3
            (vertices[2] - vertices[0]).cross(&(vertices[3] - vertices[0])),
            // Face 1-2-3
            (vertices[2] - vertices[1]).cross(&(vertices[3] - vertices[1])),
        ];
        
        // Calculate minimum angle between face normals
        let mut min_orthogonality = T::one();
        
        for i in 0..face_normals.len() {
            for j in i+1..face_normals.len() {
                let norm_i = face_normals[i].norm();
                let norm_j = face_normals[j].norm();
                
                if norm_i > T::zero() && norm_j > T::zero() {
                    let dot = face_normals[i].dot(&face_normals[j]);
                    let cos_angle = dot / (norm_i * norm_j);
                    // Orthogonality: 1 when perpendicular (cos = 0), 0 when parallel (|cos| = 1)
                    let orthogonality = T::one() - ComplexField::abs(cos_angle);
                    if orthogonality < min_orthogonality {
                        min_orthogonality = orthogonality;
                    }
                }
            }
        }
        
        // Return value between 0 (parallel faces) and 1 (perpendicular faces)
        if min_orthogonality < T::zero() {
            T::zero()
        } else if min_orthogonality > T::one() {
            T::one()
        } else {
            min_orthogonality
        }
    }

    /// Calculate cell volume for tetrahedron
    fn calculate_volume(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() < 4 {
            return T::zero();
        }

        // Volume of tetrahedron: V = |det(v1, v2, v3)| / 6
        // where v1, v2, v3 are edge vectors from vertex 0
        let v1 = vertices[1] - vertices[0];
        let v2 = vertices[2] - vertices[0];
        let v3 = vertices[3] - vertices[0];

        let six = T::from_f64(6.0).unwrap_or(T::from_usize(6).unwrap());
        ComplexField::abs(v1.cross(&v2).dot(&v3)) / six
    }

    /// Compute statistics for a set of values using advanced iterator patterns
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

        // Use advanced iterator combinators for zero-copy statistics computation
        use cfd_math::MathIteratorExt;

        let min = values.iter().cloned().fold(T::infinity(), RealField::min);
        let max = values.iter().cloned().fold(T::neg_infinity(), RealField::max);
        let count = values.len();

        // Use zero-copy slice operations for better performance
        let mean = values.iter().cloned().mean().unwrap_or_else(T::zero);
        let variance = values.iter().cloned().variance().unwrap_or_else(T::zero);
        let std_dev = ComplexField::sqrt(variance);

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
