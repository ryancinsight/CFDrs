//! Coarsening quality analysis via algebraic distance measures.
//!
//! Reference: Cleary, A. J., Falgout, R. D., Henson, V. E., & Jones, J. E. (2001)

use super::CoarseningResult;
use crate::SparseMatrix;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

/// Algebraic distance-based coarsening quality measures
#[derive(Debug, Clone)]
pub struct AlgebraicDistances<T: RealField + Copy> {
    /// Algebraic distance from each point to nearest coarse point
    pub distances: Vec<T>,
    /// Average algebraic distance
    pub average_distance: T,
    /// Maximum algebraic distance
    pub max_distance: T,
    /// Points with high algebraic distance (problematic areas)
    pub high_distance_points: Vec<(usize, T)>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> AlgebraicDistances<T> {
    /// Compute algebraic distances for a coarsening result.
    ///
    /// Fine points should be algebraically "close" to coarse points for
    /// effective interpolation. Large distances indicate poor coarsening.
    pub fn compute(coarsening: &CoarseningResult<T>, matrix: &SparseMatrix<T>) -> Self {
        let n = matrix.nrows();
        let mut distances = vec![T::max_value().unwrap_or_else(T::zero); n];
        let mut total_distance = T::zero();
        let mut max_distance = T::zero();
        let mut high_distance_points = Vec::new();

        let mut is_coarse = vec![false; n];
        for &idx in &coarsening.coarse_points {
            if idx < n {
                is_coarse[idx] = true;
            }
        }

        for i in 0..n {
            if is_coarse[i] {
                distances[i] = T::zero();
            } else {
                // Check for direct strong connections first (distance = 1)
                let mut found_strong_coarse = false;
                let row = coarsening.strength_matrix.row(i);
                for &j in row.col_indices() {
                    if is_coarse[j] {
                        distances[i] = T::one();
                        found_strong_coarse = true;
                        break;
                    }
                }

                if !found_strong_coarse {
                    let distance = compute_distance_to_nearest_coarse(
                        i,
                        &is_coarse,
                        matrix,
                        &coarsening.strength_matrix,
                    );
                    distances[i] = distance;
                }

                let min_distance = distances[i];
                total_distance += min_distance;
                if min_distance > max_distance {
                    max_distance = min_distance;
                }

                if min_distance > T::from_f64(2.0).unwrap_or_else(T::zero) {
                    high_distance_points.push((i, min_distance));
                }
            }
        }

        let average_distance = if n > 0 {
            total_distance / T::from_usize(n).unwrap_or_else(T::one)
        } else {
            T::zero()
        };

        high_distance_points.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        Self {
            distances,
            average_distance,
            max_distance,
            high_distance_points,
        }
    }

    /// Assess coarsening quality based on algebraic distances
    pub fn quality_assessment(&self) -> CoarseningQuality {
        let mut quality = CoarseningQuality {
            coarsening_ratio: 0.0,
            assignment_ratio: 0.0,
            avg_interpolation_points: 0.0,
            max_interpolation_points: 0,
            coarse_points: 0,
            total_points: self.distances.len(),
            average_distance: self.average_distance.to_f64().unwrap_or(0.0),
            max_distance: self.max_distance.to_f64().unwrap_or(0.0),
            high_distance_ratio: f64::from(self.high_distance_points.len() as u32)
                / self.distances.len() as f64,
            quality_score: 0.0,
            recommendations: Vec::new(),
        };

        let distance_score: f64 = if quality.average_distance <= 1.5 {
            1.0
        } else if quality.average_distance <= 2.0 {
            0.8
        } else if quality.average_distance <= 3.0 {
            0.6
        } else {
            0.3
        };

        let max_score: f64 = if quality.max_distance <= 3.0 {
            1.0
        } else if quality.max_distance <= 5.0 {
            0.8
        } else if quality.max_distance <= 8.0 {
            0.5
        } else {
            0.2
        };

        quality.quality_score = distance_score.midpoint(max_score);

        if quality.quality_score < 0.7 {
            quality.recommendations.push(
                "Consider using Falgout or PMIS coarsening for better connectivity".to_string(),
            );
        }

        if self.high_distance_points.len() > self.distances.len() / 10 {
            quality.recommendations.push(
                "High number of poorly connected points - consider adjusting strength threshold"
                    .to_string(),
            );
        }

        if quality.max_distance > 5.0 {
            quality.recommendations.push(
                "Some points are very far from coarse levels - may cause slow convergence"
                    .to_string(),
            );
        }

        quality
    }
}

/// Compute algebraic distance to the nearest coarse point via BFS
fn compute_distance_to_nearest_coarse<T: RealField + Copy + FromPrimitive>(
    start: usize,
    is_coarse: &[bool],
    matrix: &SparseMatrix<T>,
    strength_matrix: &SparseMatrix<T>,
) -> T {
    let n = matrix.nrows();
    let mut distances = vec![T::max_value().unwrap_or_else(T::zero); n];
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();

    queue.push_back(start);
    distances[start] = T::zero();
    visited[start] = true;

    while let Some(node) = queue.pop_front() {
        let current_dist = distances[node];

        if is_coarse[node] {
            return current_dist;
        }

        // Explore neighbours via strong connections
        let row = strength_matrix.row(node);
        for (&j, &val) in row.col_indices().iter().zip(row.values().iter()) {
            if val > T::zero() && !visited[j] {
                let edge_weight = if val > T::from_f64(0.5).unwrap_or_else(T::zero) {
                    T::one()
                } else {
                    T::from_f64(2.0).unwrap_or_else(T::zero)
                };

                let new_dist = current_dist + edge_weight;
                if new_dist < distances[j] {
                    distances[j] = new_dist;
                    visited[j] = true;
                    queue.push_back(j);
                }
            }
        }

        // Also consider direct matrix connections (weaker)
        let matrix_row = matrix.row(node);
        for (&j, &val) in matrix_row
            .col_indices()
            .iter()
            .zip(matrix_row.values().iter())
        {
            if j != node && val.abs() > T::from_f64(1e-12).unwrap_or_else(T::zero) && !visited[j] {
                let edge_weight = T::from_f64(3.0).unwrap_or_else(T::zero);
                let new_dist = current_dist + edge_weight;
                if new_dist < distances[j] {
                    distances[j] = new_dist;
                    visited[j] = true;
                    queue.push_back(j);
                }
            }
        }
    }

    T::from_f64(100.0).unwrap_or_else(T::zero)
}

/// Analyze coarsening quality comprehensively
pub fn analyze_coarsening_quality<T: RealField + Copy + FromPrimitive + ToPrimitive>(
    result: &CoarseningResult<T>,
    matrix: &SparseMatrix<T>,
) -> CoarseningQuality {
    let total_points = result.fine_to_coarse_map.len();
    let coarse_points = result.coarse_points.len();
    let assigned_points = result
        .fine_to_coarse_map
        .iter()
        .filter(|x| x.is_some())
        .count();

    let coarsening_ratio = coarse_points as f64 / total_points as f64;
    let assignment_ratio = assigned_points as f64 / total_points as f64;

    let mut is_coarse = vec![false; total_points];
    for &cp in &result.coarse_points {
        if cp < total_points {
            is_coarse[cp] = true;
        }
    }

    let mut interpolation_points = Vec::new();
    for i in 0..total_points {
        if !is_coarse[i] {
            let row = result.strength_matrix.row(i);
            let connections = row.col_indices().iter().filter(|&&j| is_coarse[j]).count();
            interpolation_points.push(connections);
        }
    }

    let avg_interpolation_points = if interpolation_points.is_empty() {
        0.0
    } else {
        interpolation_points.iter().sum::<usize>() as f64 / interpolation_points.len() as f64
    };

    let max_interpolation_points = interpolation_points.iter().copied().max().unwrap_or(0);

    let algebraic = AlgebraicDistances::compute(result, matrix);
    let mut quality = algebraic.quality_assessment();
    quality.coarsening_ratio = coarsening_ratio;
    quality.assignment_ratio = assignment_ratio;
    quality.avg_interpolation_points = avg_interpolation_points;
    quality.max_interpolation_points = max_interpolation_points;
    quality.coarse_points = coarse_points;
    quality.total_points = total_points;
    quality
}

/// Quality metrics for coarsening
#[derive(Debug, Clone)]
pub struct CoarseningQuality {
    /// Ratio of coarse points to total points
    pub coarsening_ratio: f64,
    /// Ratio of assigned points to total points
    pub assignment_ratio: f64,
    /// Average number of interpolation points per F-point
    pub avg_interpolation_points: f64,
    /// Maximum number of interpolation points for any F-point
    pub max_interpolation_points: usize,
    /// Number of coarse points
    pub coarse_points: usize,
    /// Total number of points
    pub total_points: usize,
    /// Average algebraic distance from fine to coarse points
    pub average_distance: f64,
    /// Maximum algebraic distance from fine to coarse points
    pub max_distance: f64,
    /// Ratio of points with high algebraic distance (>2.0)
    pub high_distance_ratio: f64,
    /// Overall quality score (0.0-1.0, higher is better)
    pub quality_score: f64,
    /// Recommendations for improving coarsening
    pub recommendations: Vec<String>,
}

impl CoarseningQuality {
    /// Check if coarsening quality is acceptable
    pub fn is_acceptable(&self) -> bool {
        self.assignment_ratio > 0.9
            && self.coarsening_ratio > 0.1
            && self.coarsening_ratio < 0.8
            && self.max_interpolation_points <= 10
    }
}
