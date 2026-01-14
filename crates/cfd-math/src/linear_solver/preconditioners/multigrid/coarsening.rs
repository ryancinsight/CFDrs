//! Coarsening strategies for AMG hierarchy construction

use crate::SparseMatrix;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

/// Result of coarsening operation
#[derive(Debug, Clone)]
pub struct CoarseningResult<T: RealField + Copy> {
    /// Indices of coarse points (C-points)
    pub coarse_points: Vec<usize>,
    /// Mapping from fine points to coarse points (None for F-points)
    pub fine_to_coarse_map: Vec<Option<usize>>,
    /// Strength of connection matrix
    pub strength_matrix: SparseMatrix<T>,
}

/// Ruge-Stüben coarsening algorithm
pub fn ruge_stueben_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();
    let mut coarse_points = Vec::new();
    let mut fine_to_coarse_map = vec![None; n];

    // Step 1: Compute strength of connection matrix
    let strength_matrix = compute_strength_matrix(matrix, strength_threshold)?;

    // Transpose strength matrix to get S^T (influence graph)
    let strength_transpose = strength_matrix.transpose();

    // Step 2: Initialize lambda (measure of importance)
    let mut lambda = vec![0; n];
    let st_offsets = strength_transpose.row_offsets();
    for i in 0..n {
        lambda[i] = st_offsets[i + 1] - st_offsets[i];
    }

    // Status: 0 = Undecided, 1 = C-point, 2 = F-point
    let mut status = vec![0; n];
    let mut undecided_count = n;

    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();
    let st_indices = strength_transpose.col_indices();

    // Step 3: C/F Splitting (Standard Ruge-Stueben First Pass)
    while undecided_count > 0 {
        // Pick point with max lambda among undecided
        let mut max_lambda = -1;
        let mut i_opt = None;

        for k in 0..n {
            if status[k] == 0 && lambda[k] as i32 > max_lambda {
                max_lambda = lambda[k] as i32;
                i_opt = Some(k);
            }
        }

        if let Some(i) = i_opt {
            // Make i a C-point
            status[i] = 1;
            undecided_count -= 1;
            coarse_points.push(i);
            fine_to_coarse_map[i] = Some(coarse_points.len() - 1);

            // For all undecided j that are strongly influenced by i (j in S_i^T)
            for k in st_offsets[i]..st_offsets[i + 1] {
                let j = st_indices[k];
                if status[j] == 0 {
                    // Make j an F-point
                    status[j] = 2;
                    undecided_count -= 1;

                    // For all undecided k that strongly influence j (k in S_j)
                    for m in s_offsets[j]..s_offsets[j + 1] {
                        let kk = s_indices[m];
                        if status[kk] == 0 {
                            lambda[kk] += 1;
                        }
                    }
                }
            }
        } else {
            break;
        }
    }

    // Step 4: Map F-points to their strongest connected C-point
    for i in 0..n {
        if fine_to_coarse_map[i].is_none() {
            let mut max_strength = T::from_f64(-1.0).unwrap_or_else(T::zero);
            let mut best_coarse_idx = None;

            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                // Check if j is a C-point
                if status[j] == 1 {
                    let strength = strength_matrix.values()[k];
                    if strength > max_strength {
                        max_strength = strength;
                        // Find the index of j in coarse_points
                        best_coarse_idx = fine_to_coarse_map[j];
                    }
                }
            }

            if let Some(idx) = best_coarse_idx {
                fine_to_coarse_map[i] = Some(idx);
            }
        }
    }

    Ok(CoarseningResult {
        coarse_points,
        fine_to_coarse_map,
        strength_matrix,
    })
}

/// Aggregation-based coarsening
pub fn aggregation_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    max_aggregate_size: usize,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();
    let mut coarse_points = Vec::new();
    let mut fine_to_coarse_map = vec![None; n];
    let mut aggregated = vec![false; n];

    // Compute strength matrix for aggregation decisions
    let strength_matrix =
        compute_strength_matrix(matrix, T::from_f64(0.5).unwrap_or_else(T::zero))?;
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    let mut aggregate_id = 0;

    for i in 0..n {
        if !aggregated[i] {
            // Start new aggregate with point i as seed
            let mut aggregate = vec![i];
            aggregated[i] = true;

            // Find strongly connected neighbors to add to aggregate
            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                if !aggregated[j] && aggregate.len() < max_aggregate_size {
                    aggregate.push(j);
                    aggregated[j] = true;
                }
            }

            // Choose representative (first point) as coarse point
            coarse_points.push(aggregate[0]);

            // Map all points in aggregate to this coarse point
            for &point in &aggregate {
                fine_to_coarse_map[point] = Some(aggregate_id);
            }

            aggregate_id += 1;
        }
    }

    Ok(CoarseningResult {
        coarse_points,
        fine_to_coarse_map,
        strength_matrix,
    })
}

/// Hybrid coarsening (Ruge-Stüben with aggregation fallback)
pub fn hybrid_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
    max_aggregate_size: usize,
) -> Result<CoarseningResult<T>> {
    // Try Ruge-Stüben first
    match ruge_stueben_coarsening(matrix, strength_threshold) {
        Ok(result) => {
            // Check if result is reasonable (not too many isolated points)
            let assigned_points = result
                .fine_to_coarse_map
                .iter()
                .filter(|x| x.is_some())
                .count();
            let assignment_ratio = assigned_points as f64 / matrix.nrows() as f64;

            if assignment_ratio > 0.8 {
                Ok(result)
            } else {
                // Fall back to aggregation
                aggregation_coarsening(matrix, max_aggregate_size)
            }
        }
        Err(_) => {
            // Fall back to aggregation
            aggregation_coarsening(matrix, max_aggregate_size)
        }
    }
}

/// Falgout coarsening algorithm (CLJP method)
/// Reference: Falgout, R. D. (2006). An introduction to algebraic multigrid
pub fn falgout_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();
    let mut coarse_points = Vec::new();
    let mut fine_to_coarse_map = vec![None; n];

    // Step 1: Compute strength of connection matrix
    let strength_matrix = compute_strength_matrix(matrix, strength_threshold)?;
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    // Step 2: Compute measures for each point
    let mut measures = vec![0.0; n];
    for i in 0..n {
        // Measure based on strong connections count
        measures[i] = (s_offsets[i + 1] - s_offsets[i]) as f64;
    }

    // Step 3: Sort points by decreasing measure
    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_by(|&a, &b| measures[b].partial_cmp(&measures[a]).unwrap());

    // Step 4: CLJP (Compatible Relaxation with Lambda = 4/3)
    let lambda = 4.0 / 3.0;
    let mut unassigned: Vec<usize> = sorted_indices.clone();
    let mut status = vec![0; n]; // 0: unassigned, 1: coarse, 2: fine

    while !unassigned.is_empty() {
        // Take first point from sorted list
        let i = unassigned[0];
        unassigned.remove(0);

        if status[i] != 0 {
            continue;
        }

        // Check if point should be made coarse
        let mut should_be_coarse = true;

        // Count strongly connected coarse neighbors
        let mut coarse_neighbors = 0;
        let total_strong_connections = s_offsets[i + 1] - s_offsets[i];

        for k in s_offsets[i]..s_offsets[i + 1] {
            let j = s_indices[k];
            if status[j] == 1 {
                coarse_neighbors += 1;
            }
        }

        // Apply lambda criterion
        if total_strong_connections > 0 {
            let ratio = f64::from(coarse_neighbors) / total_strong_connections as f64;
            if ratio >= lambda {
                should_be_coarse = false;
            }
        }

        if should_be_coarse {
            // Make point coarse
            status[i] = 1;
            coarse_points.push(i);
            fine_to_coarse_map[i] = Some(coarse_points.len() - 1);

            // Mark strongly connected unassigned points as fine
            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                if status[j] == 0 {
                    status[j] = 2;
                    fine_to_coarse_map[j] = Some(coarse_points.len() - 1);
                }
            }
        } else {
            status[i] = 2;
            // Map to strongest coarse neighbor if exists, else it stays None for now
            // (Standard CLJP would handle this more rigorously)
        }

        // Re-filter unassigned
        unassigned.retain(|&idx| status[idx] == 0);
    }

    Ok(CoarseningResult {
        coarse_points,
        fine_to_coarse_map,
        strength_matrix,
    })
}

/// PMIS (Parallel Modified Independent Set) coarsening
/// Reference: Luby's algorithm adapted for parallel coarsening
pub fn pmis_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();
    let mut coarse_points = Vec::new();
    let mut fine_to_coarse_map = vec![None; n];

    // Step 1: Compute strength of connection matrix
    let strength_matrix = compute_strength_matrix(matrix, strength_threshold)?;
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    // Step 2: Initialize all points as undecided
    let mut status = vec![0; n]; // 0: undecided, 1: coarse, 2: fine
    let mut undecided: Vec<usize> = (0..n).collect();

    // Step 3: PMIS algorithm
    while !undecided.is_empty() {
        // Random permutation for tie-breaking (simplified deterministic version)
        let mut candidates: Vec<usize> = undecided.clone();
        candidates.sort_by_key(|&i| {
            // Use degree as tie-breaker (higher degree first)
            let degree = s_offsets[i + 1] - s_offsets[i];
            -(degree as i32) // Negative for descending order
        });

        for &i in &candidates {
            if status[i] != 0 {
                continue; // Already decided
            }

            // Check if all strongly connected undecided neighbors have lower priority
            let mut can_be_coarse = true;
            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                if status[j] == 0 {
                    // Compare priorities (using indices as deterministic priority)
                    if j < i {
                        can_be_coarse = false;
                        break;
                    }
                }
            }

            if can_be_coarse {
                // Make point coarse
                status[i] = 1; // coarse
                coarse_points.push(i);
                fine_to_coarse_map[i] = Some(coarse_points.len() - 1);

                // Mark strongly connected neighbors as fine
                for k in s_offsets[i]..s_offsets[i + 1] {
                    let j = s_indices[k];
                    if status[j] == 0 {
                        status[j] = 2; // fine
                        fine_to_coarse_map[j] = Some(coarse_points.len() - 1);
                    }
                }
            }
        }

        // Remove decided points from undecided list
        undecided.retain(|&i| status[i] == 0);
    }

    Ok(CoarseningResult {
        coarse_points,
        fine_to_coarse_map,
        strength_matrix,
    })
}

/// HMIS (Hybrid Modified Independent Set) coarsening
/// Combines PMIS with aggressive coarsening for better parallel performance
pub fn hmis_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
    aggressive_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();

    // First try PMIS
    match pmis_coarsening(matrix, strength_threshold) {
        Ok(result) => {
            // Check coarsening ratio
            let coarsening_ratio =
                T::from_f64(result.coarse_points.len() as f64 / n as f64).unwrap_or_else(T::zero);

            // If coarsening ratio is too low, apply aggressive coarsening
            if coarsening_ratio < aggressive_threshold {
                aggressive_coarsening(matrix, strength_threshold, aggressive_threshold)
            } else {
                Ok(result)
            }
        }
        Err(_) => {
            // Fall back to aggressive coarsening
            aggressive_coarsening(matrix, strength_threshold, aggressive_threshold)
        }
    }
}

/// Aggressive coarsening strategy for difficult matrices
fn aggressive_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    _strength_threshold: T,
    _target_ratio: T,
) -> Result<CoarseningResult<T>> {
    // Use modified aggregation with larger aggregates
    let max_aggregate_size = 8; // Larger than standard
    aggregation_coarsening(matrix, max_aggregate_size)
}

/// Algebraic distance-based coarsening quality measures
/// Reference: Cleary, A. J., Falgout, R. D., Henson, V. E., & Jones, J. E. (2001)
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
    /// Compute algebraic distances for a coarsening result
    /// Based on the concept that fine points should be "close" to coarse points
    pub fn compute(coarsening: &CoarseningResult<T>, matrix: &SparseMatrix<T>) -> Self {
        let n = matrix.nrows();
        let mut distances = vec![T::max_value().unwrap_or_else(T::zero); n];
        let mut total_distance = T::zero();
        let mut max_distance = T::zero();
        let mut high_distance_points = Vec::new();

        // For each fine point, find minimum algebraic distance to coarse points
        for i in 0..n {
            if coarsening.coarse_points.contains(&i) {
                // Coarse points have distance 0
                distances[i] = T::zero();
            } else {
                // Find minimum distance to coarse points via strongly connected paths
                let mut min_distance = T::max_value().unwrap_or_else(T::zero);

                for &coarse_idx in &coarsening.coarse_points {
                    let is_strong =
                        if let Some(val) = coarsening.strength_matrix.get_entry(i, coarse_idx) {
                            val.into_value() > T::zero()
                        } else {
                            false
                        };

                    if is_strong {
                        // Direct strong connection - distance 1
                        min_distance = if min_distance < T::one() {
                            min_distance
                        } else {
                            T::one()
                        };
                    } else {
                        // Compute algebraic distance through intermediate points
                        let distance = compute_algebraic_distance(
                            i,
                            coarse_idx,
                            matrix,
                            &coarsening.strength_matrix,
                        );
                        min_distance = if min_distance < distance {
                            min_distance
                        } else {
                            distance
                        };
                    }
                }

                distances[i] = min_distance;
                total_distance += min_distance;
                if min_distance > max_distance {
                    max_distance = min_distance;
                }

                // Track points with high algebraic distance
                if min_distance > T::from_f64(2.0).unwrap_or_else(T::zero) {
                    // Threshold for "high distance"
                    high_distance_points.push((i, min_distance));
                }
            }
        }

        let average_distance = if n > 0 {
            total_distance / T::from_usize(n).unwrap_or_else(T::one)
        } else {
            T::zero()
        };

        // Sort by distance (highest first)
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
            coarsening_ratio: 0.0,              // Not computed in this method
            assignment_ratio: 0.0,              // Not computed in this method
            avg_interpolation_points: 0.0,      // Not computed in this method
            max_interpolation_points: 0,        // Not computed in this method
            coarse_points: 0,                   // Not computed in this method
            total_points: self.distances.len(), // Use distances length as total points
            average_distance: self.average_distance.to_f64().unwrap_or(0.0),
            max_distance: self.max_distance.to_f64().unwrap_or(0.0),
            high_distance_ratio: f64::from(self.high_distance_points.len() as u32)
                / self.distances.len() as f64,
            quality_score: 0.0,
            recommendations: Vec::new(),
        };

        // Quality score based on distance metrics
        // Lower distances = better quality
        let distance_score: f64 = if quality.average_distance <= 1.5 {
            1.0 // Excellent
        } else if quality.average_distance <= 2.0 {
            0.8 // Good
        } else if quality.average_distance <= 3.0 {
            0.6 // Fair
        } else {
            0.3 // Poor
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

        // Generate recommendations
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

/// Compute algebraic distance between two points
/// Uses a simplified breadth-first search on the strength graph
fn compute_algebraic_distance<T: RealField + Copy + FromPrimitive>(
    start: usize,
    target: usize,
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

        if node == target {
            return current_dist;
        }

        // Explore neighbors via strong connections
        let row = strength_matrix.row(node);
        for (&j, &val) in row.col_indices().iter().zip(row.values().iter()) {
            if val > T::zero() && !visited[j] {
                let edge_weight = if val > T::from_f64(0.5).unwrap_or_else(T::zero) {
                    T::one() // Strong connection
                } else {
                    T::from_f64(2.0).unwrap_or_else(T::zero) // Weaker connection
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
                let edge_weight = T::from_f64(3.0).unwrap_or_else(T::zero); // Direct matrix connection (weaker)
                let new_dist = current_dist + edge_weight;
                if new_dist < distances[j] {
                    distances[j] = new_dist;
                    visited[j] = true;
                    queue.push_back(j);
                }
            }
        }
    }

    // If no path found, return large distance
    T::from_f64(100.0).unwrap_or_else(T::zero)
}

/// Compute strength of connection matrix
fn compute_strength_matrix<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<SparseMatrix<T>> {
    let n = matrix.nrows();
    let mut row_offsets = vec![0; n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    let m_offsets = matrix.row_offsets();
    let m_indices = matrix.col_indices();
    let m_values = matrix.values();

    for i in 0..n {
        // Find maximum off-diagonal element in row i
        let mut max_off_diag: T = T::zero();
        for k in m_offsets[i]..m_offsets[i + 1] {
            let j = m_indices[k];
            if i != j {
                max_off_diag = if max_off_diag > m_values[k].abs() {
                    max_off_diag
                } else {
                    m_values[k].abs()
                };
            }
        }

        // Add strong connections
        for k in m_offsets[i]..m_offsets[i + 1] {
            let j = m_indices[k];
            if i != j {
                let strength = m_values[k].abs();
                if strength >= strength_threshold * max_off_diag {
                    col_indices.push(j);
                    values.push(strength);
                }
            }
        }
        row_offsets[i + 1] = col_indices.len();
    }

    SparseMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values).map_err(|e| {
        Error::Numerical(NumericalErrorKind::InvalidValue {
            value: format!("Failed to create strength matrix: {e}"),
        })
    })
}

/// Analyze coarsening quality
pub fn analyze_coarsening_quality<T: RealField + Copy + FromPrimitive>(
    result: &CoarseningResult<T>,
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

    // Create fast lookup for coarse points
    let mut is_coarse = vec![false; total_points];
    for &cp in &result.coarse_points {
        if cp < total_points {
            is_coarse[cp] = true;
        }
    }

    // Analyze interpolation complexity
    let mut interpolation_points = Vec::new();
    for i in 0..total_points {
        // If point is NOT coarse (i.e., it's an F-point), calculate interpolation
        if !is_coarse[i] {
            // Count how many coarse points this fine point connects to
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

    CoarseningQuality {
        coarsening_ratio,
        assignment_ratio,
        avg_interpolation_points,
        max_interpolation_points,
        coarse_points,
        total_points,
        average_distance: 0.0,       // Not computed in this basic analysis
        max_distance: 0.0,           // Not computed in this basic analysis
        high_distance_ratio: 0.0,    // Not computed in this basic analysis
        quality_score: 0.0,          // Not computed in this basic analysis
        recommendations: Vec::new(), // Not computed in this basic analysis
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    use nalgebra_sparse::CsrMatrix;

    fn create_test_matrix() -> CsrMatrix<f64> {
        // Create a simple 2D Poisson matrix
        let n = 4;
        let mut matrix = DMatrix::zeros(n * n, n * n);

        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                matrix[(idx, idx)] = 4.0; // Main diagonal

                // Neighbors
                if i > 0 {
                    matrix[(idx, (i - 1) * n + j)] = -1.0;
                }
                if i < n - 1 {
                    matrix[(idx, (i + 1) * n + j)] = -1.0;
                }
                if j > 0 {
                    matrix[(idx, i * n + (j - 1))] = -1.0;
                }
                if j < n - 1 {
                    matrix[(idx, i * n + (j + 1))] = -1.0;
                }
            }
        }

        CsrMatrix::from(&matrix)
    }

    #[test]
    fn test_ruge_stueben_coarsening() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();

        assert!(!result.coarse_points.is_empty());
        assert_eq!(result.fine_to_coarse_map.len(), matrix.nrows());
        assert_eq!(result.strength_matrix.nrows(), matrix.nrows());

        let quality = analyze_coarsening_quality(&result);
        assert!(quality.coarsening_ratio > 0.0);
        assert!(quality.assignment_ratio > 0.0);
    }

    #[test]
    fn test_aggregation_coarsening() {
        let matrix = create_test_matrix();
        let result = aggregation_coarsening(&matrix, 4).unwrap();

        assert!(!result.coarse_points.is_empty());
        assert_eq!(result.fine_to_coarse_map.len(), matrix.nrows());

        let quality = analyze_coarsening_quality(&result);
        assert!(quality.coarsening_ratio > 0.0);
    }

    #[test]
    fn test_coarsening_quality_analysis() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();
        let quality = analyze_coarsening_quality(&result);

        assert!(quality.coarsening_ratio > 0.0 && quality.coarsening_ratio <= 1.0);
        assert!(quality.assignment_ratio >= 0.0 && quality.assignment_ratio <= 1.0);
        assert!(quality.avg_interpolation_points >= 0.0);
    }

    #[test]
    fn test_hybrid_coarsening() {
        let matrix = create_test_matrix();
        let result = hybrid_coarsening(&matrix, 0.25, 4).unwrap();

        assert!(!result.coarse_points.is_empty());
        assert_eq!(result.fine_to_coarse_map.len(), matrix.nrows());
    }

    #[test]
    fn test_mapping_correctness() {
        let matrix = create_test_matrix();

        // Test all major coarsening algorithms
        let results = vec![
            (
                "Ruge-Stueben",
                ruge_stueben_coarsening(&matrix, 0.25).unwrap(),
            ),
            ("Aggregation", aggregation_coarsening(&matrix, 4).unwrap()),
            ("Falgout", falgout_coarsening(&matrix, 0.25).unwrap()),
            ("PMIS", pmis_coarsening(&matrix, 0.25).unwrap()),
        ];

        for (name, result) in results {
            // Verify that all mapped points point to valid coarse indices
            for (i, &map) in result.fine_to_coarse_map.iter().enumerate() {
                if let Some(coarse_idx) = map {
                    assert!(
                        coarse_idx < result.coarse_points.len(),
                        "[{name}] Mapped index out of bounds"
                    );

                    // If i is a coarse point, it should map to its own index in coarse_points
                    if result.coarse_points.contains(&i) {
                        assert_eq!(
                            result.coarse_points[coarse_idx], i,
                            "[{name}] Coarse point mapped to wrong index"
                        );
                    }
                }
            }

            // Verify that every coarse point is correctly mapped to its index in coarse_points
            for (idx, &cp) in result.coarse_points.iter().enumerate() {
                assert_eq!(
                    result.fine_to_coarse_map[cp],
                    Some(idx),
                    "[{name}] Coarse point {cp} not correctly mapped"
                );
            }
        }
    }

    #[test]
    fn test_coarsening_ratio_bounds() {
        // Create larger test matrix for more realistic coarsening
        // 10x10 grid = 100 points (versus 4x4 = 16 which is too small)
        let n = 10;
        let mut matrix = DMatrix::zeros(n * n, n * n);

        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                matrix[(idx, idx)] = 4.0; // Main diagonal

                // Neighbors (5-point stencil)
                if i > 0 {
                    matrix[(idx, (i - 1) * n + j)] = -1.0;
                }
                if i < n - 1 {
                    matrix[(idx, (i + 1) * n + j)] = -1.0;
                }
                if j > 0 {
                    matrix[(idx, i * n + (j - 1))] = -1.0;
                }
                if j < n - 1 {
                    matrix[(idx, i * n + (j + 1))] = -1.0;
                }
            }
        }

        let sparse_matrix = CsrMatrix::from(&matrix);
        let result = ruge_stueben_coarsening(&sparse_matrix, 0.25).unwrap();

        let n_total = sparse_matrix.nrows();
        let n_coarse = result.coarse_points.len();
        let ratio = n_coarse as f64 / n_total as f64;

        println!("Coarsening ratio: {ratio} ({n_coarse} coarse / {n_total} total)");

        // For a 2D Laplacian, Ruge-Stüben typically produces ratio ~0.25 - 0.5
        // Standard coarsening picks roughly every other point in each direction
        assert!(ratio >= 0.15, "Coarsening ratio too low: {ratio}");
        assert!(ratio <= 0.65, "Coarsening ratio too high: {ratio}");
    }

    #[test]
    fn test_interpolation_operator_shape() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();

        let n_fine = matrix.nrows();
        let n_coarse = result.coarse_points.len();

        // Verify fine_to_coarse_map length
        assert_eq!(result.fine_to_coarse_map.len(), n_fine);

        // Verify that all mapped indices are within [0, n_coarse)
        for &map in &result.fine_to_coarse_map {
            if let Some(idx) = map {
                assert!(
                    idx < n_coarse,
                    "Coarse index {idx} out of bounds (n_coarse={n_coarse})"
                );
            }
        }

        // Verify that every coarse index is used at least once (by the coarse point itself)
        let mut used_indices = vec![false; n_coarse];
        for &map in &result.fine_to_coarse_map {
            if let Some(idx) = map {
                used_indices[idx] = true;
            }
        }
        assert!(
            used_indices.iter().all(|&x| x),
            "Not all coarse indices are utilized"
        );
    }

    #[test]
    fn test_coarsening_convergence_behavior() {
        // This test verifies that the coarsening produces a hierarchy that
        // *should* lead to convergence. It checks:
        // 1. Coarsening ratio is within reasonable bounds (0.1 < ratio < 0.8)
        // 2. All fine points are either C-points or strongly connected to a C-point
        // 3. No isolated points (unless the matrix has them)

        let n = 8;
        let mut matrix = DMatrix::zeros(n * n, n * n);
        // Create 8x8 Laplacian
        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                matrix[(idx, idx)] = 4.0;
                if i > 0 {
                    matrix[(idx, (i - 1) * n + j)] = -1.0;
                }
                if i < n - 1 {
                    matrix[(idx, (i + 1) * n + j)] = -1.0;
                }
                if j > 0 {
                    matrix[(idx, i * n + (j - 1))] = -1.0;
                }
                if j < n - 1 {
                    matrix[(idx, i * n + (j + 1))] = -1.0;
                }
            }
        }

        let sparse_matrix = CsrMatrix::from(&matrix);
        let result = ruge_stueben_coarsening(&sparse_matrix, 0.25).unwrap();
        let quality = analyze_coarsening_quality(&result);
        let coarsening_ratio = quality.coarsening_ratio;
        let assignment_ratio = quality.assignment_ratio;
        let avg_interpolation_points = quality.avg_interpolation_points;

        // Check 1: Coarsening ratio
        assert!(
            coarsening_ratio > 0.1,
            "Coarsening too aggressive: {coarsening_ratio}"
        );
        assert!(
            coarsening_ratio < 0.8,
            "Coarsening too weak: {coarsening_ratio}"
        );

        // Check 2: Assignment ratio (should be 1.0 for this connected problem)
        assert!(
            assignment_ratio > 0.95,
            "Not all points assigned to coarse grid: {assignment_ratio}"
        );

        // Check 3: Interpolation points
        // Every F-point should interpolate from at least one C-point
        // (Average should be > 0)
        assert!(
            avg_interpolation_points > 0.0,
            "F-points have no interpolation sources: {avg_interpolation_points}"
        );
    }
}
