//! Coarsening algorithms for AMG hierarchy construction.
//!
//! Provides Ruge-Stüben, aggregation, Falgout (CLJP), PMIS, HMIS, and hybrid
//! coarsening strategies.

use super::CoarseningResult;
use crate::SparseMatrix;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use rand::Rng;

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
            if max_lambda <= 0 {
                break;
            }

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

    // Second Pass - Classify remaining undecided points as C-points
    for i in 0..n {
        if status[i] == 0 {
            status[i] = 1;
            coarse_points.push(i);
            fine_to_coarse_map[i] = Some(coarse_points.len() - 1);
        }
    }

    // Step 4: Map F-points to their strongest connected C-point
    assign_closest_coarse_points(&mut fine_to_coarse_map, &status, &strength_matrix);

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

    let strength_matrix =
        compute_strength_matrix(matrix, T::from_f64(0.5).unwrap_or_else(T::zero))?;
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    let mut aggregate_id = 0;

    for i in 0..n {
        if !aggregated[i] {
            let mut aggregate = vec![i];
            aggregated[i] = true;

            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                if !aggregated[j] && aggregate.len() < max_aggregate_size {
                    aggregate.push(j);
                    aggregated[j] = true;
                }
            }

            coarse_points.push(aggregate[0]);

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
    match ruge_stueben_coarsening(matrix, strength_threshold) {
        Ok(result) => {
            let assigned_points = result
                .fine_to_coarse_map
                .iter()
                .filter(|x| x.is_some())
                .count();
            let assignment_ratio = assigned_points as f64 / matrix.nrows() as f64;

            if assignment_ratio > 0.8 {
                Ok(result)
            } else {
                aggregation_coarsening(matrix, max_aggregate_size)
            }
        }
        Err(_) => aggregation_coarsening(matrix, max_aggregate_size),
    }
}

/// Assign unmapped F-points to their strongest connected C-point
fn assign_closest_coarse_points<T: RealField + Copy + FromPrimitive>(
    fine_to_coarse_map: &mut [Option<usize>],
    status: &[i32],
    strength_matrix: &SparseMatrix<T>,
) {
    let n = strength_matrix.nrows();
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    for i in 0..n {
        if fine_to_coarse_map[i].is_none() {
            let mut max_strength = T::from_f64(-1.0).unwrap_or_else(T::zero);
            let mut best_coarse_idx = None;

            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                if status[j] == 1 {
                    let strength = strength_matrix.values()[k];
                    if strength > max_strength {
                        max_strength = strength;
                        best_coarse_idx = fine_to_coarse_map[j];
                    }
                }
            }

            if let Some(idx) = best_coarse_idx {
                fine_to_coarse_map[i] = Some(idx);
            }
        }
    }
}

/// Falgout coarsening algorithm (CLJP method)
///
/// Reference: Falgout, R. D. (2006). An introduction to algebraic multigrid
pub fn falgout_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();
    let mut coarse_points = Vec::new();
    let mut fine_to_coarse_map = vec![None; n];

    let strength_matrix = compute_strength_matrix(matrix, strength_threshold)?;
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    let mut measures = vec![0.0; n];
    for i in 0..n {
        measures[i] = (s_offsets[i + 1] - s_offsets[i]) as f64;
    }

    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_by(|&a, &b| measures[b].partial_cmp(&measures[a]).unwrap());

    let lambda = 4.0 / 3.0;
    let mut status = vec![0; n];

    for &i in &sorted_indices {
        if status[i] != 0 {
            continue;
        }

        let mut should_be_coarse = true;
        let mut coarse_neighbors = 0;
        let total_strong_connections = s_offsets[i + 1] - s_offsets[i];

        for k in s_offsets[i]..s_offsets[i + 1] {
            let j = s_indices[k];
            if status[j] == 1 {
                coarse_neighbors += 1;
            }
        }

        if total_strong_connections > 0 {
            let ratio = f64::from(coarse_neighbors) / total_strong_connections as f64;
            if ratio >= lambda {
                should_be_coarse = false;
            }
        }

        if should_be_coarse {
            status[i] = 1;
            coarse_points.push(i);
            fine_to_coarse_map[i] = Some(coarse_points.len() - 1);

            for k in s_offsets[i]..s_offsets[i + 1] {
                let j = s_indices[k];
                if status[j] == 0 {
                    status[j] = 2;
                    fine_to_coarse_map[j] = Some(coarse_points.len() - 1);
                }
            }
        } else {
            status[i] = 2;
        }
    }

    assign_closest_coarse_points(&mut fine_to_coarse_map, &status, &strength_matrix);

    Ok(CoarseningResult {
        coarse_points,
        fine_to_coarse_map,
        strength_matrix,
    })
}

/// PMIS (Parallel Modified Independent Set) coarsening
///
/// Reference: Luby's algorithm adapted for parallel coarsening
pub fn pmis_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();
    let mut coarse_points = Vec::new();
    let mut fine_to_coarse_map = vec![None; n];

    let strength_matrix = compute_strength_matrix(matrix, strength_threshold)?;
    let s_offsets = strength_matrix.row_offsets();
    let s_indices = strength_matrix.col_indices();

    let mut status = vec![0; n];

    let mut rng = rand::thread_rng();
    let priorities: Vec<f64> = (0..n).map(|_| rng.gen()).collect();

    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_by(|&a, &b| priorities[b].partial_cmp(&priorities[a]).unwrap());

    for &i in &sorted_indices {
        if status[i] != 0 {
            continue;
        }

        status[i] = 1;
        coarse_points.push(i);
        fine_to_coarse_map[i] = Some(coarse_points.len() - 1);

        for k in s_offsets[i]..s_offsets[i + 1] {
            let j = s_indices[k];
            if status[j] == 0 {
                status[j] = 2;
            }
        }
    }

    assign_closest_coarse_points(&mut fine_to_coarse_map, &status, &strength_matrix);

    Ok(CoarseningResult {
        coarse_points,
        fine_to_coarse_map,
        strength_matrix,
    })
}

/// HMIS (Hybrid Modified Independent Set) coarsening
///
/// Combines PMIS with aggressive coarsening for better parallel performance.
pub fn hmis_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
    aggressive_threshold: T,
) -> Result<CoarseningResult<T>> {
    let n = matrix.nrows();

    match pmis_coarsening(matrix, strength_threshold) {
        Ok(result) => {
            let coarsening_ratio =
                T::from_f64(result.coarse_points.len() as f64 / n as f64).unwrap_or_else(T::zero);

            if coarsening_ratio < aggressive_threshold {
                aggressive_coarsening(matrix, strength_threshold, aggressive_threshold)
            } else {
                Ok(result)
            }
        }
        Err(_) => aggressive_coarsening(matrix, strength_threshold, aggressive_threshold),
    }
}

/// Aggressive coarsening strategy for difficult matrices
fn aggressive_coarsening<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    _strength_threshold: T,
    _target_ratio: T,
) -> Result<CoarseningResult<T>> {
    let max_aggregate_size = 8;
    aggregation_coarsening(matrix, max_aggregate_size)
}

/// Compute strength of connection matrix
pub(super) fn compute_strength_matrix<T: RealField + Copy + FromPrimitive>(
    matrix: &SparseMatrix<T>,
    strength_threshold: T,
) -> Result<SparseMatrix<T>> {
    use cfd_core::error::{Error, NumericalErrorKind};

    let n = matrix.nrows();
    let mut row_offsets = vec![0; n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    let m_offsets = matrix.row_offsets();
    let m_indices = matrix.col_indices();
    let m_values = matrix.values();

    for i in 0..n {
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
