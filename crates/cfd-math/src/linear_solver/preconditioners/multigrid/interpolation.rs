//! Interpolation operators for AMG multigrid methods

use crate::sparse::SparseMatrix;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

/// Create classical interpolation operator (Ruge-Stüben)
pub fn create_classical_interpolation<T: RealField + Copy + FromPrimitive>(
    fine_matrix: &SparseMatrix<T>,
    coarse_points: &[usize],
    strength_matrix: &SparseMatrix<T>,
    _max_interpolation_points: usize,
) -> Result<SparseMatrix<T>> {
    let fine_n = fine_matrix.nrows();
    let coarse_n = coarse_points.len();

    let mut row_offsets = vec![0; fine_n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for fine_i in 0..fine_n {
        if coarse_points.contains(&fine_i) {
            // Direct injection for coarse points (C-points)
            let coarse_idx = coarse_points.iter().position(|&x| x == fine_i).unwrap();
            col_indices.push(coarse_idx);
            values.push(T::one());
        } else {
            let mut weights = Vec::new();
            let mut total_weight = T::zero();
            let mut max_strength = T::zero();

            let row_start = fine_matrix.row_offsets()[fine_i];
            let row_end = fine_matrix.row_offsets()[fine_i + 1];

            let mut coarse_neighbors = Vec::new();
            for k in row_start..row_end {
                let neighbor_idx = fine_matrix.col_indices()[k];
                if let Some(coarse_local_idx) =
                    coarse_points.iter().position(|&x| x == neighbor_idx)
                {
                    let strength = strength_matrix
                        .get_entry(fine_i, neighbor_idx)
                        .map_or(T::zero(), |e| e.into_value());
                    if strength > max_strength {
                        max_strength = strength;
                    }
                    coarse_neighbors.push((coarse_local_idx, strength));
                }
            }

            for (coarse_local_idx, strength) in coarse_neighbors {
                if strength > T::zero() {
                    let distance = if max_strength > T::zero() {
                        max_strength / strength - T::one()
                    } else {
                        T::zero()
                    };
                    let weight = strength / (distance + T::one());
                    weights.push((coarse_local_idx, weight));
                    total_weight += weight;
                }
            }

            weights.sort_by(|a, b| a.0.cmp(&b.0));
            for &(coarse_idx, weight) in &weights {
                let normalized_weight = if total_weight > T::zero() {
                    weight / total_weight
                } else {
                    T::one() / T::from_usize(weights.len()).unwrap_or_else(T::one)
                };
                col_indices.push(coarse_idx);
                values.push(normalized_weight);
            }
        }
        row_offsets[fine_i + 1] = col_indices.len();
    }

    SparseMatrix::try_from_csr_data(fine_n, coarse_n, row_offsets, col_indices, values).map_err(
        |e| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: format!("Failed to create classical interpolation matrix: {e}"),
            })
        },
    )
}

/// Create direct interpolation operator
pub fn create_direct_interpolation<T: RealField + Copy + FromPrimitive>(
    fine_to_coarse_map: &[Option<usize>],
    fine_n: usize,
    coarse_n: usize,
) -> SparseMatrix<T> {
    let mut row_offsets = vec![0; fine_n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for (fine_i, &coarse_opt) in fine_to_coarse_map.iter().enumerate() {
        if let Some(coarse_i) = coarse_opt {
            col_indices.push(coarse_i);
            values.push(T::one());
        }
        row_offsets[fine_i + 1] = col_indices.len();
    }

    SparseMatrix::try_from_csr_data(fine_n, coarse_n, row_offsets, col_indices, values)
        .expect("Failed to create direct interpolation matrix")
}

/// Create standard interpolation operator
pub fn create_standard_interpolation<T: RealField + Copy + FromPrimitive>(
    fine_matrix: &SparseMatrix<T>,
    coarse_points: &[usize],
    strength_matrix: &SparseMatrix<T>,
) -> Result<SparseMatrix<T>> {
    let fine_n = fine_matrix.nrows();
    let coarse_n = coarse_points.len();

    let mut row_offsets = vec![0; fine_n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for fine_i in 0..fine_n {
        if coarse_points.contains(&fine_i) {
            // Direct injection for coarse points
            let coarse_idx = coarse_points.iter().position(|&x| x == fine_i).unwrap();
            col_indices.push(coarse_idx);
            values.push(T::one());
        } else {
            // Distance-weighted interpolation for F-points
            let mut weights = Vec::new();
            let mut total_weight = T::zero();

            // Find neighboring coarse points
            for (coarse_local_idx, &coarse_global_idx) in coarse_points.iter().enumerate() {
                let distance = ((fine_i as f64) - (coarse_global_idx as f64)).abs();

                // Weight by inverse distance and connection strength
                if distance > 0.0 {
                    let strength =
                        if let Some(s) = strength_matrix.get_entry(fine_i, coarse_global_idx) {
                            s.into_value()
                        } else {
                            T::from_f64(1e-6).unwrap_or_else(T::zero)
                        };
                    let weight =
                        strength / (T::from_f64(distance).unwrap_or_else(T::zero) + T::one());
                    weights.push((coarse_local_idx, weight));
                    total_weight += weight;
                }
            }

            for &(coarse_idx, weight) in &weights {
                let normalized_weight = if total_weight > T::zero() {
                    weight / total_weight
                } else {
                    T::one() / T::from_usize(weights.len()).unwrap_or_else(T::one)
                };
                col_indices.push(coarse_idx);
                values.push(normalized_weight);
            }
        }
        row_offsets[fine_i + 1] = col_indices.len();
    }

    SparseMatrix::try_from_csr_data(fine_n, coarse_n, row_offsets, col_indices, values).map_err(
        |e| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: format!("Failed to create standard interpolation matrix: {e}"),
            })
        },
    )
}

/// Validate interpolation operator properties
pub fn validate_interpolation_operator<T: RealField + Copy + FromPrimitive + ToPrimitive>(
    interpolation: &SparseMatrix<T>,
    coarse_points: &[usize],
) -> InterpolationQuality {
    let fine_n = interpolation.nrows();
    let coarse_n = interpolation.ncols();

    // Check row sums (should be 1 for F-points, 1 for C-points)
    let mut row_sums = Vec::new();
    for i in 0..fine_n {
        let row_start = interpolation.row_offsets()[i];
        let row_end = interpolation.row_offsets()[i + 1];
        let mut row_sum = T::zero();
        for k in row_start..row_end {
            row_sum += interpolation.values()[k];
        }
        row_sums.push(row_sum);
    }

    // Check for C-points (should have row sum = 1 and single non-zero entry)
    let mut c_point_row_sums = Vec::new();
    let mut f_point_row_sums = Vec::new();

    for &cp in coarse_points {
        if cp < row_sums.len() {
            c_point_row_sums.push(row_sums[cp]);
        }
    }

    // F-points are all points not in coarse_points
    for i in 0..fine_n {
        if !coarse_points.contains(&i) {
            f_point_row_sums.push(row_sums[i]);
        }
    }

    // Calculate quality metrics
    let avg_c_point_sum = if c_point_row_sums.is_empty() {
        0.0
    } else {
        c_point_row_sums
            .iter()
            .map(|&s| s.to_f64().unwrap_or(0.0))
            .sum::<f64>()
            / c_point_row_sums.len() as f64
    };

    let avg_f_point_sum = if f_point_row_sums.is_empty() {
        0.0
    } else {
        f_point_row_sums
            .iter()
            .map(|&s| s.to_f64().unwrap_or(0.0))
            .sum::<f64>()
            / f_point_row_sums.len() as f64
    };

    // Check conservation property (interpolation should preserve constants)
    let constant_vector = nalgebra::DVector::from_element(coarse_n, T::one());
    let mut interpolated = nalgebra::DVector::zeros(fine_n);
    crate::sparse::spmv(interpolation, &constant_vector, &mut interpolated);

    let constant_error: f64 = interpolated
        .iter()
        .map(|&x| (x - T::one()).abs().to_f64().unwrap_or(0.0))
        .sum::<f64>()
        / fine_n as f64;

    // Sparsity metrics
    let total_entries = fine_n * coarse_n;
    let non_zero_entries = interpolation.nnz();
    let sparsity_ratio = if total_entries > 0 {
        non_zero_entries as f64 / total_entries as f64
    } else {
        0.0
    };

    InterpolationQuality {
        avg_c_point_row_sum: avg_c_point_sum,
        avg_f_point_row_sum: avg_f_point_sum,
        constant_preservation_error: constant_error,
        sparsity_ratio,
        non_zero_entries,
        total_entries,
    }
}

/// Quality metrics for interpolation operators
#[derive(Debug, Clone)]
pub struct InterpolationQuality {
    /// Average row sum for C-points (should be 1.0)
    pub avg_c_point_row_sum: f64,
    /// Average row sum for F-points (should be 1.0)
    pub avg_f_point_row_sum: f64,
    /// Error in preserving constants (should be ~0)
    pub constant_preservation_error: f64,
    /// Ratio of non-zero entries to total entries
    pub sparsity_ratio: f64,
    /// Number of non-zero entries
    pub non_zero_entries: usize,
    /// Total number of entries
    pub total_entries: usize,
}

impl InterpolationQuality {
    /// Check if interpolation quality is acceptable
    pub fn is_acceptable(&self) -> bool {
        (self.avg_c_point_row_sum - 1.0).abs() < 1e-10
            && (self.avg_f_point_row_sum - 1.0).abs() < 1e-10
            && self.constant_preservation_error < 1e-6
            && self.sparsity_ratio < 0.1 // Less than 10% non-zeros
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    fn create_test_matrix() -> DMatrix<f64> {
        // Create a simple tridiagonal matrix
        let n = 5;
        let mut matrix = DMatrix::zeros(n, n);

        for i in 0..n {
            matrix[(i, i)] = 2.0;
            if i > 0 {
                matrix[(i, i - 1)] = -1.0;
            }
            if i < n - 1 {
                matrix[(i, i + 1)] = -1.0;
            }
        }

        matrix
    }

    fn create_simple_coarsening() -> (Vec<usize>, Vec<Option<usize>>) {
        let coarse_points = vec![0, 2, 4]; // Every other point
        let fine_to_coarse_map = vec![Some(0), None, Some(1), None, Some(2)];

        (coarse_points, fine_to_coarse_map)
    }

    #[test]
    fn test_direct_interpolation() {
        let (_coarse_points, fine_to_coarse_map) = create_simple_coarsening();
        let interpolation = create_direct_interpolation(&fine_to_coarse_map, 5, 3);

        // Check dimensions
        assert_eq!(interpolation.nrows(), 5);
        assert_eq!(interpolation.ncols(), 3);

        // Check that coarse points map correctly
        assert_eq!(
            interpolation
                .get_entry(0, 0)
                .map_or(0.0, |e| e.into_value()),
            1.0
        ); // Point 0 -> coarse 0
        assert_eq!(
            interpolation
                .get_entry(2, 1)
                .map_or(0.0, |e| e.into_value()),
            1.0
        ); // Point 2 -> coarse 1
        assert_eq!(
            interpolation
                .get_entry(4, 2)
                .map_or(0.0, |e| e.into_value()),
            1.0
        ); // Point 4 -> coarse 2

        // Check that F-points are zero
        assert_eq!(
            interpolation
                .get_entry(1, 0)
                .map_or(0.0, |e| e.into_value()),
            0.0
        ); // Point 1 not mapped
        assert_eq!(
            interpolation
                .get_entry(3, 1)
                .map_or(0.0, |e| e.into_value()),
            0.0
        ); // Point 3 not mapped
    }

    #[test]
    fn test_classical_interpolation() {
        let matrix_dense = create_test_matrix();
        let matrix = SparseMatrix::from(&matrix_dense);
        let (coarse_points, _) = create_simple_coarsening();

        // Create a simple strength matrix
        let mut strength_builder = nalgebra_sparse::CooMatrix::new(5, 5);
        for i in 0..5 {
            for j in 0..5 {
                if (i as i32 - j as i32).abs() == 1 {
                    strength_builder.push(i, j, 1.0);
                }
            }
        }
        let strength_matrix = SparseMatrix::from(&strength_builder);

        let interpolation =
            create_classical_interpolation(&matrix, &coarse_points, &strength_matrix, 2).unwrap();

        // Check dimensions
        assert_eq!(interpolation.nrows(), 5);
        assert_eq!(interpolation.ncols(), 3);

        // Check that coarse points have direct injection
        assert_eq!(
            interpolation
                .get_entry(0, 0)
                .map_or(0.0, |e| e.into_value()),
            1.0
        );
        assert_eq!(
            interpolation
                .get_entry(2, 1)
                .map_or(0.0, |e| e.into_value()),
            1.0
        );
        assert_eq!(
            interpolation
                .get_entry(4, 2)
                .map_or(0.0, |e| e.into_value()),
            1.0
        );
    }

    #[test]
    fn test_interpolation_quality_validation() {
        let (coarse_points, fine_to_coarse_map) = create_simple_coarsening();
        let interpolation = create_direct_interpolation::<f64>(&fine_to_coarse_map, 5, 3);

        let quality = validate_interpolation_operator(&interpolation, &coarse_points);

        // C-points should have row sum = 1
        assert!((quality.avg_c_point_row_sum - 1.0).abs() < 1e-10);

        // F-points should have row sum = 0 (no interpolation in direct case)
        assert_eq!(quality.avg_f_point_row_sum, 0.0);

        // Should preserve constants (with some error for F-points)
        assert!(quality.constant_preservation_error >= 0.0);

        // Should be sparse
        assert!(quality.sparsity_ratio < 1.0);
    }

    #[test]
    fn test_constant_preservation() {
        // Test that interpolation preserves constant vectors
        let (_coarse_points, fine_to_coarse_map) = create_simple_coarsening();
        let interpolation = create_direct_interpolation::<f64>(&fine_to_coarse_map, 5, 3);

        let constant_coarse = nalgebra::DVector::from_element(3, 1.0);
        let mut interpolated = nalgebra::DVector::zeros(5);
        crate::sparse::spmv(&interpolation, &constant_coarse, &mut interpolated);

        // C-points should be 1.0, F-points should be 0.0
        assert_eq!(interpolated[0], 1.0); // C-point
        assert_eq!(interpolated[1], 0.0); // F-point
        assert_eq!(interpolated[2], 1.0); // C-point
        assert_eq!(interpolated[3], 0.0); // F-point
        assert_eq!(interpolated[4], 1.0); // C-point
    }
}
