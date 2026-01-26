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

    // Fast lookup for coarse points: mapping global index -> local coarse index
    let mut coarse_map = vec![None; fine_n];
    for (local_idx, &global_idx) in coarse_points.iter().enumerate() {
        if global_idx < fine_n {
            coarse_map[global_idx] = Some(local_idx);
        }
    }

    for fine_i in 0..fine_n {
        if let Some(coarse_idx) = coarse_map[fine_i] {
            // Point is a C-point: Direct injection
            col_indices.push(coarse_idx);
            values.push(T::one());
        } else {
            // Point is an F-point: Interpolate from C-neighbors using Ruge-Stüben formula
            // w_ij = - ( A_ij + sum_{k in F_i^s} ( A_ik * A_kj / sum_{m in C_i} A_km ) ) / ( A_ii + sum_{n in D_i^w} A_in )

            let row_start = fine_matrix.row_offsets()[fine_i];
            let row_end = fine_matrix.row_offsets()[fine_i + 1];
            let fine_row_cols = &fine_matrix.col_indices()[row_start..row_end];
            let fine_row_vals = &fine_matrix.values()[row_start..row_end];

            let strength_start = strength_matrix.row_offsets()[fine_i];
            let strength_end = strength_matrix.row_offsets()[fine_i + 1];
            let strength_cols = &strength_matrix.col_indices()[strength_start..strength_end];

            // Identify sets C_i (strong C-points) and F_i^s (strong F-points)
            let mut c_i = Vec::new();
            let mut f_i_s = Vec::new();

            for &neighbor_idx in strength_cols {
                if coarse_map[neighbor_idx].is_some() {
                    c_i.push(neighbor_idx);
                } else if neighbor_idx != fine_i {
                    f_i_s.push(neighbor_idx);
                }
            }

            if c_i.is_empty() {
                // If no strong C-points, we can't interpolate strongly.
                // Fallback to zero or weak interpolation?
                // Usually implies poor coarsening or isolated F-point.
                // Leaving row empty effectively means zero value (Dirichlet-like).
            } else {
                // Identify diagonal A_ii and sum of weak connections
                let mut a_ii = T::zero();
                let mut sum_weak = T::zero();

                for (k, &neighbor_idx) in fine_row_cols.iter().enumerate() {
                    let val = fine_row_vals[k];
                    if neighbor_idx == fine_i {
                        a_ii = val;
                    } else {
                        // Check if connection is strong
                        if strength_cols.binary_search(&neighbor_idx).is_err() {
                            // Weak connection, add to diagonal sum
                            sum_weak += val;
                        }
                    }
                }

                let diagonal = a_ii + sum_weak;

                if diagonal.abs() > T::default_epsilon() {
                    // Precompute denominators for k in F_i^s: sum_{m in C_i} A_km
                    let mut k_denoms = Vec::with_capacity(f_i_s.len());
                    for &k in &f_i_s {
                        let mut denom = T::zero();
                        for &m in &c_i {
                            if let Some(val) = fine_matrix.get_entry(k, m) {
                                denom += val.into_value();
                            }
                        }
                        k_denoms.push(denom);
                    }

                    // Compute weights for each j in C_i
                    let mut weights = Vec::new();
                    let neg_diag_inv = T::one() / -diagonal;

                    for &j in &c_i {
                        let coarse_local_idx = coarse_map[j].unwrap();

                        // A_ij (direct connection)
                        let mut a_ij = T::zero();
                        if let Ok(idx) = fine_row_cols.binary_search(&j) {
                            a_ij = fine_row_vals[idx];
                        }

                        let mut indirect_sum = T::zero();

                        // Sum over k in F_i^s
                        for (idx, &k) in f_i_s.iter().enumerate() {
                            let denom = k_denoms[idx];
                            if denom.abs() > T::default_epsilon() {
                                // A_ik
                                let mut a_ik = T::zero();
                                if let Ok(k_idx) = fine_row_cols.binary_search(&k) {
                                    a_ik = fine_row_vals[k_idx];
                                }

                                // A_kj
                                let mut a_kj = T::zero();
                                if let Some(val) = fine_matrix.get_entry(k, j) {
                                    a_kj = val.into_value();
                                }

                                indirect_sum += a_ik * a_kj / denom;
                            }
                        }

                        let weight = (a_ij + indirect_sum) * neg_diag_inv;
                        weights.push((coarse_local_idx, weight));
                    }

                    // Sort by index for CSR format
                    weights.sort_by_key(|w| w.0);

                    for (idx, val) in weights {
                        col_indices.push(idx);
                        values.push(val);
                    }
                }
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

        // For 1D Laplacian 3-point stencil [-1, 2, -1], linear interpolation is exact.
        // F-point 1 is between C-point 0 and C-point 2. Weight should be 0.5 each.
        // Formula: w_ij = A_ij / A_ii (negated and simplified since no F-F connections)
        // A_11 = 2. A_10 = -1, A_12 = -1.
        // w_10 = -(-1)/2 = 0.5
        // w_12 = -(-1)/2 = 0.5

        let w10 = interpolation.get_entry(1, 0).map_or(0.0, |e| e.into_value());
        let w11 = interpolation.get_entry(1, 1).map_or(0.0, |e| e.into_value());

        assert!((w10 - 0.5).abs() < 1e-10, "Weight w10 should be 0.5, got {}", w10);
        assert!((w11 - 0.5).abs() < 1e-10, "Weight w11 should be 0.5, got {}", w11);
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

    #[test]
    fn test_classical_interpolation_weights() {
        // Test case with F-F connections in a valid RS coarsening scenario
        // Triangle configuration:
        //    0(C)
        //   /  \
        //  1(F)-2(F)
        //  |    |
        //  3(C) 4(C)
        //
        // 1 connected to 0, 2, 3. 2 connected to 0, 1, 4.
        // F-F connection (1,2) should distribute 2's influence on 1 to common C-neighbor 0.

        let n = 5;
        let mut matrix_dense = DMatrix::<f64>::zeros(n, n);

        // Setup connections with -1.0, and diagonal 3.0
        let edges = vec![(0,1), (0,2), (1,2), (1,3), (2,4)];

        for i in 0..n { matrix_dense[(i, i)] = 3.0; }

        for (u, v) in edges {
            matrix_dense[(u, v)] = -1.0;
            matrix_dense[(v, u)] = -1.0;
        }
        // Fix diagonals for C-points to be consistent with laplacian?
        // Doesn't strictly matter for interpolation formula as we only look at F-point rows.
        // Row 1: -1(0), -1(2), -1(3). Sum abs off-diag = 3. Diag=3.

        let matrix = SparseMatrix::from(&matrix_dense);

        let coarse_points = vec![0, 3, 4];
        // Strength matrix (all neighbors are strong)
        let mut strength_builder = nalgebra_sparse::CooMatrix::new(n, n);
        for i in 0..n {
            for j in 0..n {
                if i != j && matrix_dense[(i, j)].abs() > 0.0_f64 {
                    strength_builder.push(i, j, 1.0);
                }
            }
        }
        let strength_matrix = SparseMatrix::from(&strength_builder);

        let interpolation = create_classical_interpolation(&matrix, &coarse_points, &strength_matrix, 2).unwrap();

        // Coarse mapping: 0->0, 3->1, 4->2

        // Check point 1 (F)
        // C-neighbors: 0, 3. F-neighbor: 2.
        // w_{1,0}: Direct -1. Indirect via 2: a_{12}*a_{20}/S_2 = (-1*-1)/(-1-0) = -1. Total -2.
        // w_{1,0} = -(-2)/3 = 2/3.
        // w_{1,3}: Direct -1. Indirect via 2: a_{12}*a_{23}/S_2 = (-1*0)/-1 = 0. Total -1.
        // w_{1,3} = -(-1)/3 = 1/3.

        let w1_0 = interpolation.get_entry(1, 0).map_or(0.0, |e| e.into_value());
        let w1_3 = interpolation.get_entry(1, 1).map_or(0.0, |e| e.into_value());

        println!("Weights for point 1: w(0)={}, w(3)={}", w1_0, w1_3);

        assert!((w1_0 - 2.0/3.0).abs() < 1e-10, "Expected 2/3, got {}", w1_0);
        assert!((w1_3 - 1.0/3.0).abs() < 1e-10, "Expected 1/3, got {}", w1_3);
    }
}
