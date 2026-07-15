//! Restriction operators for AMG multigrid methods

use leto::{Array1, Array2};
use leto_ops::MatrixProduct;

type RestrictionMatrix = Array2<f64>;
type RestrictionVector = Array1<f64>;

#[inline]
fn matrix_entry(matrix: &RestrictionMatrix, row: usize, col: usize) -> f64 {
    *matrix
        .get([row, col])
        .expect("invariant: restriction matrix index is in bounds")
}

#[inline]
fn vector_entry(vector: &RestrictionVector, row: usize) -> f64 {
    *vector
        .get([row])
        .expect("invariant: restriction vector index is in bounds")
}

/// Create restriction operator from interpolation operator
///
/// The restriction operator is typically the transpose of the interpolation operator.
/// This ensures that the Galerkin condition (coarse_matrix = R * fine_matrix * P) is satisfied.
pub fn create_restriction_from_interpolation(
    interpolation: &RestrictionMatrix,
) -> RestrictionMatrix {
    interpolation
        .transpose([1, 0])
        .expect("invariant: rank-2 transpose axes are valid")
        .to_contiguous()
}

/// Create injection restriction operator
///
/// Simple injection: each fine point contributes only to its corresponding coarse point.
/// This is the transpose of direct interpolation.
pub fn create_injection_restriction(
    fine_to_coarse_map: &[Option<usize>],
    fine_n: usize,
    coarse_n: usize,
) -> RestrictionMatrix {
    let mut restriction = RestrictionMatrix::zeros([coarse_n, fine_n]);

    for (fine_i, &coarse_opt) in fine_to_coarse_map.iter().enumerate() {
        if let Some(coarse_i) = coarse_opt {
            *restriction
                .get_mut([coarse_i, fine_i])
                .expect("invariant: fine-to-coarse map points inside restriction shape") = 1.0;
        }
    }

    restriction
}

/// Create full weighting restriction operator
///
/// Full weighting averages all fine values in a neighborhood.
/// This is commonly used for geometric multigrid.
pub fn create_full_weighting_restriction(fine_n: usize, coarse_n: usize) -> RestrictionMatrix {
    // Assuming 1D grid for simplicity
    // In practice, this would depend on the grid topology
    let mut restriction = RestrictionMatrix::zeros([coarse_n, fine_n]);

    for coarse_i in 0..coarse_n {
        let fine_start = coarse_i * 2;
        let fine_end = (fine_start + 2).min(fine_n);

        let weight = 1.0 / (fine_end - fine_start) as f64;

        for fine_i in fine_start..fine_end {
            *restriction
                .get_mut([coarse_i, fine_i])
                .expect("invariant: full-weighting index is inside restriction shape") = weight;
        }
    }

    restriction
}

/// Create half weighting restriction operator
///
/// Half weighting uses a weighted average with emphasis on direct neighbors.
pub fn create_half_weighting_restriction(fine_n: usize, coarse_n: usize) -> RestrictionMatrix {
    // Assuming 1D grid for simplicity
    let mut restriction = RestrictionMatrix::zeros([coarse_n, fine_n]);

    for coarse_i in 0..coarse_n {
        let fine_center = coarse_i * 2;

        // Weight direct point and neighbors
        *restriction
            .get_mut([coarse_i, fine_center])
            .expect("invariant: half-weighting center index is inside restriction shape") = 0.5;

        if fine_center > 0 {
            *restriction
                .get_mut([coarse_i, fine_center - 1])
                .expect("invariant: half-weighting left index is inside restriction shape") = 0.25;
        }

        if fine_center + 1 < fine_n {
            *restriction
                .get_mut([coarse_i, fine_center + 1])
                .expect("invariant: half-weighting right index is inside restriction shape") = 0.25;
        }
    }

    restriction
}

/// Validate restriction operator properties
pub fn validate_restriction_operator(
    restriction: &RestrictionMatrix,
    interpolation: &RestrictionMatrix,
) -> RestrictionQuality {
    let [coarse_n, fine_n] = restriction.shape();

    // Check that restriction is transpose of interpolation (within tolerance)
    let transpose_error: f64 = (0..coarse_n)
        .flat_map(|row| {
            (0..fine_n).map(move |col| {
                (matrix_entry(restriction, row, col) - matrix_entry(interpolation, col, row)).abs()
            })
        })
        .sum::<f64>()
        / (coarse_n * fine_n) as f64;

    // Check row sums (should be reasonable)
    let mut row_sums = Vec::new();
    for i in 0..coarse_n {
        let row_sum: f64 = (0..fine_n).map(|j| matrix_entry(restriction, i, j)).sum();
        row_sums.push(row_sum);
    }

    let avg_row_sum = row_sums.iter().sum::<f64>() / row_sums.len() as f64;
    let max_row_sum = row_sums.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let min_row_sum = row_sums.iter().copied().fold(f64::INFINITY, f64::min);

    // Check column sums (related to interpolation properties)
    let mut col_sums = Vec::new();
    for j in 0..fine_n {
        let col_sum: f64 = (0..coarse_n).map(|i| matrix_entry(restriction, i, j)).sum();
        col_sums.push(col_sum);
    }

    let avg_col_sum = col_sums.iter().sum::<f64>() / col_sums.len() as f64;

    // Sparsity metrics
    let total_entries = coarse_n * fine_n;
    let non_zero_entries = restriction.iter().filter(|&&x| x.abs() > 1e-12).count();
    let sparsity_ratio = non_zero_entries as f64 / total_entries as f64;

    RestrictionQuality {
        transpose_error,
        avg_row_sum,
        max_row_sum,
        min_row_sum,
        avg_col_sum,
        sparsity_ratio,
        non_zero_entries,
        total_entries,
    }
}

/// Quality metrics for restriction operators
#[derive(Debug, Clone)]
pub struct RestrictionQuality {
    /// Error in transpose relationship with interpolation (should be ~0)
    pub transpose_error: f64,
    /// Average row sum
    pub avg_row_sum: f64,
    /// Maximum row sum
    pub max_row_sum: f64,
    /// Minimum row sum
    pub min_row_sum: f64,
    /// Average column sum
    pub avg_col_sum: f64,
    /// Ratio of non-zero entries to total entries
    pub sparsity_ratio: f64,
    /// Number of non-zero entries
    pub non_zero_entries: usize,
    /// Total number of entries
    pub total_entries: usize,
}

impl RestrictionQuality {
    /// Check if restriction quality is acceptable
    pub fn is_acceptable(&self) -> bool {
        self.transpose_error < 1e-10 && self.avg_row_sum > 0.0 && self.sparsity_ratio < 0.1
        // Less than 10% non-zeros
    }
}

/// Apply restriction operator to a vector
pub fn restrict_vector(
    restriction: &RestrictionMatrix,
    fine_vector: &RestrictionVector,
) -> RestrictionVector {
    let [coarse_n, fine_n] = restriction.shape();
    Array1::from_shape_fn([coarse_n], |[row]| {
        (0..fine_n)
            .map(|col| matrix_entry(restriction, row, col) * vector_entry(fine_vector, col))
            .sum()
    })
}

/// Apply restriction operator to a matrix (Galerkin projection)
pub fn restrict_matrix(
    restriction: &RestrictionMatrix,
    fine_matrix: &RestrictionMatrix,
    interpolation: &RestrictionMatrix,
) -> Result<RestrictionMatrix, &'static str> {
    // Coarse matrix = R * A_fine * P
    let restricted = restriction
        .matmul(fine_matrix)
        .map_err(|_| "restriction/fine matrix shape mismatch")?;
    restricted
        .matmul(interpolation)
        .map_err(|_| "restricted/interpolation matrix shape mismatch")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_interpolation() -> RestrictionMatrix {
        // Create a simple interpolation matrix
        let mut interpolation = RestrictionMatrix::zeros([5, 3]);
        *interpolation.get_mut([0, 0]).unwrap() = 1.0; // C-point
        *interpolation.get_mut([1, 0]).unwrap() = 0.5; // F-point interpolated from C0
        *interpolation.get_mut([1, 1]).unwrap() = 0.5; // F-point interpolated from C1
        *interpolation.get_mut([2, 1]).unwrap() = 1.0; // C-point
        *interpolation.get_mut([3, 1]).unwrap() = 0.3; // F-point
        *interpolation.get_mut([3, 2]).unwrap() = 0.7; // F-point
        *interpolation.get_mut([4, 2]).unwrap() = 1.0; // C-point

        interpolation
    }

    #[test]
    fn test_restriction_from_interpolation() {
        let interpolation = create_test_interpolation();
        let restriction = create_restriction_from_interpolation(&interpolation);

        // Check dimensions
        assert_eq!(restriction.shape(), [3, 5]);
        assert_eq!(interpolation.shape(), [5, 3]);

        // Check transpose property
        for row in 0..restriction.shape()[0] {
            for col in 0..restriction.shape()[1] {
                assert!(
                    (matrix_entry(&restriction, row, col) - matrix_entry(&interpolation, col, row))
                        .abs()
                        < 1e-15
                );
            }
        }
    }

    #[test]
    fn test_injection_restriction() {
        let fine_to_coarse_map = vec![Some(0), None, Some(1), None, Some(2)];

        let restriction = create_injection_restriction(&fine_to_coarse_map, 5, 3);

        // Check dimensions
        assert_eq!(restriction.shape(), [3, 5]);

        // Check that only mapped points contribute
        assert_eq!(matrix_entry(&restriction, 0, 0), 1.0); // Fine 0 -> coarse 0
        assert_eq!(matrix_entry(&restriction, 1, 2), 1.0); // Fine 2 -> coarse 1
        assert_eq!(matrix_entry(&restriction, 2, 4), 1.0); // Fine 4 -> coarse 2

        // Check that unmapped points don't contribute
        assert_eq!(matrix_entry(&restriction, 0, 1), 0.0); // Fine 1 not mapped
        assert_eq!(matrix_entry(&restriction, 1, 3), 0.0); // Fine 3 not mapped
    }

    #[test]
    fn test_restriction_quality_validation() {
        let interpolation = create_test_interpolation();
        let restriction = create_restriction_from_interpolation(&interpolation);

        let quality = validate_restriction_operator(&restriction, &interpolation);

        // Should have very low transpose error
        assert!(quality.transpose_error < 1e-10);

        // Should have reasonable row/column sums
        assert!(quality.avg_row_sum > 0.0);
        assert!(quality.avg_col_sum > 0.0);

        // Should be sparse
        assert!(quality.sparsity_ratio < 1.0);
    }

    #[test]
    fn test_restrict_vector() {
        let interpolation = create_test_interpolation();
        let restriction = create_restriction_from_interpolation(&interpolation);

        let fine_vector = RestrictionVector::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0])
            .expect("valid test vector");
        let coarse_vector = restrict_vector(&restriction, &fine_vector);

        // Check dimensions
        assert_eq!(coarse_vector.shape(), [restriction.shape()[0]]);

        assert_eq!(vector_entry(&coarse_vector, 0), 2.0);
        assert_eq!(vector_entry(&coarse_vector, 1), 5.2);
        assert_eq!(vector_entry(&coarse_vector, 2), 7.8);
    }

    #[test]
    fn test_full_weighting_restriction() {
        let restriction = create_full_weighting_restriction(8, 4);

        // Check dimensions
        assert_eq!(restriction.shape(), [4, 8]);

        // Check that each coarse point gets contributions from 2 fine points
        for i in 0..4 {
            let row_sum: f64 = (0..8).map(|j| matrix_entry(&restriction, i, j)).sum();
            assert!((row_sum - 1.0).abs() < 1e-10); // Should sum to 1
        }
    }

    #[test]
    fn test_half_weighting_restriction() {
        let restriction = create_half_weighting_restriction(8, 4);

        // Check dimensions
        assert_eq!(restriction.shape(), [4, 8]);

        // Check that row sums are reasonable
        for i in 0..4 {
            let row_sum: f64 = (0..8).map(|j| matrix_entry(&restriction, i, j)).sum();
            assert!(row_sum > 0.0 && row_sum <= 1.0);
        }
    }

    #[test]
    fn test_restrict_matrix_uses_leto_galerkin_product() {
        let interpolation = create_test_interpolation();
        let restriction = create_restriction_from_interpolation(&interpolation);
        let fine_matrix = RestrictionMatrix::from_shape_fn([5, 5], |[row, col]| {
            if row == col {
                2.0
            } else if row.abs_diff(col) == 1 {
                -1.0
            } else {
                0.0
            }
        });

        let coarse_matrix = restrict_matrix(&restriction, &fine_matrix, &interpolation)
            .expect("Galerkin projection dimensions are compatible");

        assert_eq!(coarse_matrix.shape(), [3, 3]);
        assert_eq!(matrix_entry(&coarse_matrix, 0, 0), 1.5);
        assert_eq!(matrix_entry(&coarse_matrix, 0, 1), -0.5);
        assert!((matrix_entry(&coarse_matrix, 1, 1) - 1.08).abs() < 1e-12);
        assert!((matrix_entry(&coarse_matrix, 1, 2) + 0.58).abs() < 1e-12);
        assert!((matrix_entry(&coarse_matrix, 2, 2) - 1.58).abs() < 1e-12);
    }
}
