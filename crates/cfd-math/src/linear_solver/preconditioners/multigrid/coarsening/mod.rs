//! Coarsening strategies for AMG hierarchy construction.
//!
//! Provides multiple coarsening algorithms (Ruge-Stüben, aggregation, Falgout,
//! PMIS, HMIS) and quality analysis via algebraic distance measures.

mod algorithms;
mod quality;

pub use algorithms::*;
pub use quality::*;

use crate::SparseMatrix;
use nalgebra::RealField;

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

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    use nalgebra_sparse::CsrMatrix;

    fn create_test_matrix() -> CsrMatrix<f64> {
        let n = 4;
        let mut matrix = DMatrix::zeros(n * n, n * n);

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

        CsrMatrix::from(&matrix)
    }

    #[test]
    fn test_ruge_stueben_coarsening() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();

        assert!(!result.coarse_points.is_empty());
        assert_eq!(result.fine_to_coarse_map.len(), matrix.nrows());
        assert_eq!(result.strength_matrix.nrows(), matrix.nrows());

        let quality = analyze_coarsening_quality(&result, &matrix);
        assert!(quality.coarsening_ratio > 0.0);
        assert!(quality.assignment_ratio > 0.0);
    }

    #[test]
    fn test_aggregation_coarsening() {
        let matrix = create_test_matrix();
        let result = aggregation_coarsening(&matrix, 4).unwrap();

        assert!(!result.coarse_points.is_empty());
        assert_eq!(result.fine_to_coarse_map.len(), matrix.nrows());

        let quality = analyze_coarsening_quality(&result, &matrix);
        assert!(quality.coarsening_ratio > 0.0);
    }

    #[test]
    fn test_coarsening_quality_analysis() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();
        let quality = analyze_coarsening_quality(&result, &matrix);

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
            for (i, &map) in result.fine_to_coarse_map.iter().enumerate() {
                if let Some(coarse_idx) = map {
                    assert!(
                        coarse_idx < result.coarse_points.len(),
                        "[{name}] Mapped index out of bounds"
                    );

                    if result.coarse_points.contains(&i) {
                        assert_eq!(
                            result.coarse_points[coarse_idx], i,
                            "[{name}] Coarse point mapped to wrong index"
                        );
                    }
                }
            }

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
        let n = 10;
        let mut matrix = DMatrix::zeros(n * n, n * n);

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

        let n_total = sparse_matrix.nrows();
        let n_coarse = result.coarse_points.len();
        let ratio = n_coarse as f64 / n_total as f64;

        assert!(ratio >= 0.15, "Coarsening ratio too low: {ratio}");
        assert!(ratio <= 0.65, "Coarsening ratio too high: {ratio}");
    }

    #[test]
    fn test_interpolation_operator_shape() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();

        let n_fine = matrix.nrows();
        let n_coarse = result.coarse_points.len();

        assert_eq!(result.fine_to_coarse_map.len(), n_fine);

        for &map in &result.fine_to_coarse_map {
            if let Some(idx) = map {
                assert!(
                    idx < n_coarse,
                    "Coarse index {idx} out of bounds (n_coarse={n_coarse})"
                );
            }
        }

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
    fn test_algebraic_distances_compute() {
        let matrix = create_test_matrix();
        let result = ruge_stueben_coarsening(&matrix, 0.25).unwrap();
        let distances = AlgebraicDistances::compute(&result, &matrix);

        assert_eq!(distances.distances.len(), matrix.nrows());
        for &cp in &result.coarse_points {
            assert_eq!(distances.distances[cp], 0.0);
        }
        for (i, &d) in distances.distances.iter().enumerate() {
            if !result.coarse_points.contains(&i) {
                assert!(d > 0.0);
                assert!(d < 100.0);
            }
        }
    }

    #[test]
    fn test_coarsening_convergence_behavior() {
        let n = 8;
        let mut matrix = DMatrix::zeros(n * n, n * n);
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
        let quality = analyze_coarsening_quality(&result, &sparse_matrix);
        let coarsening_ratio = quality.coarsening_ratio;
        let assignment_ratio = quality.assignment_ratio;
        let avg_interpolation_points = quality.avg_interpolation_points;

        assert!(
            coarsening_ratio > 0.1,
            "Coarsening too aggressive: {coarsening_ratio}"
        );
        assert!(
            coarsening_ratio < 0.8,
            "Coarsening too weak: {coarsening_ratio}"
        );

        assert!(
            assignment_ratio > 0.95,
            "Not all points assigned to coarse grid: {assignment_ratio}"
        );

        assert!(
            avg_interpolation_points > 0.0,
            "F-points have no interpolation sources: {avg_interpolation_points}"
        );
    }
}
