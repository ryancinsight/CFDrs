//! Tests for sparse matrix module

#[cfg(test)]
mod tests {
    use super::super::*;
    use approx::assert_relative_eq;
    use cfd_core::error::Result;

    #[test]
    fn test_matrix_builder() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(3, 3);
        builder.add_triplets(vec![(0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0)])?;

        let matrix = builder.build()?;

        // Check row 0 entries
        let row0 = matrix.row(0);
        let row0_entries: Vec<_> = row0
            .col_indices()
            .iter()
            .zip(row0.values().iter())
            .map(|(&col, &val)| (col, val))
            .collect();
        assert_eq!(row0_entries.len(), 2);
        assert!(row0_entries.contains(&(0, 1.0)));
        assert!(row0_entries.contains(&(2, 2.0)));

        // Check row 1 entries
        let row1 = matrix.row(1);
        let row1_entries: Vec<_> = row1
            .col_indices()
            .iter()
            .zip(row1.values().iter())
            .map(|(&col, &val)| (col, val))
            .collect();
        assert_eq!(row1_entries.len(), 1);
        assert!(row1_entries.contains(&(1, 3.0)));

        Ok(())
    }

    #[test]
    fn test_tridiagonal_pattern() -> Result<()> {
        let matrix = SparsePatterns::tridiagonal(4, -1.0, 2.0, -1.0)?;

        assert_eq!(matrix.nrows(), 4);
        assert_eq!(matrix.ncols(), 4);
        assert_eq!(matrix.nnz(), 10); // 3*4 - 2 = 10

        // Check diagonal values
        let diag = matrix.diagonal();
        for i in 0..4 {
            assert_relative_eq!(diag[i], 2.0, epsilon = 1e-10);
        }

        Ok(())
    }

    #[test]
    fn test_five_point_stencil() -> Result<()> {
        let matrix = SparsePatterns::five_point_stencil(3, 3, 1.0, 1.0)?;

        // Check that it creates a 9x9 matrix (3x3 grid)
        assert_eq!(matrix.nrows(), 9);
        assert_eq!(matrix.ncols(), 9);

        // Center point should have 4 connections + diagonal
        let center = 4; // (1,1) in 3x3 grid
        let center_row = matrix.row(center);
        assert_eq!(center_row.nnz(), 5); // center + 4 neighbors

        Ok(())
    }

    #[test]
    fn test_frobenius_norm() -> Result<()> {
        let matrix = SparsePatterns::tridiagonal(3, 1.0, 2.0, 1.0)?;
        let norm = matrix.frobenius_norm();

        // Expected: sqrt(2^2 * 3 + 1^2 * 4) = sqrt(12 + 4) = 4.0
        assert_relative_eq!(norm, 4.0, epsilon = 1e-10);

        Ok(())
    }

    #[test]
    fn test_diagonal_dominance() -> Result<()> {
        // Create a diagonally dominant matrix
        let dominant = SparsePatterns::tridiagonal(3, 1.0, 4.0, 1.0)?;
        assert!(dominant.is_diagonally_dominant());

        // Create a non-diagonally dominant matrix
        let non_dominant = SparsePatterns::tridiagonal(3, 2.0, 3.0, 2.0)?;
        assert!(!non_dominant.is_diagonally_dominant());

        Ok(())
    }

    #[test]
    fn test_spmv_basic() {
        use nalgebra::DVector;
        use nalgebra_sparse::CsrMatrix;
        use crate::sparse::spmv;

        // Create a simple 3x3 matrix:
        // [2  0  1]
        // [0  3  0]
        // [1  0  2]
        let row_offsets = vec![0, 2, 3, 5];
        let col_indices = vec![0, 2, 1, 0, 2];
        let values = vec![2.0f64, 1.0, 3.0, 1.0, 2.0];
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values).unwrap();

        // Test with x = [1, 2, 3]
        let x = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut y = DVector::zeros(3);

        spmv(&a, &x, &mut y);

        // Expected: [2*1 + 1*3, 3*2, 1*1 + 2*3] = [5, 6, 7]
        assert_relative_eq!(y[0], 5.0, epsilon = 1e-10);
        assert_relative_eq!(y[1], 6.0, epsilon = 1e-10);
        assert_relative_eq!(y[2], 7.0, epsilon = 1e-10);
    }

    #[test]
    fn test_spmv_f32_simd_correctness() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_f32_simd};

        // Create a larger matrix to test SIMD (20x20)
        let n = 20;
        let mut builder = SparseMatrixBuilder::new(n, n);
        
        // Create a tridiagonal matrix with some extra entries
        for i in 0..n {
            builder.add_triplets(vec![(i, i, 4.0f32)]).unwrap();
            if i > 0 {
                builder.add_triplets(vec![(i, i - 1, -1.0f32)]).unwrap();
            }
            if i < n - 1 {
                builder.add_triplets(vec![(i, i + 1, -1.0f32)]).unwrap();
            }
        }
        let a = builder.build().unwrap();

        // Test vector
        let x = DVector::from_fn(n, |i, _| (i + 1) as f32);
        
        // Compute with scalar version
        let mut y_scalar = DVector::zeros(n);
        spmv(&a, &x, &mut y_scalar);

        // Compute with SIMD version
        let mut y_simd = DVector::zeros(n);
        spmv_f32_simd(&a, &x, &mut y_simd);

        // Compare results
        for i in 0..n {
            assert_relative_eq!(y_simd[i], y_scalar[i], epsilon = 1e-5);
        }
    }

    #[test]
    fn test_spmv_f32_simd_sparse_rows() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_f32_simd};

        // Test with very sparse rows (< 4 non-zeros per row)
        let mut builder = SparseMatrixBuilder::new(10, 10);
        for i in 0..10 {
            builder.add_triplets(vec![(i, i, 1.0f32)]).unwrap();
            if i < 5 {
                builder.add_triplets(vec![(i, i + 5, 0.5f32)]).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let x = DVector::from_element(10, 1.0f32);
        
        let mut y_scalar = DVector::zeros(10);
        spmv(&a, &x, &mut y_scalar);

        let mut y_simd = DVector::zeros(10);
        spmv_f32_simd(&a, &x, &mut y_simd);

        for i in 0..10 {
            assert_relative_eq!(y_simd[i], y_scalar[i], epsilon = 1e-5);
        }
    }

    #[test]
    fn test_spmv_f32_simd_dense_rows() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_f32_simd};

        // Test with dense rows (many non-zeros, good for SIMD)
        let n = 16;
        let mut builder = SparseMatrixBuilder::new(n, n);
        for i in 0..n {
            for j in 0..n {
                if (i as i32 - j as i32).abs() <= 3 {
                    let val = 1.0f32 / ((i as i32 - j as i32).abs() + 1) as f32;
                    builder.add_triplets(vec![(i, j, val)]).unwrap();
                }
            }
        }
        let a = builder.build().unwrap();

        let x = DVector::from_fn(n, |i, _| (i % 3) as f32 + 0.5);
        
        let mut y_scalar = DVector::zeros(n);
        spmv(&a, &x, &mut y_scalar);

        let mut y_simd = DVector::zeros(n);
        spmv_f32_simd(&a, &x, &mut y_simd);

        for i in 0..n {
            assert_relative_eq!(y_simd[i], y_scalar[i], epsilon = 1e-4);
        }
    }

    #[test]
    fn test_spmv_parallel_correctness() {
        use nalgebra::DVector;
        use nalgebra_sparse::CsrMatrix;
        use crate::sparse::{spmv, spmv_parallel};

        // Create a simple 3x3 matrix:
        // [2  0  1]
        // [0  3  0]
        // [1  0  2]
        let row_offsets = vec![0, 2, 3, 5];
        let col_indices = vec![0, 2, 1, 0, 2];
        let values = vec![2.0f64, 1.0, 3.0, 1.0, 2.0];
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values).unwrap();

        // Test with x = [1, 2, 3]
        let x = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        // Compute with scalar version
        let mut y_scalar = DVector::zeros(3);
        spmv(&a, &x, &mut y_scalar);

        // Compute with parallel version
        let mut y_parallel = DVector::zeros(3);
        spmv_parallel(&a, &x, &mut y_parallel);

        // Expected: [2*1 + 1*3, 3*2, 1*1 + 2*3] = [5, 6, 7]
        for i in 0..3 {
            assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_spmv_parallel_large_matrix() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_parallel};

        // Create a larger matrix to benefit from parallelization (1000x1000)
        let n = 1000;
        let mut builder = SparseMatrixBuilder::new(n, n);
        
        // Create a tridiagonal matrix
        for i in 0..n {
            builder.add_triplets(vec![(i, i, 4.0f64)]).unwrap();
            if i > 0 {
                builder.add_triplets(vec![(i, i - 1, -1.0f64)]).unwrap();
            }
            if i < n - 1 {
                builder.add_triplets(vec![(i, i + 1, -1.0f64)]).unwrap();
            }
        }
        let a = builder.build().unwrap();

        // Test vector
        let x = DVector::from_fn(n, |i, _| (i + 1) as f64);
        
        // Compute with scalar version
        let mut y_scalar = DVector::zeros(n);
        spmv(&a, &x, &mut y_scalar);

        // Compute with parallel version
        let mut y_parallel = DVector::zeros(n);
        spmv_parallel(&a, &x, &mut y_parallel);

        // Compare results
        for i in 0..n {
            assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_spmv_parallel_five_point_stencil() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_parallel};

        // Create a 50x50 five-point stencil (2500x2500 matrix, ~12k non-zeros)
        let nx = 50;
        let ny = 50;
        let matrix = SparsePatterns::five_point_stencil(nx, ny, 1.0, 1.0).unwrap();

        // Test vector with varying values
        let n = nx * ny;
        let x = DVector::from_fn(n, |i, _| ((i % 10) as f64) * 0.1 + 1.0);
        
        // Compute with scalar version
        let mut y_scalar = DVector::zeros(n);
        spmv(&matrix, &x, &mut y_scalar);

        // Compute with parallel version
        let mut y_parallel = DVector::zeros(n);
        spmv_parallel(&matrix, &x, &mut y_parallel);

        // Compare results
        for i in 0..n {
            assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_spmv_parallel_sparse_pattern() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_parallel};

        // Test with very sparse rows (edge case for parallel overhead)
        let mut builder = SparseMatrixBuilder::new(100, 100);
        for i in 0..100 {
            builder.add_triplets(vec![(i, i, 1.0f64)]).unwrap();
            // Add one off-diagonal entry every 5 rows
            if i % 5 == 0 && i < 95 {
                builder.add_triplets(vec![(i, i + 5, 0.5f64)]).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let x = DVector::from_element(100, 2.0f64);
        
        let mut y_scalar = DVector::zeros(100);
        spmv(&a, &x, &mut y_scalar);

        let mut y_parallel = DVector::zeros(100);
        spmv_parallel(&a, &x, &mut y_parallel);

        for i in 0..100 {
            assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_spmv_parallel_dense_block() {
        use nalgebra::DVector;
        use crate::sparse::{spmv, spmv_parallel};

        // Test with denser structure (pentadiagonal-like pattern)
        let n = 500;
        let mut builder = SparseMatrixBuilder::new(n, n);
        for i in 0..n {
            // Main diagonal
            builder.add_triplets(vec![(i, i, 4.0f64)]).unwrap();
            // Off-diagonals (bandwidth = 2)
            if i > 0 {
                builder.add_triplets(vec![(i, i - 1, -1.0f64)]).unwrap();
            }
            if i > 1 {
                builder.add_triplets(vec![(i, i - 2, -0.5f64)]).unwrap();
            }
            if i < n - 1 {
                builder.add_triplets(vec![(i, i + 1, -1.0f64)]).unwrap();
            }
            if i < n - 2 {
                builder.add_triplets(vec![(i, i + 2, -0.5f64)]).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let x = DVector::from_fn(n, |i, _| (i as f64).sin());
        
        let mut y_scalar = DVector::zeros(n);
        spmv(&a, &x, &mut y_scalar);

        let mut y_parallel = DVector::zeros(n);
        spmv_parallel(&a, &x, &mut y_parallel);

        for i in 0..n {
            assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
        }
    }
}
