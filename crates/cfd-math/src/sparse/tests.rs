//! Tests for sparse matrix module

#[cfg(test)]
use super::*;
use approx::assert_relative_eq;
use cfd_core::error::Result;
use leto::Array1;

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
fn matrix_builder_dirichlet_column_elimination_updates_leto_rhs() -> Result<()> {
    let mut builder = SparseMatrixBuilder::new(3, 3);
    builder.add_triplets(vec![
        (0, 0, 4.0),
        (0, 1, 2.0),
        (1, 0, 2.0),
        (1, 1, 5.0),
        (1, 2, 1.0),
        (2, 2, 3.0),
    ])?;
    builder.set_dirichlet_row(1, 10.0, 7.0);

    let mut rhs = Array1::from_shape_vec([3], vec![1.0, 70.0, 3.0]).unwrap();
    let matrix = builder.build_with_rhs(&mut rhs)?;

    assert_eq!(rhs[0], -13.0);
    assert_eq!(rhs[1], 70.0);
    assert_eq!(rhs[2], 3.0);

    let row0 = matrix.row(0);
    let row0_entries: Vec<_> = row0
        .col_indices()
        .iter()
        .zip(row0.values().iter())
        .map(|(&col, &val)| (col, val))
        .collect();
    assert_eq!(row0_entries, vec![(0, 4.0)]);

    let row1 = matrix.row(1);
    let row1_entries: Vec<_> = row1
        .col_indices()
        .iter()
        .zip(row1.values().iter())
        .map(|(&col, &val)| (col, val))
        .collect();
    assert_eq!(row1_entries, vec![(1, 10.0)]);

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
fn test_sparse_extension_scaling_and_condition_use_leto_provider() -> Result<()> {
    use leto_ops::CsrMatrix;

    let mut matrix = CsrMatrix::from_parts(
        vec![4.0f64, -1.0, 5.0, 1.5, -2.0, 6.0],
        vec![0, 2, 1, 0, 1, 2],
        vec![0, 2, 3, 6],
        3,
        3,
    )
    .unwrap();

    assert_relative_eq!(
        SparseMatrixExt::condition_estimate(&matrix)?,
        19.0 / 12.0,
        epsilon = 1e-12
    );

    matrix.scale(0.5);
    assert_eq!(matrix.values(), &[2.0, -0.5, 2.5, 0.75, -1.0, 3.0]);

    SparseMatrixExt::scale_rows(
        &mut matrix,
        &Array1::from_shape_vec([3], vec![2.0, 3.0, 4.0]).unwrap(),
    )?;
    assert_eq!(matrix.values(), &[4.0, -1.0, 7.5, 3.0, -4.0, 12.0]);

    SparseMatrixExt::scale_columns(
        &mut matrix,
        &Array1::from_shape_vec([3], vec![10.0, 20.0, 30.0]).unwrap(),
    )?;
    assert_eq!(matrix.values(), &[40.0, -30.0, 150.0, 30.0, -80.0, 360.0]);
    assert_eq!(matrix.row_ptr(), &[0, 2, 3, 6]);
    assert_eq!(matrix.col_indices(), &[0, 2, 1, 0, 1, 2]);

    assert!(SparseMatrixExt::scale_rows(
        &mut matrix,
        &Array1::from_shape_vec([2], vec![1.0, 2.0]).unwrap(),
    )
    .is_err());
    assert!(SparseMatrixExt::scale_columns(
        &mut matrix,
        &Array1::from_shape_vec([2], vec![1.0, 2.0]).unwrap(),
    )
    .is_err());

    let rectangular =
        CsrMatrix::from_parts(vec![1.0f64, 2.0], vec![0, 2], vec![0, 1, 2], 2, 3).unwrap();
    assert!(SparseMatrixExt::condition_estimate(&rectangular).is_err());

    Ok(())
}

#[test]
fn test_spmv_basic() {
    use crate::sparse::{spmv, try_spmv};
    use leto_ops::CsrMatrix;

    // Create a simple 3x3 matrix:
    // [2  0  1]
    // [0  3  0]
    // [1  0  2]
    let row_offsets = vec![0, 2, 3, 5];
    let col_indices = vec![0, 2, 1, 0, 2];
    let values = vec![2.0f64, 1.0, 3.0, 1.0, 2.0];
    let a = CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).unwrap();

    // Test with x = [1, 2, 3]
    let x = Array1::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();
    let mut y = Array1::zeros([3]);

    spmv(&a, &x, &mut y);

    // Expected: [2*1 + 1*3, 3*2, 1*1 + 2*3] = [5, 6, 7]
    assert_relative_eq!(y[0], 5.0, epsilon = 1e-10);
    assert_relative_eq!(y[1], 6.0, epsilon = 1e-10);
    assert_relative_eq!(y[2], 7.0, epsilon = 1e-10);

    let bad_x = Array1::from_shape_vec([2], vec![1.0, 2.0]).unwrap();
    assert!(try_spmv(&a, &bad_x, &mut y).is_err());
}

#[test]
fn test_sparse_sparse_mul_uses_leto_provider() {
    use crate::sparse::try_sparse_sparse_mul;
    use leto_ops::CsrMatrix;

    // A =
    // [ 2  0 -1 ]
    // [ 0  3  0 ]
    let a =
        CsrMatrix::from_parts(vec![2.0f64, -1.0, 3.0], vec![0, 2, 1], vec![0, 2, 3], 2, 3).unwrap();

    // B =
    // [ 0  4 ]
    // [ 5  0 ]
    // [ 6  7 ]
    let b = CsrMatrix::from_parts(
        vec![4.0f64, 5.0, 6.0, 7.0],
        vec![1, 0, 0, 1],
        vec![0, 1, 2, 4],
        3,
        2,
    )
    .unwrap();

    let product = try_sparse_sparse_mul(&a, &b).unwrap();

    assert_eq!(product.nrows(), 2);
    assert_eq!(product.ncols(), 2);
    assert_eq!(product.row_ptr(), &[0, 2, 3]);
    assert_eq!(product.col_indices(), &[0, 1, 0]);
    assert_relative_eq!(product.values()[0], -6.0, epsilon = 1e-12);
    assert_relative_eq!(product.values()[1], 1.0, epsilon = 1e-12);
    assert_relative_eq!(product.values()[2], 15.0, epsilon = 1e-12);

    let mismatched = CsrMatrix::from_parts(vec![], vec![], vec![0, 0, 0, 0, 0], 4, 1).unwrap();
    assert!(try_sparse_sparse_mul(&a, &mismatched).is_err());
}

#[test]
fn test_sparse_transpose_uses_leto_provider() {
    use crate::sparse::try_sparse_transpose;
    use leto_ops::CsrMatrix;

    let matrix = CsrMatrix::from_parts(
        vec![1.0f64, 2.0, 3.0, 4.0, 5.0],
        vec![0, 3, 1, 0, 2],
        vec![0, 2, 3, 5],
        3,
        4,
    )
    .unwrap();

    let transposed = try_sparse_transpose(&matrix).unwrap();

    assert_eq!(transposed.nrows(), 4);
    assert_eq!(transposed.ncols(), 3);
    assert_eq!(transposed.row_ptr(), &[0, 2, 3, 4, 5]);
    assert_eq!(transposed.col_indices(), &[0, 2, 1, 2, 0]);
    assert_relative_eq!(transposed.values()[0], 1.0, epsilon = 1e-12);
    assert_relative_eq!(transposed.values()[1], 4.0, epsilon = 1e-12);
    assert_relative_eq!(transposed.values()[2], 3.0, epsilon = 1e-12);
    assert_relative_eq!(transposed.values()[3], 5.0, epsilon = 1e-12);
    assert_relative_eq!(transposed.values()[4], 2.0, epsilon = 1e-12);

    let round_trip = try_sparse_transpose(&transposed).unwrap();
    assert_eq!(round_trip.row_ptr(), matrix.row_ptr());
    assert_eq!(round_trip.col_indices(), matrix.col_indices());
    assert_eq!(round_trip.values(), matrix.values());
}

#[test]
fn test_spmv_parallel_correctness() {
    use crate::sparse::{spmv, spmv_parallel};
    use leto_ops::CsrMatrix;

    // Create a simple 3x3 matrix:
    // [2  0  1]
    // [0  3  0]
    // [1  0  2]
    let row_offsets = vec![0, 2, 3, 5];
    let col_indices = vec![0, 2, 1, 0, 2];
    let values = vec![2.0f64, 1.0, 3.0, 1.0, 2.0];
    let a = CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).unwrap();

    // Test with x = [1, 2, 3]
    let x = Array1::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();

    // Compute with scalar version
    let mut y_scalar = Array1::zeros([3]);
    spmv(&a, &x, &mut y_scalar);

    // Compute with parallel version
    let mut y_parallel = Array1::zeros([3]);
    spmv_parallel(&a, &x, &mut y_parallel);

    // Expected: [2*1 + 1*3, 3*2, 1*1 + 2*3] = [5, 6, 7]
    for i in 0..3 {
        assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
    }
}

#[test]
fn test_spmv_parallel_large_matrix() {
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
    let x = Array1::from_shape_vec([n], (0..n).map(|i| (i + 1) as f64).collect()).unwrap();

    // Compute with scalar version
    let mut y_scalar = Array1::zeros([n]);
    spmv(&a, &x, &mut y_scalar);

    // Compute with parallel version
    let mut y_parallel = Array1::zeros([n]);
    spmv_parallel(&a, &x, &mut y_parallel);

    // Compare results
    for i in 0..n {
        assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
    }
}

#[test]
fn test_spmv_parallel_five_point_stencil() {
    use crate::sparse::{spmv, spmv_parallel};

    // Create a 50x50 five-point stencil (2500x2500 matrix, ~12k non-zeros)
    let nx = 50;
    let ny = 50;
    let matrix = SparsePatterns::five_point_stencil(nx, ny, 1.0, 1.0).unwrap();

    // Test vector with varying values
    let n = nx * ny;
    let x = Array1::from_shape_vec([n], (0..n).map(|i| ((i % 10) as f64) * 0.1 + 1.0).collect())
        .unwrap();

    // Compute with scalar version
    let mut y_scalar = Array1::zeros([n]);
    spmv(&matrix, &x, &mut y_scalar);

    // Compute with parallel version
    let mut y_parallel = Array1::zeros([n]);
    spmv_parallel(&matrix, &x, &mut y_parallel);

    // Compare results
    for i in 0..n {
        assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
    }
}

#[test]
fn test_spmv_parallel_sparse_pattern() {
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

    let x = Array1::from_shape_vec([100], vec![2.0f64; 100]).unwrap();

    let mut y_scalar = Array1::zeros([100]);
    spmv(&a, &x, &mut y_scalar);

    let mut y_parallel = Array1::zeros([100]);
    spmv_parallel(&a, &x, &mut y_parallel);

    for i in 0..100 {
        assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
    }
}

#[test]
fn test_spmv_parallel_dense_block() {
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

    let x = Array1::from_shape_vec([n], (0..n).map(|i| (i as f64).sin()).collect()).unwrap();

    let mut y_scalar = Array1::zeros([n]);
    spmv(&a, &x, &mut y_scalar);

    let mut y_parallel = Array1::zeros([n]);
    spmv_parallel(&a, &x, &mut y_parallel);

    for i in 0..n {
        assert_relative_eq!(y_parallel[i], y_scalar[i], epsilon = 1e-10);
    }
}
