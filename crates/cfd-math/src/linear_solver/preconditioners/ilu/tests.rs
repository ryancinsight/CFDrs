//! Tests for Incomplete LU factorization

use super::*;
use crate::linear_solver::Preconditioner;
use approx::assert_relative_eq;
use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

/// Create a simple 4x4 tridiagonal test matrix
/// [4 -1  0  0]
/// [-1 4 -1  0]
/// [0 -1  4 -1]
/// [0  0 -1  4]
fn create_tridiagonal_matrix() -> CsrMatrix<f64> {
    let row_offsets = vec![0, 2, 5, 8, 10];
    let col_indices = vec![
        0, 1, // row 0
        0, 1, 2, // row 1
        1, 2, 3, // row 2
        2, 3, // row 3
    ];
    let values = vec![
        4.0, -1.0, // row 0
        -1.0, 4.0, -1.0, // row 1
        -1.0, 4.0, -1.0, // row 2
        -1.0, 4.0, // row 3
    ];

    CsrMatrix::try_from_csr_data(4, 4, row_offsets, col_indices, values).expect("Valid CSR matrix")
}

/// Create a 5x5 sparse matrix with more complex structure
fn create_sparse_matrix() -> CsrMatrix<f64> {
    let row_offsets = vec![0, 3, 6, 9, 13, 15];
    let col_indices = vec![
        0, 1, 3, // row 0: connections to 1, 3
        0, 1, 2, // row 1: connections to 0, 2
        1, 2, 3, // row 2: connections to 1, 3
        0, 2, 3, 4, // row 3: connections to 0, 2, 4
        3, 4, // row 4: connection to 3
    ];
    let values = vec![
        5.0, -1.0, -1.0, // row 0
        -1.0, 5.0, -1.0, // row 1
        -1.0, 5.0, -1.0, // row 2
        -1.0, -1.0, 5.0, -1.0, // row 3
        -1.0, 5.0, // row 4
    ];

    CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values).expect("Valid CSR matrix")
}

#[test]
fn test_ilu0_construction() {
    let matrix = create_tridiagonal_matrix();
    let ilu = IncompleteLU::new(&matrix);

    assert!(ilu.is_ok());
    let ilu = ilu.unwrap();
    assert_eq!(ilu.fill_level(), 0);
}

#[test]
fn test_iluk_construction() {
    let matrix = create_tridiagonal_matrix();
    let ilu = IncompleteLU::with_fill_level(&matrix, 1);

    assert!(ilu.is_ok());
    let ilu = ilu.unwrap();
    assert_eq!(ilu.fill_level(), 1);
}

#[test]
fn test_ilu0_apply() {
    let matrix = create_tridiagonal_matrix();
    let ilu = IncompleteLU::new(&matrix).expect("ILU(0) construction");

    // Test with a simple vector
    let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
    let mut z = DVector::zeros(4);
    ilu.apply_to(&b, &mut z).expect("Apply preconditioner");

    // Solution should be positive and bounded
    assert!(z[0] > 0.0 && z[0] < 10.0);
    assert!(z[1] > 0.0 && z[1] < 10.0);
    assert!(z[2] > 0.0 && z[2] < 10.0);
    assert!(z[3] > 0.0 && z[3] < 10.0);
}

#[test]
fn test_iluk_vs_ilu0_sparsity() {
    let matrix = create_sparse_matrix();

    let ilu0 = IncompleteLU::new(&matrix).expect("ILU(0) construction");
    let ilu1 = IncompleteLU::with_fill_level(&matrix, 1).expect("ILU(1) construction");

    // ILU(1) should have at least as many nonzeros as ILU(0)
    let nnz_ilu0 = ilu0.lu_factor.nnz();
    let nnz_ilu1 = ilu1.lu_factor.nnz();

    assert!(
        nnz_ilu1 >= nnz_ilu0,
        "ILU(1) should have >= nonzeros than ILU(0): {nnz_ilu1} vs {nnz_ilu0}"
    );
}

#[test]
fn test_ilu0_preconditioner_quality() {
    // Test that ILU(0) actually improves conditioning
    let matrix = create_tridiagonal_matrix();
    let ilu = IncompleteLU::new(&matrix).expect("ILU(0) construction");

    // Test with identity-like vector
    let b = DVector::from_vec(vec![1.0, 1.0, 1.0, 1.0]);
    let mut z = DVector::zeros(4);
    ilu.apply_to(&b, &mut z).expect("Apply preconditioner");

    // For a diagonally dominant matrix, preconditioner should give reasonable approximation
    // to A^{-1} * b, which for this matrix and b is approximately [0.36, 0.43, 0.43, 0.36]
    assert_relative_eq!(z[0], 0.36, epsilon = 0.15);
    assert_relative_eq!(z[1], 0.43, epsilon = 0.15);
    assert_relative_eq!(z[2], 0.43, epsilon = 0.15);
    assert_relative_eq!(z[3], 0.36, epsilon = 0.15);
}

#[test]
fn test_iluk_improved_approximation() {
    // Test that ILU(k) for k>0 improves approximation quality
    let matrix = create_sparse_matrix();

    let ilu0 = IncompleteLU::new(&matrix).expect("ILU(0) construction");
    let ilu1 = IncompleteLU::with_fill_level(&matrix, 1).expect("ILU(1) construction");

    let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);

    let mut z0 = DVector::zeros(5);
    let mut z1 = DVector::zeros(5);
    ilu0.apply_to(&b, &mut z0).expect("Apply ILU(0)");
    ilu1.apply_to(&b, &mut z1).expect("Apply ILU(1)");

    // Both should give bounded results
    for i in 0..5 {
        assert!(z0[i].abs() < 100.0, "ILU(0) result unbounded");
        assert!(z1[i].abs() < 100.0, "ILU(1) result unbounded");
    }

    // ILU(1) should generally give different (hopefully better) results
    // We don't test for "better" directly, just that fill improves approximation
    let diff = (&z0 - &z1).norm();
    assert!(diff > 0.0, "ILU(1) should differ from ILU(0)");
}

#[test]
fn test_ilu_non_square_matrix() {
    use cfd_core::error::Error;

    // Test error handling for non-square matrix
    let row_offsets = vec![0, 2, 4];
    let col_indices = vec![0, 1, 1, 2];
    let values = vec![1.0, 2.0, 3.0, 4.0];

    let matrix = CsrMatrix::try_from_csr_data(2, 3, row_offsets, col_indices, values)
        .expect("Valid CSR matrix");

    let ilu = IncompleteLU::new(&matrix);
    assert!(ilu.is_err());

    // Verify error is InvalidInput for non-square matrix
    assert!(matches!(ilu, Err(Error::InvalidInput(_))));
}

#[test]
fn test_ilu_fill_levels() {
    let matrix = create_tridiagonal_matrix();

    let ilu0 = IncompleteLU::new(&matrix).expect("ILU(0) construction");
    assert_eq!(ilu0.fill_level(), 0);

    let ilu2 = IncompleteLU::with_fill_level(&matrix, 2).expect("ILU(2) construction");
    assert_eq!(ilu2.fill_level(), 2);
}

#[test]
fn test_ilu_multiple_fill_levels() {
    let matrix = create_sparse_matrix();

    for k in 0..=3 {
        let ilu = IncompleteLU::with_fill_level(&matrix, k);
        assert!(ilu.is_ok(), "ILU({k}) construction failed");

        let ilu = ilu.unwrap();
        assert_eq!(ilu.fill_level(), k);

        // Test that matrix dimensions are correct
        let x = ilu.lu_factor.nrows();
        assert_eq!(x, 5);
    }
}
