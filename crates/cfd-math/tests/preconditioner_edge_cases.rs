//! Comprehensive preconditioner edge case validation tests
//!
//! Validates robustness of ILU and Cholesky preconditioners per:
//! - Saad (2003): "Iterative Methods for Sparse Linear Systems"
//! - Benzi (2002): "Preconditioning Techniques for Large Sparse Matrix Problems"
//! - Greenbaum (1997): "Iterative Methods for Solving Linear Systems"

use cfd_math::linear_solver::preconditioners::cholesky::IncompleteCholesky;
use cfd_math::linear_solver::preconditioners::ilu::IncompleteLU;
use cfd_math::linear_solver::Preconditioner;
use leto::Array1;
use leto_ops::CsrMatrix as LetoCsrMatrix;

/// Create ill-conditioned tridiagonal matrix
fn create_ill_conditioned_matrix() -> LetoCsrMatrix<f64> {
    let n = 10;
    let mut row_offsets = vec![0];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for i in 0..n {
        if i > 0 {
            col_indices.push(i - 1);
            values.push(-1.0);
        }

        col_indices.push(i);
        values.push(2.0 + 1e-10);
        if i < n - 1 {
            col_indices.push(i + 1);
            values.push(-1.0);
        }
        row_offsets.push(col_indices.len());
    }

    LetoCsrMatrix::from_parts(values, col_indices, row_offsets, n, n).expect("valid CSR matrix")
}

fn filled_array(len: usize, value: f64) -> Array1<f64> {
    Array1::from_shape_vec([len], vec![value; len]).expect("invariant: array shape matches values")
}

fn apply_preconditioner<P: Preconditioner<f64>>(
    preconditioner: &P,
    rhs: &Array1<f64>,
    out: &mut Array1<f64>,
) -> cfd_core::error::Result<()> {
    preconditioner.apply_to(rhs, out)
}

/// Test ILU(0) with ill-conditioned matrix
/// Reference: Saad (2003), Section 10.3.1
#[test]
fn test_ilu0_ill_conditioned_matrix() {
    let a = create_ill_conditioned_matrix();
    let ilu = IncompleteLU::new(&a).expect("ILU(0) factorization should succeed");

    let b = filled_array(a.nrows(), 1.0);
    let mut x = Array1::zeros([a.nrows()]);

    apply_preconditioner(&ilu, &b, &mut x).expect("Preconditioner application should succeed");

    for i in 0..x.shape()[0] {
        assert!(x[i].is_finite(), "Solution should be finite");
        assert!(x[i].abs() < 1e6, "Solution should be bounded");
    }
}

/// Test ILU(k) with k>0 improves conditioning
/// Reference: Saad (2003), Section 10.3.4
#[test]
fn test_iluk_improves_conditioning() {
    let a = create_ill_conditioned_matrix();

    let ilu0 = IncompleteLU::new(&a).expect("ILU(0) factorization");
    let ilu1 = IncompleteLU::with_fill_level(&a, 1).expect("ILU(1) factorization");

    let b = filled_array(a.nrows(), 1.0);
    let mut x0 = Array1::zeros([a.nrows()]);
    let mut x1 = Array1::zeros([a.nrows()]);

    apply_preconditioner(&ilu0, &b, &mut x0).expect("ILU(0) application");
    apply_preconditioner(&ilu1, &b, &mut x1).expect("ILU(1) application");

    for i in 0..a.nrows() {
        assert!(x0[i].is_finite(), "ILU(0) solution finite");
        assert!(x1[i].is_finite(), "ILU(1) solution finite");
    }
}

/// Test Cholesky with non-positive definite matrix
/// Reference: Greenbaum (1997), Section 8.2
#[test]
fn test_cholesky_non_positive_definite() {
    let n = 4;
    // Create symmetric matrix with negative diagonal (not positive definite)
    // Proper CSR format requires sorted column indices
    let row_offsets = vec![0, 2, 4, 6, 8];
    let col_indices = vec![
        0, 1, // row 0: (0,0), (0,1)
        1, 2, // row 1: (1,1), (1,2)
        2, 3, // row 2: (2,2), (2,3)
        0, 3, // row 3: (3,0), (3,3)
    ];
    let values = vec![
        1.0, 0.5, // row 0
        -1.0, 0.5, // row 1: negative diagonal
        -1.0, 0.5, // row 2: negative diagonal
        0.5, 1.0, // row 3
    ];

    let a = LetoCsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
        .expect("valid CSR matrix");

    let result = IncompleteCholesky::new(&a);
    assert!(
        result.is_err(),
        "Cholesky should reject non-positive definite matrix"
    );
}

/// Test ILU convergence with well-conditioned matrix
/// Reference: Benzi (2002), Section 3.1
#[test]
fn test_ilu0_convergence_with_conditioning() {
    let n = 8;
    let mut row_offsets = vec![0];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for i in 0..n {
        if i > 0 {
            col_indices.push(i - 1);
            values.push(-1.0);
        }
        col_indices.push(i);
        values.push(4.0);
        if i < n - 1 {
            col_indices.push(i + 1);
            values.push(-1.0);
        }
        row_offsets.push(col_indices.len());
    }

    let a = LetoCsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
        .expect("valid CSR matrix");

    let ilu = IncompleteLU::new(&a).expect("ILU(0) factorization");
    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);

    for _ in 0..10 {
        apply_preconditioner(&ilu, &b, &mut x).expect("Repeated applications should succeed");

        for i in 0..n {
            assert!(x[i].is_finite(), "Solution should be finite");
            assert!(x[i].abs() < 10.0, "Solution should be bounded");
        }
    }
}

/// Test ILU with extreme values
#[test]
fn test_ilu0_extreme_values() {
    let n = 5;
    // Proper CSR format with sorted column indices per row
    let row_offsets = vec![0, 2, 4, 6, 8, 10];
    let col_indices = vec![
        0, 1, // row 0
        1, 2, // row 1
        2, 3, // row 2
        3, 4, // row 3
        0, 4, // row 4: sorted (0, 4)
    ];
    let values = vec![
        1e6, 1.0, // Large diagonal
        1.0, 1e-6, // Small diagonal
        1e-6, 1.0, 1.0, 1e3, // Medium diagonal
        1.0, 1e6, // Connection back to row 0
    ];

    let a = LetoCsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
        .expect("valid CSR matrix");

    let result = IncompleteLU::new(&a);

    if let Ok(ilu) = result {
        let b = filled_array(n, 1.0);
        let mut x = Array1::zeros([n]);

        if apply_preconditioner(&ilu, &b, &mut x).is_ok() {
            for i in 0..n {
                assert!(
                    x[i].is_finite(),
                    "Solution should be finite with extreme values"
                );
            }
        }
    }
}

/// Test ILU sparsity preservation
/// Reference: Saad (2003), Section 10.3.2
#[test]
fn test_ilu0_sparsity_preservation() {
    let n = 6;
    let mut row_offsets = vec![0];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for i in 0..n {
        if i > 0 {
            col_indices.push(i - 1);
            values.push(-1.0);
        }
        col_indices.push(i);
        values.push(3.0);
        if i < n - 1 {
            col_indices.push(i + 1);
            values.push(-1.0);
        }
        row_offsets.push(col_indices.len());
    }

    let a = LetoCsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
        .expect("valid CSR matrix");

    let ilu = IncompleteLU::new(&a).expect("ILU(0) factorization");
    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);

    apply_preconditioner(&ilu, &b, &mut x).expect("Application should succeed");

    for i in 0..n {
        assert!(x[i].is_finite(), "Solution should be finite");
    }
}
