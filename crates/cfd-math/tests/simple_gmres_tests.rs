//! Basic GMRES solver tests to validate core functionality

use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{IterativeLinearSolver, LinearOperator, GMRES};
use leto::Array1;
use leto_ops::CsrMatrix;

fn array(values: Vec<f64>) -> Array1<f64> {
    Array1::from_shape_vec([values.len()], values).expect("valid Leto vector shape")
}

fn filled_array(len: usize, value: f64) -> Array1<f64> {
    Array1::from_shape_vec([len], vec![value; len]).expect("valid Leto vector shape")
}

fn csr_from_entries(nrows: usize, ncols: usize, entries: &[(usize, usize, f64)]) -> CsrMatrix<f64> {
    let mut rows = vec![Vec::<(usize, f64)>::new(); nrows];
    for &(row, col, value) in entries {
        rows[row].push((col, value));
    }

    let mut row_offsets = Vec::with_capacity(nrows + 1);
    let mut col_indices = Vec::with_capacity(entries.len());
    let mut values = Vec::with_capacity(entries.len());
    row_offsets.push(0);
    for row in &mut rows {
        row.sort_by_key(|&(col, _)| col);
        for &(col, value) in row.iter() {
            col_indices.push(col);
            values.push(value);
        }
        row_offsets.push(col_indices.len());
    }

    CsrMatrix::from_parts(values, col_indices, row_offsets, nrows, ncols)
        .expect("valid Leto CSR matrix")
}

fn residual_norm(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let mut ax = Array1::zeros([b.shape()[0]]);
    a.apply(x, &mut ax).expect("Leto CSR operator application");
    let mut sum = 0.0;
    for idx in 0..b.shape()[0] {
        let diff = ax[idx] - b[idx];
        sum += diff * diff;
    }
    sum.sqrt()
}

fn assert_residual_below(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>, threshold: f64) {
    let residual = residual_norm(a, x, b);
    assert!(
        residual < threshold,
        "GMRES residual {residual} exceeded {threshold}"
    );
}

/// Test basic GMRES functionality
#[test]
fn test_gmres_basic() {
    let n = 5;
    let a = csr_from_entries(
        n,
        n,
        &[
            (0, 0, 2.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 2.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 3.0),
            (2, 3, -1.0),
            (3, 2, -1.0),
            (3, 3, 2.0),
            (3, 4, -1.0),
            (4, 3, -1.0),
            (4, 4, 2.0),
        ],
    );

    let b = array(vec![1.0, 2.0, 3.0, 2.0, 1.0]);
    let mut x = Array1::zeros([n]);
    let config = cfd_math::linear_solver::IterativeSolverConfig::new(1e-8).with_max_iterations(100);
    let solver = GMRES::new(config, 4); // Restart after 4 iterations

    let result = solver.solve(&a, &b, &mut x, None::<&IdentityPreconditioner>);
    assert!(result.is_ok(), "GMRES should converge");

    assert_residual_below(&a, &x, &b, 1.5);
}

/// Test GMRES restart mechanism
#[test]
fn test_gmres_restart() {
    let n = 8;
    let mut entries = Vec::with_capacity(3 * n - 2);
    for i in 0..n {
        entries.push((i, i, 2.0));
        if i > 0 {
            entries.push((i, i - 1, -1.0));
        }
        if i < n - 1 {
            entries.push((i, i + 1, -1.0));
        }
    }
    let a = csr_from_entries(n, n, &entries);

    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);

    // Use small restart dimension to force restart testing
    let config = cfd_math::linear_solver::IterativeSolverConfig::new(1e-8).with_max_iterations(200);
    let solver = GMRES::new(config, 3); // Small restart dimension

    let result = solver.solve(&a, &b, &mut x, None::<&IdentityPreconditioner>);
    assert!(result.is_ok(), "GMRES with restart should converge");
    assert_residual_below(&a, &x, &b, 1e-6);
}

/// Test GMRES with preconditioner
#[test]
fn test_gmres_with_preconditioner() {
    let n = 6;
    let mut entries = Vec::with_capacity(n);
    for i in 0..n {
        entries.push((i, i, 1.0 + (i as f64) * 0.1)); // Diagonally dominant
    }
    let a = csr_from_entries(n, n, &entries);

    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);
    let config = cfd_math::linear_solver::IterativeSolverConfig::new(1e-10).with_max_iterations(50);
    let solver = GMRES::new(config, n);
    let precond = IdentityPreconditioner;

    let result = solver.solve(&a, &b, &mut x, Some(&precond));
    assert!(result.is_ok(), "GMRES with preconditioner should work");

    assert_residual_below(&a, &x, &b, 1e-8);
}
