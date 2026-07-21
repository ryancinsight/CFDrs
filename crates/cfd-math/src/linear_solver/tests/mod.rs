//! Tests for linear solvers

#[cfg(test)]
use crate::linear_solver::preconditioners::{
    IdentityPreconditioner, IncompleteLU, JacobiPreconditioner, SORPreconditioner,
};
use crate::linear_solver::traits::IterativeLinearSolver;
use crate::linear_solver::IterativeSolverConfig;
use crate::linear_solver::{BiCGSTAB, ConjugateGradient, Preconditioner, GMRES};
use crate::sparse::{self, SparseMatrix, SparseMatrixBuilder};
use eunomia::assert_relative_eq;
use cfd_core::error::Error;
use leto::Array1;
use leto_ops::CsrMatrix;

fn array(values: Vec<f64>) -> Array1<f64> {
    Array1::from_shape_vec([values.len()], values).expect("valid Leto vector shape")
}

fn array_from_fn(len: usize, mut value_at: impl FnMut(usize) -> f64) -> Array1<f64> {
    array((0..len).map(&mut value_at).collect())
}

fn filled_array(len: usize, value: f64) -> Array1<f64> {
    Array1::from_shape_vec([len], vec![value; len]).expect("valid Leto vector shape")
}

fn fill_array(values: &mut Array1<f64>, value: f64) {
    for idx in 0..values.shape()[0] {
        values[idx] = value;
    }
}

fn matrix_from_entries(
    n: usize,
    entries: impl IntoIterator<Item = (usize, usize, f64)>,
) -> CsrMatrix<f64> {
    let mut builder = SparseMatrixBuilder::new(n, n);
    builder
        .add_triplets(entries)
        .expect("invariant: test matrix entries are in bounds");
    builder
        .build()
        .expect("invariant: test CSR matrix is valid")
}

fn diagonal_matrix(diagonal: &[f64]) -> CsrMatrix<f64> {
    matrix_from_entries(
        diagonal.len(),
        diagonal
            .iter()
            .enumerate()
            .map(|(idx, &value)| (idx, idx, value)),
    )
}

fn vector_norm(values: &Array1<f64>) -> f64 {
    (0..values.shape()[0])
        .map(|idx| values[idx] * values[idx])
        .sum::<f64>()
        .sqrt()
}

fn vector_distance(lhs: &Array1<f64>, rhs: &Array1<f64>) -> f64 {
    assert_eq!(lhs.shape(), rhs.shape(), "vector shapes must match");
    (0..lhs.shape()[0])
        .map(|idx| {
            let diff = lhs[idx] - rhs[idx];
            diff * diff
        })
        .sum::<f64>()
        .sqrt()
}

fn residual_norm(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let mut ax = Array1::zeros([b.shape()[0]]);
    sparse::spmv(a, x, &mut ax);
    vector_distance(&ax, b)
}

fn relative_residual_norm(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>) -> f64 {
    residual_norm(a, x, b) / vector_norm(b)
}

fn assert_spmv_matches(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>, epsilon: f64) {
    let mut ax = Array1::zeros([b.shape()[0]]);
    sparse::spmv(a, x, &mut ax);
    for idx in 0..b.shape()[0] {
        assert_relative_eq!(ax[idx], b[idx], epsilon = epsilon);
    }
}

fn assert_array_close(actual: &Array1<f64>, expected: &Array1<f64>, epsilon: f64) {
    assert_eq!(actual.shape(), expected.shape(), "vector shapes must match");
    for idx in 0..actual.shape()[0] {
        assert_relative_eq!(actual[idx], expected[idx], epsilon = epsilon);
    }
}

fn all_finite(values: &Array1<f64>) -> bool {
    (0..values.shape()[0]).all(|idx| values[idx].is_finite())
}

fn all_not_nan(values: &Array1<f64>) -> bool {
    (0..values.shape()[0]).all(|idx| !values[idx].is_nan())
}

fn create_tridiagonal_matrix(
    n: usize,
) -> std::result::Result<SparseMatrix<f64>, Box<dyn std::error::Error>> {
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        builder.add_entry(i, i, 2.0)?;

        if i > 0 {
            builder.add_entry(i, i - 1, -1.0)?;
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, -1.0)?;
        }
    }

    Ok(builder.build()?)
}

fn basic_preconditioner_matrix(matrix: &CsrMatrix<f64>) -> CsrMatrix<f64> {
    matrix.clone()
}

fn assert_invalid_configuration(result: cfd_core::error::Result<()>, expected_message: &str) {
    match result {
        Err(Error::InvalidConfiguration(message)) => assert_eq!(message, expected_message),
        Err(error) => panic!("expected invalid configuration, got {error:?}"),
        Ok(()) => panic!("expected invalid configuration error"),
    }
}

#[test]
fn test_conjugate_gradient() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let b = filled_array(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

    let solver = ConjugateGradient::new(config);
    let mut x = Array1::zeros([n]);
    let identity_precond = IdentityPreconditioner;
    solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

    assert_spmv_matches(&a, &x, &b, 1e-8);
    Ok(())
}

#[test]
fn test_bicgstab() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let b = filled_array(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

    let solver = BiCGSTAB::new(config);
    let mut x = Array1::zeros([n]);
    let identity_precond = IdentityPreconditioner;
    solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

    assert_spmv_matches(&a, &x, &b, 1e-8);
    Ok(())
}

#[test]
fn test_jacobi_preconditioner() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let preconditioner_matrix = basic_preconditioner_matrix(&a);
    let precond = JacobiPreconditioner::new(&preconditioner_matrix)?;

    let r = filled_array(n, 1.0);
    let mut z = Array1::zeros([n]);
    precond.apply_to(&r, &mut z)?;

    for i in 0..n {
        assert_relative_eq!(z[i], r[i] / 2.0, epsilon = 1e-10);
    }
    Ok(())
}

#[test]
fn test_sor_preconditioner() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let preconditioner_matrix = basic_preconditioner_matrix(&a);

    let sor = SORPreconditioner::with_omega_for_1d_poisson(&preconditioner_matrix)?;

    assert!(
        sor.omega() > 1.0 && sor.omega() < 2.0,
        "Omega {} should be in (1,2)",
        sor.omega()
    );

    let r = filled_array(n, 1.0);
    let mut z = Array1::zeros([n]);
    sor.apply_to(&r, &mut z)?;

    for i in 0..n {
        assert!(z[i] > 0.0, "SOR result component {i} should be positive");
    }

    for i in 0..n {
        assert!(
            z[i] <= 1.0,
            "SOR result magnitude should be bounded by RHS magnitude"
        );
        assert!(
            z[i] >= 0.1,
            "SOR result should maintain reasonable lower bound"
        );
    }

    let uniform_scaling = r[0] / 2.0;
    let mut has_variation = false;
    for i in 0..n {
        if (z[i] - uniform_scaling).abs() > 1e-6 {
            has_variation = true;
            break;
        }
    }
    assert!(
        has_variation,
        "SOR should produce non-uniform result, not simple scaling"
    );
    Ok(())
}

#[test]
fn basic_preconditioners_reject_mismatched_leto_vector_lengths(
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let preconditioner_matrix = basic_preconditioner_matrix(&a);
    let identity = IdentityPreconditioner;
    let jacobi = JacobiPreconditioner::new(&preconditioner_matrix)?;
    let sor = SORPreconditioner::with_omega_for_1d_poisson(&preconditioner_matrix)?;

    let r_short = leto::Array1::zeros([n - 1]);
    let mut z = leto::Array1::zeros([n]);
    assert_invalid_configuration(
        jacobi.apply_to(&r_short, &mut z),
        "Jacobi residual length mismatch: expected 5, got 4",
    );
    assert_invalid_configuration(
        sor.apply_to(&r_short, &mut z),
        "SOR residual length mismatch: expected 5, got 4",
    );

    let r = leto::Array1::zeros([n]);
    let mut z_short = leto::Array1::zeros([n - 1]);
    assert_invalid_configuration(
        identity.apply_to(&r, &mut z_short),
        "identity preconditioner output length mismatch: expected 5, got 4",
    );
    assert_invalid_configuration(
        jacobi.apply_to(&r, &mut z_short),
        "Jacobi output length mismatch: expected 5, got 4",
    );
    assert_invalid_configuration(
        sor.apply_to(&r, &mut z_short),
        "SOR output length mismatch: expected 5, got 4",
    );

    Ok(())
}

#[test]
fn test_ilu_works_on_non_tridiagonal() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 3;
    let mut builder = SparseMatrixBuilder::new(n, n);

    builder.add_entry(0, 0, 2.0)?;
    builder.add_entry(0, 2, 1.0)?;
    builder.add_entry(1, 1, 2.0)?;
    builder.add_entry(2, 2, 2.0)?;
    let non_tridiag = builder.build()?;

    let preconditioner_matrix = basic_preconditioner_matrix(&non_tridiag);
    let result = IncompleteLU::new(&preconditioner_matrix);
    assert!(result.is_ok());
    Ok(())
}

#[test]
fn test_preconditioned_cg() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let b = filled_array(n, 1.0);
    let preconditioner_matrix = basic_preconditioner_matrix(&a);
    let precond = JacobiPreconditioner::new(&preconditioner_matrix)?;

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

    let solver = ConjugateGradient::new(config);
    let mut x = Array1::zeros([n]);
    solver.solve(&a, &b, &mut x, Some(&precond))?;

    assert_spmv_matches(&a, &x, &b, 1e-8);
    Ok(())
}

#[test]
fn test_gauss_seidel() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let preconditioner_matrix = basic_preconditioner_matrix(&a);
    let precond = SORPreconditioner::new(&preconditioner_matrix, 1.0)?;

    let r = filled_array(n, 1.0);
    let mut z = Array1::zeros([n]);
    precond.apply_to(&r, &mut z)?;

    assert_relative_eq!(z[0], 0.5, epsilon = 1e-10);

    for i in 0..n {
        assert!(
            z[i] > 0.0,
            "Gauss-Seidel result component {i} should be positive"
        );
    }

    let mut has_variation = false;
    for i in 1..n {
        if (z[i] - z[0]).abs() > 1e-6 {
            has_variation = true;
            break;
        }
    }
    assert!(
        has_variation,
        "Gauss-Seidel should produce non-uniform solution"
    );

    for i in 0..n {
        assert!(
            z[i] <= 1.0,
            "Gauss-Seidel components should be reasonably bounded"
        );
    }
    Ok(())
}

#[test]
fn test_convergence_with_different_tolerances(
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 10;
    let a = create_tridiagonal_matrix(n)?;
    let b = filled_array(n, 1.0);

    let tolerances = vec![1e-4, 1e-6, 1e-8];

    for tol in tolerances {
        let config = IterativeSolverConfig::new(tol).with_max_iterations(1000);

        let solver = ConjugateGradient::new(config);
        let mut x = Array1::zeros([n]);
        let identity_precond = IdentityPreconditioner;
        solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

        let relative_residual = relative_residual_norm(&a, &x, &b);
        assert!(relative_residual < tol * 10.0);
    }
    Ok(())
}

#[test]
fn test_gmres_diagonal_matrix() {
    let n = 4;
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);

    let b = array(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&precond)).unwrap();

    let expected = array(vec![1.0, 2.0, 3.0, 4.0]);
    let error = vector_distance(&x, &expected);
    assert!(error < 1e-9, "Solution error: {error}");
}

#[test]
fn test_gmres_non_symmetric() {
    let n = 3;
    let a = matrix_from_entries(
        n,
        [
            (0, 0, 4.0),
            (0, 1, 1.0),
            (1, 0, 2.0),
            (1, 1, 5.0),
            (1, 2, 1.0),
            (2, 1, 3.0),
            (2, 2, 6.0),
        ],
    );

    let b = array(vec![5.0, 8.0, 9.0]);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&precond)).unwrap();

    let residual = residual_norm(&a, &x, &b);
    assert!(residual < 1e-9, "Residual: {residual}");
}

#[test]
fn test_gmres_zero_rhs() {
    let n = 4;
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);

    let b = Array1::zeros([n]);
    let mut x = filled_array(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&precond)).unwrap();

    assert!(
        vector_norm(&x) < 1e-14,
        "Solution should be zero for zero RHS"
    );
}

#[test]
fn test_gmres_restart_dimension() {
    let n = 5;
    let a = matrix_from_entries(
        n,
        (0..n).flat_map(|i| {
            let mut entries = vec![(i, i, 4.0)];
            if i > 0 {
                entries.push((i, i - 1, 1.0));
            }
            if i < n - 1 {
                entries.push((i, i + 1, 1.0));
            }
            entries
        }),
    );

    let b = array(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 3);
    let precond = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&precond)).unwrap();

    let residual = residual_norm(&a, &x, &b);
    assert!(residual < 1e-9, "Residual: {residual}");
}

#[test]
fn test_gmres_with_initial_guess() {
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);

    let b = array(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = array(vec![0.5, 1.5, 2.5, 3.5]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&precond)).unwrap();

    let expected = array(vec![1.0, 2.0, 3.0, 4.0]);
    let error = vector_distance(&x, &expected);
    assert!(error < 1e-9, "Solution error: {error}");
}

#[test]
fn test_gmres_larger_nonsymmetric() {
    let n = 5;
    let a = matrix_from_entries(
        n,
        (0..n).flat_map(|i| {
            let mut entries = vec![(i, i, 5.0)];
            if i > 0 {
                entries.push((i, i - 1, 2.0));
            }
            if i < n - 1 {
                entries.push((i, i + 1, 1.0));
            }
            entries
        }),
    );

    let b = array(vec![6.0, 11.0, 11.0, 11.0, 8.0]);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&precond)).unwrap();

    let residual = residual_norm(&a, &x, &b);
    assert!(residual < 1e-9, "Residual: {residual}");
}

#[test]
fn test_gmres_convergence_tight_tolerance() {
    let n = 4;
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);

    let b = array(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-14).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    let result = solver.solve(&a, &b, &mut x, Some(&precond));
    assert!(result.is_ok());
}

#[test]
fn test_gmres_max_iterations_exceeded() {
    let n = 4;
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);

    let b = array(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-14).with_max_iterations(1);
    let solver = GMRES::new(config, 1);
    let precond = IdentityPreconditioner;

    let result = solver.solve(&a, &b, &mut x, Some(&precond));
    assert!(result.is_err());
}

#[test]
fn test_gmres_dimension_mismatch() {
    let a = diagonal_matrix(&[1.0, 2.0, 3.0, 4.0]);

    let b = array(vec![1.0, 4.0]);
    let mut x = Array1::zeros([2]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    let result = solver.solve(&a, &b, &mut x, Some(&precond));
    assert!(result.is_err());
}

#[test]
fn test_gmres_configurable_trait() {
    use crate::linear_solver::traits::Configurable;

    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(200);
    let solver = GMRES::<f64>::new(config, 20);

    let retrieved = solver.config();
    assert!((retrieved.tolerance - 1e-8).abs() < 1e-10);
    assert_eq!(retrieved.max_iterations, 200);
}

#[test]
#[should_panic(expected = "GMRES restart dimension must be positive")]
fn test_gmres_zero_restart_dim_panics() {
    let config = IterativeSolverConfig::default();
    let _solver = GMRES::<f64>::new(config, 0);
}

mod adversarial_solver_tests;
mod edge_case_tests;
mod extended_edge_case_tests;
