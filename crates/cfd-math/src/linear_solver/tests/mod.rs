//! Tests for linear solvers

#[cfg(test)]
use crate::linear_solver::preconditioners::{
    IdentityPreconditioner, IncompleteLU, JacobiPreconditioner, SORPreconditioner,
};
use crate::linear_solver::traits::IterativeLinearSolver;
use crate::linear_solver::IterativeSolverConfig;
use crate::linear_solver::{BiCGSTAB, ConjugateGradient, Preconditioner, GMRES};
use crate::sparse::{SparseMatrix, SparseMatrixBuilder};
use approx::assert_relative_eq;
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};

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

#[test]
fn test_conjugate_gradient() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let b = DVector::from_element(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

    let solver = ConjugateGradient::new(config);
    let mut x = DVector::zeros(n);
    let identity_precond = IdentityPreconditioner;
    solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

    let ax = &a * &x;
    for i in 0..n {
        assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
    }
    Ok(())
}

#[test]
fn test_bicgstab() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let b = DVector::from_element(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

    let solver = BiCGSTAB::new(config);
    let mut x = DVector::zeros(n);
    let identity_precond = IdentityPreconditioner;
    solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

    let ax = &a * &x;
    for i in 0..n {
        assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
    }
    Ok(())
}

#[test]
fn test_jacobi_preconditioner() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let precond = JacobiPreconditioner::new(&a)?;

    let r = DVector::from_element(n, 1.0);
    let mut z = DVector::zeros(n);
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

    let sor = SORPreconditioner::with_omega_for_1d_poisson(&a)?;

    assert!(
        sor.omega() > 1.0 && sor.omega() < 2.0,
        "Omega {} should be in (1,2)",
        sor.omega()
    );

    let r = DVector::from_element(n, 1.0);
    let mut z = DVector::zeros(n);
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
fn test_ilu_works_on_non_tridiagonal() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 3;
    let mut builder = SparseMatrixBuilder::new(n, n);

    builder.add_entry(0, 0, 2.0)?;
    builder.add_entry(0, 2, 1.0)?;
    builder.add_entry(1, 1, 2.0)?;
    builder.add_entry(2, 2, 2.0)?;
    let non_tridiag = builder.build()?;

    let result = IncompleteLU::new(&non_tridiag);
    assert!(result.is_ok());
    Ok(())
}

#[test]
fn test_preconditioned_cg() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let b = DVector::from_element(n, 1.0);
    let precond = JacobiPreconditioner::new(&a)?;

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

    let solver = ConjugateGradient::new(config);
    let mut x = DVector::zeros(n);
    solver.solve(&a, &b, &mut x, Some(&precond))?;

    let ax = &a * &x;
    for i in 0..n {
        assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
    }
    Ok(())
}

#[test]
fn test_gauss_seidel() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let a = create_tridiagonal_matrix(n)?;
    let precond = SORPreconditioner::new(&a, 1.0)?;

    let r = DVector::from_element(n, 1.0);
    let mut z = DVector::zeros(n);
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
    let b = DVector::from_element(n, 1.0);

    let tolerances = vec![1e-4, 1e-6, 1e-8];

    for tol in tolerances {
        let config = IterativeSolverConfig::new(tol).with_max_iterations(1000);

        let solver = ConjugateGradient::new(config);
        let mut x = DVector::zeros(n);
        let identity_precond = IdentityPreconditioner;
        solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

        let ax = &a * &x;
        let residual = &b - &ax;
        let relative_residual = residual.norm() / b.norm();
        assert!(relative_residual < tol * 10.0);
    }
    Ok(())
}

#[test]
fn test_gmres_diagonal_matrix() {
    let n = 4;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, (i + 1) as f64);
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver
        .solve_preconditioned(&a, &b, &precond, &mut x)
        .unwrap();

    let expected = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
    let error = (&x - &expected).norm();
    assert!(error < 1e-9, "Solution error: {error}");
}

#[test]
fn test_gmres_non_symmetric() {
    let n = 3;
    let mut coo = CooMatrix::new(n, n);
    coo.push(0, 0, 4.0);
    coo.push(0, 1, 1.0);
    coo.push(1, 0, 2.0);
    coo.push(1, 1, 5.0);
    coo.push(1, 2, 1.0);
    coo.push(2, 1, 3.0);
    coo.push(2, 2, 6.0);
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![5.0, 8.0, 9.0]);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver
        .solve_preconditioned(&a, &b, &precond, &mut x)
        .unwrap();

    let ax = &a * &x;
    let residual = (&ax - &b).norm();
    assert!(residual < 1e-9, "Residual: {residual}");
}

#[test]
fn test_gmres_zero_rhs() {
    let n = 4;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, (i + 1) as f64);
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::zeros(n);
    let mut x = DVector::from_element(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver
        .solve_preconditioned(&a, &b, &precond, &mut x)
        .unwrap();

    assert!(x.norm() < 1e-14, "Solution should be zero for zero RHS");
}

#[test]
fn test_gmres_restart_dimension() {
    let n = 5;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 4.0);
        if i > 0 {
            coo.push(i, i - 1, 1.0);
        }
        if i < n - 1 {
            coo.push(i, i + 1, 1.0);
        }
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 3);
    let precond = IdentityPreconditioner;

    solver
        .solve_preconditioned(&a, &b, &precond, &mut x)
        .unwrap();

    let ax = &a * &x;
    let residual = (&ax - &b).norm();
    assert!(residual < 1e-9, "Residual: {residual}");
}

#[test]
fn test_gmres_with_initial_guess() {
    let n = 4;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, (i + 1) as f64);
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = DVector::from_vec(vec![0.5, 1.5, 2.5, 3.5]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver
        .solve_preconditioned(&a, &b, &precond, &mut x)
        .unwrap();

    let expected = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
    let error = (&x - &expected).norm();
    assert!(error < 1e-9, "Solution error: {error}");
}

#[test]
fn test_gmres_larger_nonsymmetric() {
    let n = 5;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 5.0);
        if i > 0 {
            coo.push(i, i - 1, 2.0);
        }
        if i < n - 1 {
            coo.push(i, i + 1, 1.0);
        }
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![6.0, 11.0, 11.0, 11.0, 8.0]);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    solver
        .solve_preconditioned(&a, &b, &precond, &mut x)
        .unwrap();

    let ax = &a * &x;
    let residual = (&ax - &b).norm();
    assert!(residual < 1e-9, "Residual: {residual}");
}

#[test]
fn test_gmres_convergence_tight_tolerance() {
    let n = 4;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, (i + 1) as f64);
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-14).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
    assert!(result.is_ok());
}

#[test]
fn test_gmres_max_iterations_exceeded() {
    let n = 4;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, (i + 1) as f64);
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-14).with_max_iterations(1);
    let solver = GMRES::new(config, 1);
    let precond = IdentityPreconditioner;

    let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
    assert!(result.is_err());
}

#[test]
fn test_gmres_dimension_mismatch() {
    let n = 4;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, (i + 1) as f64);
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 4.0]);
    let mut x = DVector::zeros(2);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 10);
    let precond = IdentityPreconditioner;

    let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
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

mod edge_case_tests;
mod extended_edge_case_tests;
