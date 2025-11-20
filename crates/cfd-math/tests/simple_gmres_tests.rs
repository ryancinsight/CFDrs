//! Basic GMRES solver tests to validate core functionality

use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{IterativeLinearSolver, GMRES};
use cfd_math::sparse;
use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

/// Test basic GMRES functionality
#[test]
fn test_gmres_basic() {
    let n = 5;
    let mut coo = nalgebra_sparse::CooMatrix::new(n, n);
    coo.push(0, 0, 2.0);
    coo.push(0, 1, -1.0);
    coo.push(1, 0, -1.0);
    coo.push(1, 1, 2.0);
    coo.push(1, 2, -1.0);
    coo.push(2, 1, -1.0);
    coo.push(2, 2, 3.0);
    coo.push(2, 3, -1.0);
    coo.push(3, 2, -1.0);
    coo.push(3, 3, 2.0);
    coo.push(3, 4, -1.0);
    coo.push(4, 3, -1.0);
    coo.push(4, 4, 2.0);
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 2.0, 1.0]);
    let mut x = DVector::zeros(n);
    let config = cfd_math::linear_solver::IterativeSolverConfig::new(1e-8).with_max_iterations(100);
    let solver = GMRES::new(config, 4); // Restart after 4 iterations

    let result = solver.solve::<IdentityPreconditioner>(&a, &b, &mut x, None);
    assert!(result.is_ok(), "GMRES should converge");

    // Verify residual is reasonable
    let mut ax = DVector::zeros(n);
    sparse::spmv(&a, &x, &mut ax);
    let residual_norm = (&ax - &b).norm();
    assert!(residual_norm < 1.5, "Residual {} too large", residual_norm);
}

/// Test GMRES restart mechanism
#[test]
fn test_gmres_restart() {
    let n = 8;
    let mut coo = nalgebra_sparse::CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 2.0);
        if i > 0 {
            coo.push(i, i - 1, -1.0);
        }
        if i < n - 1 {
            coo.push(i, i + 1, -1.0);
        }
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);

    // Use small restart dimension to force restart testing
    let config = cfd_math::linear_solver::IterativeSolverConfig::new(1e-8).with_max_iterations(200);
    let solver = GMRES::new(config, 3); // Small restart dimension

    let result = solver.solve::<IdentityPreconditioner>(&a, &b, &mut x, None);
    assert!(result.is_ok(), "GMRES with restart should converge");
}

/// Test GMRES with preconditioner
#[test]
fn test_gmres_with_preconditioner() {
    let n = 6;
    let mut coo = nalgebra_sparse::CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 1.0 + (i as f64) * 0.1); // Diagonally dominant
    }
    let a = CsrMatrix::from(&coo);

    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);
    let config = cfd_math::linear_solver::IterativeSolverConfig::new(1e-10).with_max_iterations(50);
    let solver = GMRES::new(config, n);
    let precond = IdentityPreconditioner;

    let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
    assert!(result.is_ok(), "GMRES with preconditioner should work");

    let final_residual = (&a * &x - &b).norm();
    assert!(
        final_residual < 1e-8,
        "Final residual {} too large",
        final_residual
    );
}
