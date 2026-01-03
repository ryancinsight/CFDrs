//! Tests for matrix-free linear solvers.

use crate::linear_solver::operators::{IdentityOperator, ScaledOperator};
use crate::linear_solver::preconditioners::IdentityPreconditioner;
use crate::linear_solver::{ConjugateGradient, GMRES, IterativeSolverConfig};
use crate::linear_solver::traits::IterativeLinearSolver;
use nalgebra::DVector;
use approx::assert_relative_eq;

#[test]
fn test_matrix_free_cg_identity() {
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);
    let operator = IdentityOperator;

    let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut x = DVector::zeros(5);

    solver.solve(&operator, &b, &mut x, None::<&IdentityPreconditioner>).unwrap();

    for i in 0..5 {
        assert_relative_eq!(x[i], b[i], epsilon = 1e-8);
    }
}

#[test]
fn test_matrix_free_gmres_identity() {
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 5);
    let operator = IdentityOperator;

    let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut x = DVector::zeros(5);

    solver.solve(&operator, &b, &mut x, None::<&IdentityPreconditioner>).unwrap();

    for i in 0..5 {
        assert_relative_eq!(x[i], b[i], epsilon = 1e-8);
    }
}

#[test]
fn test_scaled_operator_integration() {
    let base_op = IdentityOperator;
    let scaled_op = ScaledOperator::new(&base_op, 2.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);

    let b = DVector::from_vec(vec![2.0, 4.0, 6.0]); // Should give solution [1, 2, 3]
    let mut x = DVector::zeros(3);

    solver.solve(&scaled_op, &b, &mut x, None::<&IdentityPreconditioner>).unwrap();

    assert_relative_eq!(x[0], 1.0, epsilon = 1e-8);
    assert_relative_eq!(x[1], 2.0, epsilon = 1e-8);
    assert_relative_eq!(x[2], 3.0, epsilon = 1e-8);
}

#[test]
fn test_operator_size_mismatch() {
    // IdentityOperator returns size 0 currently, but solvers check b.len() vs operator.size()
    // Let's use a real operator like LaplacianOperator2D
    use crate::linear_solver::operators::LaplacianOperator2D;
    let operator = LaplacianOperator2D::new(2, 2, 1.0, 1.0);
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);

    let b = DVector::from_vec(vec![1.0, 2.0]); // Wrong size (expected 4)
    let mut x = DVector::zeros(4);

    assert!(solver.solve(&operator, &b, &mut x, None::<&IdentityPreconditioner>).is_err());
}
