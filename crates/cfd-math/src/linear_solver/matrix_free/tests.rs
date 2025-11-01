//! Tests for matrix-free linear solvers.

use super::*;
use super::operator::{IdentityOperator, ScaledOperator};
use super::traits::MatrixFreeSolver;
use crate::linear_solver::config::IterativeSolverConfig;
use approx::assert_relative_eq;

#[test]
fn test_matrix_free_cg_identity() {
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = MatrixFreeCG::new(config);
    let operator = IdentityOperator::new(5);

    let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let mut x = vec![0.0; 5];

    solver.solve(&operator, &b, &mut x).unwrap();

    for i in 0..5 {
        assert_relative_eq!(x[i], b[i], epsilon = 1e-8);
    }
}

#[test]
fn test_matrix_free_gmres_identity() {
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = MatrixFreeGMRES::new(config, 5);
    let operator = IdentityOperator::new(5);

    let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let mut x = vec![0.0; 5];

    solver.solve(&operator, &b, &mut x).unwrap();

    for i in 0..5 {
        assert_relative_eq!(x[i], b[i], epsilon = 1e-8);
    }
}

#[test]
fn test_scaled_operator_integration() {
    let base_op = IdentityOperator::new(3);
    let scaled_op = ScaledOperator::new(&base_op, 2.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = MatrixFreeCG::new(config);

    let b = vec![2.0, 4.0, 6.0]; // Should give solution [1, 2, 3]
    let mut x = vec![0.0; 3];

    solver.solve(&scaled_op, &b, &mut x).unwrap();

    assert_relative_eq!(x[0], 1.0, epsilon = 1e-8);
    assert_relative_eq!(x[1], 2.0, epsilon = 1e-8);
    assert_relative_eq!(x[2], 3.0, epsilon = 1e-8);
}

#[test]
fn test_operator_size_mismatch() {
    let operator = IdentityOperator::new(3);
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = MatrixFreeCG::new(config);

    let b = vec![1.0, 2.0]; // Wrong size
    let mut x = vec![0.0; 3];

    assert!(solver.solve(&operator, &b, &mut x).is_err());
}
