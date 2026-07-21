//! Tests for matrix-free linear solvers.

use crate::linear_solver::operators::{IdentityOperator, ScaledOperator};
use crate::linear_solver::preconditioners::IdentityPreconditioner;
use crate::linear_solver::traits::IterativeLinearSolver;
use crate::linear_solver::{ConjugateGradient, IterativeSolverConfig, GMRES};
use aequitas::systems::si::{quantities::Length, units::Meter};
use approx::assert_relative_eq;
use cfd_core::error::Error;
use leto::{Array1, BoundaryCondition};

fn array(values: Vec<f64>) -> Array1<f64> {
    Array1::from_shape_vec([values.len()], values).expect("valid Leto vector shape")
}

fn assert_array_close(actual: &Array1<f64>, expected: &Array1<f64>, epsilon: f64) {
    assert_eq!(actual.shape(), expected.shape());
    for idx in 0..actual.shape()[0] {
        assert_relative_eq!(actual[idx], expected[idx], epsilon = epsilon);
    }
}

#[test]
fn test_matrix_free_cg_identity() {
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);
    let operator = IdentityOperator;

    let b = array(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut x = Array1::zeros([5]);

    solver
        .solve(&operator, &b, &mut x, None::<&IdentityPreconditioner>)
        .unwrap();

    assert_array_close(&x, &b, 1e-8);
}

#[test]
fn test_matrix_free_gmres_identity() {
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = GMRES::new(config, 5);
    let operator = IdentityOperator;

    let b = array(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let mut x = Array1::zeros([5]);

    solver
        .solve(&operator, &b, &mut x, None::<&IdentityPreconditioner>)
        .unwrap();

    assert_array_close(&x, &b, 1e-8);
}

#[test]
fn test_scaled_operator_integration() {
    let base_op = IdentityOperator;
    let scaled_op = ScaledOperator::new(&base_op, 2.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);

    let b = array(vec![2.0, 4.0, 6.0]); // Should give solution [1, 2, 3]
    let mut x = Array1::zeros([3]);

    solver
        .solve(&scaled_op, &b, &mut x, None::<&IdentityPreconditioner>)
        .unwrap();

    assert_relative_eq!(x[0], 1.0, epsilon = 1e-8);
    assert_relative_eq!(x[1], 2.0, epsilon = 1e-8);
    assert_relative_eq!(x[2], 3.0, epsilon = 1e-8);
}

#[test]
fn test_operator_size_mismatch() {
    // IdentityOperator returns size 0 currently, but solvers check b.len() vs operator.size()
    // Let's use a real operator like LaplacianOperator2D
    use crate::linear_solver::operators::LaplacianOperator2D;
    let operator = LaplacianOperator2D::new(
        2,
        2,
        Length::from_unit::<Meter>(1.0),
        Length::from_unit::<Meter>(1.0),
        BoundaryCondition::Dirichlet,
    )
    .expect("valid Laplacian grid");
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);

    let b = array(vec![1.0, 2.0]); // Wrong size (expected 4)
    let mut x = Array1::zeros([2]);

    let result = solver.solve(&operator, &b, &mut x, None::<&IdentityPreconditioner>);
    match result {
        Err(Error::InvalidConfiguration(message)) => {
            assert_eq!(message, "Operator size (4) doesn't match RHS vector (2)");
        }
        Err(error) => panic!("expected invalid configuration, got {error:?}"),
        Ok(_) => panic!("expected operator size mismatch"),
    }
}
