use super::*;
use crate::high_order::dg::{
    matrix_cols, matrix_rows, matrix_zeros, vector_from_element, vector_len,
};
use eunomia::assert_relative_eq;
use leto::Array1;

#[test]
fn test_dg_operator_new() {
    let order = 2;
    let num_components = 1;
    let dg_op = DGOperator::new(order, num_components, None).expect("expected value");

    assert_eq!(dg_op.order, order);
    assert_eq!(dg_op.num_components, num_components);
    // Orthogonal basis uses order + 2 quadrature points
    assert_eq!(vector_len(&dg_op.basis.quad_points), order + 2);
    assert_eq!(vector_len(&dg_op.basis.quad_weights), order + 2);
    assert_eq!(matrix_rows(&dg_op.basis.phi), order + 1);
    assert_eq!(matrix_cols(&dg_op.basis.phi), order + 2);
    assert_eq!(matrix_rows(&dg_op.basis.mass_matrix), order + 1);
    assert_eq!(matrix_cols(&dg_op.basis.mass_matrix), order + 1);
    assert_eq!(matrix_rows(&dg_op.basis.stiffness_matrix), order + 1);
    assert_eq!(matrix_cols(&dg_op.basis.stiffness_matrix), order + 1);
    assert_eq!(matrix_rows(&dg_op.basis.diff_matrix), order + 1);
    assert_eq!(matrix_cols(&dg_op.basis.diff_matrix), order + 1);
}

#[test]
fn test_dg_operator_compute_derivative() {
    let order = 2;
    let num_components = 1;
    let dg_op = DGOperator::new(order, num_components, None).expect("expected value");

    // u(x) = x -> u' (x) = 1
    // In Legendre basis: u(x) = P_1(x) -> u' (x) = P_0(x)
    let mut u = matrix_zeros(1, order + 1);
    u[[0, 1]] = 1.0;

    let du_dx = dg_op.compute_derivative(&u).expect("expected value");

    assert_relative_eq!(du_dx[[0, 0]], 1.0, epsilon = 1e-10);
    assert_relative_eq!(du_dx[[0, 1]], 0.0, epsilon = 1e-10);
    assert_relative_eq!(du_dx[[0, 2]], 0.0, epsilon = 1e-10);
}

#[test]
fn test_dg_operator_rhs() {
    let order = 2;
    let num_components = 1;
    let dg_op = DGOperator::new(order, num_components, None).expect("expected value");

    // Test with a constant solution u(x) = 1
    let mut u = matrix_zeros(1, order + 1);
    u[[0, 0]] = 1.0; // P_0(x) = 1

    // For a constant solution, the RHS should be zero
    let flux = |u: &Array1<f64>| u.clone();
    let bc = |_: f64, u: &Array1<f64>, _: bool| u.clone();

    let rhs = dg_op.rhs(&u, flux, bc).expect("expected value");

    for i in 0..=order {
        assert_relative_eq!(rhs[[0, i]], 0.0, epsilon = 1e-10);
    }
}

#[test]
fn test_dg_operator_project() {
    let order = 3;
    let num_components = 1;
    let dg_op = DGOperator::new(order, num_components, None).expect("expected value");

    // Project f(x) = x^2
    let f = |x: f64| vector_from_element(1, x * x);
    let u = dg_op.project(f).expect("expected value");

    // Verify projection at quadrature points
    for (i, &xi) in dg_op.basis.quad_points.iter().enumerate() {
        let mut projected = 0.0;
        for j in 0..=order {
            projected += u[[0, j]] * dg_op.basis.phi[[j, i]];
        }
        assert_relative_eq!(projected, xi * xi, epsilon = 1e-10);
    }

    // Verify L2 error
    let mut l2_error_sq = 0.0;
    for (q, &x) in dg_op.basis.quad_points.iter().enumerate() {
        let w = dg_op.basis.quad_weights[q];
        let mut projected = 0.0;
        for j in 0..=order {
            projected += u[[0, j]] * dg_op.basis.phi[[j, q]];
        }
        let error = projected - x * x;
        l2_error_sq += w * error * error;
    }
    assert!(l2_error_sq.sqrt() < 1e-10);
}

#[test]
fn test_dg_operator_numerical_flux() {
    let order = 2;
    let num_components = 1;
    let dg_op = DGOperator::new(order, num_components, None).expect("expected value");

    let u = vector_from_element(1, 1.0);
    let f = vector_from_element(1, 1.0);
    let bc = |_: f64, u: &Array1<f64>, _: bool| u.clone();

    // Central flux: f_num = 0.5 * (f_l + f_r)
    let f_num = dg_op.compute_numerical_flux(0.0, &f, &u, Some(&u), bc);
    assert_relative_eq!(f_num[0], 1.0, epsilon = 1e-10);
}
