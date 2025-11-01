//! Tests for the Discontinuous Galerkin module

use super::*;
use approx::assert_relative_eq;
use nalgebra::DVector;
use std::f64::consts::PI;

#[test]
fn test_dg_basis_legendre() {
    // Test Legendre polynomial evaluation
    let basis = DGBasis::new(3, BasisType::Legendre).unwrap();
    
    // Test at x = 0
    let x = 0.0;
    let phi = basis.evaluate(x);
    
    // P0(0) = 1/√2, P1(0) = 0, P2(0) = -√(5/8), P3(0) = 0
    assert_relative_eq!(phi[0], 1.0 / 2.0_f64.sqrt(), epsilon = 1e-10);
    assert_relative_eq!(phi[1], 0.0, epsilon = 1e-10);
    assert_relative_eq!(phi[2], - (5.0/8.0).sqrt() / 2.0, epsilon = 1e-10);
    assert_relative_eq!(phi[3], 0.0, epsilon = 1e-10);
}

#[test]
fn test_dg_basis_lagrange() {
    // Test Lagrange polynomial evaluation at nodes
    let basis = DGBasis::new(2, BasisType::Lagrange).unwrap();
    let nodes = basis.nodes();
    
    // Test each basis function at each node
    for (i, &xi) in nodes.iter().enumerate() {
        let phi = basis.evaluate(xi);
        for (j, &phi_j) in phi.iter().enumerate() {
            // Should be 1 at node i, 0 at other nodes
            let expected = if i == j { 1.0 } else { 0.0 };
            assert_relative_eq!(phi_j, expected, epsilon = 1e-10);
        }
    }
}

#[test]
fn test_dg_operator_mass_matrix() {
    // Test that the mass matrix is symmetric and positive definite
    let order = 2;
    let dg_op = DGOperator::new(order, 1, None).unwrap();
    
    // Check symmetry
    for i in 0..=order {
        for j in 0..=order {
            assert_relative_eq!(
                dg_op.mass_matrix()[(i, j)],
                dg_op.mass_matrix()[(j, i)],
                epsilon = 1e-10,
                "Mass matrix not symmetric at ({}, {})",
                i, j
            );
        }
    }
    
    // Check positive definiteness by computing the Cholesky decomposition
    let chol = dg_op.mass_matrix().cholesky();
    assert!(chol.is_some(), "Mass matrix is not positive definite");
}

#[test]
fn test_dg_operator_derivative() {
    // Test that the derivative of a polynomial is computed correctly
    let order = 3;
    let dg_op = DGOperator::new(order, 1, None).unwrap();
    
    // Test function: u(x) = x^3 - x
    // Derivative: u'(x) = 3x^2 - 1
    let u = DVector::from_fn(order + 1, |i, _| {
        let x = dg_op.nodes()[i];
        x.powi(3) - x
    });
    
    let du_dx = dg_op.compute_derivative(&u).unwrap();
    
    // Check at each node
    for i in 0..=order {
        let x = dg_op.nodes()[i];
        let expected = 3.0 * x.powi(2) - 1.0;
        assert_relative_eq!(
            du_dx[i],
            expected,
            epsilon = 1e-10,
            "Derivative mismatch at x = {}",
            x
        );
    }
}

#[test]
fn test_dg_solver_linear_advection() -> Result<(), Box<dyn std::error::Error>> {
    // Test solving the linear advection equation
    let order = 2;
    let t_final = 1.0;
    let cfl = 0.1;
    
    // Set up DG operator with upwind flux
    let params = DGOperatorParams::new()
        .with_volume_flux(FluxType::Upwind)
        .with_surface_flux(FluxType::Upwind);
        
    let dg_op = DGOperator::new(order, 1, Some(params))?;
    
    // Set up solver
    let mut solver = DGSolver::new(
        dg_op,
        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
            .with_t_final(t_final)
            .with_cfl(cfl),
    )?;
    
    // Initial condition: u(x,0) = sin(πx)
    solver.initialize(|x| DVector::from_vec(vec![(PI * x).sin()]))?;
    
    // Solve
    solver.solve(
        |t, u| {
            // For linear advection: du/dt = -du/dx
            Ok(-solver.operator.compute_derivative(u)?)
        },
        None::<fn(_, _) -> _>,
    )?;
    
    // Check the solution at t = t_final
    // Exact solution: u(x,t) = sin(π(x - t))
    let (l2_error, l2_norm) = solver.compute_error(|x| {
        DVector::from_vec(vec![(PI * (x - t_final)).sin()])
    })?;
    
    // The error should be small
    assert!(
        l2_error < 1e-4,
        "L2 error {} is too large (expected < 1e-4)",
        l2_error
    );
    
    Ok(())
}
