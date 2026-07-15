//! Extended comprehensive edge case and property-based tests for linear solvers
//!
//! This module adds property-based testing with proptest and additional edge cases
//! to complement the existing edge_case_tests.rs file.
//!
//! Standards: Saad (2003) "Iterative Methods for Sparse Linear Systems"

use crate::linear_solver::preconditioners::multigrid::AMGConfig;
#[cfg(test)]
use crate::linear_solver::preconditioners::{
    AlgebraicMultigrid, IdentityPreconditioner, JacobiPreconditioner,
};
use crate::linear_solver::traits::{IterativeLinearSolver, Preconditioner};
use crate::linear_solver::IterativeSolverConfig;
use crate::linear_solver::{BiCGSTAB, ConjugateGradient, GMRES};
use crate::sparse::SparseMatrixBuilder;
use approx::assert_relative_eq;
use leto::Array1;
use leto_ops::CsrMatrix as AtlasSparseMatrix;
use proptest::prelude::*;

use super::{
    array_from_fn, assert_array_close, assert_spmv_matches, basic_preconditioner_matrix,
    fill_array, filled_array, residual_norm,
};

fn amg_matrix_from_solver_matrix(
    matrix: &crate::sparse::SparseMatrix<f64>,
) -> AtlasSparseMatrix<f64> {
    matrix.clone()
}

#[test]
fn test_single_element_system() -> Result<(), Box<dyn std::error::Error>> {
    let n = 1;
    let mut builder = SparseMatrixBuilder::new(n, n);
    builder.add_entry(0, 0, 5.0)?;
    let a = builder.build()?;

    let b = filled_array(n, 10.0);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(10);
    let solver = ConjugateGradient::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    assert_relative_eq!(x[0], 2.0, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_poisson_2d_amg_convergence() -> Result<(), Box<dyn std::error::Error>> {
    let n = 32;
    let mut builder = SparseMatrixBuilder::new(n * n, n * n);

    for i in 0..n {
        for j in 0..n {
            let idx = i * n + j;
            builder.add_entry(idx, idx, 4.0)?;
            if i > 0 {
                builder.add_entry(idx, (i - 1) * n + j, -1.0)?;
            }
            if i < n - 1 {
                builder.add_entry(idx, (i + 1) * n + j, -1.0)?;
            }
            if j > 0 {
                builder.add_entry(idx, i * n + (j - 1), -1.0)?;
            }
            if j < n - 1 {
                builder.add_entry(idx, i * n + (j + 1), -1.0)?;
            }
        }
    }
    let a = builder.build()?;

    let b = filled_array(n * n, 1.0);
    let mut x = Array1::zeros([n * n]);

    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(100);
    let amg_config = AMGConfig {
        max_levels: 5,
        min_coarse_size: 20,
        ..Default::default()
    };
    let amg_a = amg_matrix_from_solver_matrix(&a);
    let amg = AlgebraicMultigrid::new(&amg_a, amg_config)?;

    let solver_gmres = GMRES::new(config, 30);
    solver_gmres.solve(&a, &b, &mut x, Some(&amg))?;

    let residual = residual_norm(&a, &x, &b);
    assert!(
        residual < 1e-7,
        "GMRES + AMG failed to converge, residual: {residual}"
    );

    fill_array(&mut x, 0.0);
    let solver_bicg = BiCGSTAB::new(config);
    solver_bicg.solve(&a, &b, &mut x, Some(&amg))?;

    let residual = residual_norm(&a, &x, &b);
    assert!(
        residual < 1e-7,
        "BiCGSTAB + AMG failed to converge, residual: {residual}"
    );

    Ok(())
}

#[test]
fn test_diagonal_matrix_cg() -> Result<(), Box<dyn std::error::Error>> {
    let n = 10;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        builder.add_entry(i, i, (i + 1) as f64)?;
    }
    let a = builder.build()?;

    let b = array_from_fn(n, |i| (i + 1) as f64);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-12).with_max_iterations(10);
    let solver = ConjugateGradient::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    for i in 0..n {
        assert_relative_eq!(x[i], 1.0, epsilon = 1e-10);
    }
    Ok(())
}

#[test]
fn test_identity_matrix() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        builder.add_entry(i, i, 1.0)?;
    }
    let a = builder.build()?;

    let b = array_from_fn(n, |i| (i + 1) as f64);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-12).with_max_iterations(10);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    assert_array_close(&x, &b, 1e-10);
    Ok(())
}

#[test]
fn test_strongly_dominant_matrix() -> Result<(), Box<dyn std::error::Error>> {
    let n = 10;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        builder.add_entry(i, i, 10.0)?;
        if i > 0 {
            builder.add_entry(i, i - 1, 1.0)?;
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, 1.0)?;
        }
    }
    let a = builder.build()?;

    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    assert_spmv_matches(&a, &x, &b, 1e-8);
    Ok(())
}

#[test]
fn test_all_ones_rhs() -> Result<(), Box<dyn std::error::Error>> {
    let n = 10;
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
    let a = builder.build()?;

    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(1000);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    assert_spmv_matches(&a, &x, &b, 1e-8);
    Ok(())
}

#[test]
fn test_alternating_sign_rhs() -> Result<(), Box<dyn std::error::Error>> {
    let n = 10;
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
    let a = builder.build()?;

    let b = array_from_fn(n, |i| if i % 2 == 0 { 1.0 } else { -1.0 });
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    assert_spmv_matches(&a, &x, &b, 1e-6);
    Ok(())
}

#[test]
fn test_jacobi_uniform_diagonal() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        builder.add_entry(i, i, 3.0)?;
        if i > 0 {
            builder.add_entry(i, i - 1, 0.5)?;
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, 0.5)?;
        }
    }
    let a = builder.build()?;

    let preconditioner_matrix = basic_preconditioner_matrix(&a);
    let precond = JacobiPreconditioner::new(&preconditioner_matrix)?;
    let r = filled_array(n, 6.0);
    let mut z = Array1::zeros([n]);

    precond.apply_to(&r, &mut z)?;

    for i in 0..n {
        assert_relative_eq!(z[i], 2.0, epsilon = 1e-10);
    }
    Ok(())
}

#[test]
fn test_large_sparse_system() -> Result<(), Box<dyn std::error::Error>> {
    let n = 100;
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
    let a = builder.build()?;

    let b = filled_array(n, 1.0);
    let mut x = Array1::zeros([n]);

    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
    let solver = ConjugateGradient::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    let residual = residual_norm(&a, &x, &b);
    assert!(residual < 1e-6, "Residual {residual} exceeds tolerance");
    Ok(())
}

proptest! {
    #[test]
    fn prop_cg_solves_spd_system(
        n in 5..20usize,
        diag in 3.0..10.0f64,
        off_diag in -2.0..2.0f64
    ) {
        prop_assume!(diag > 2.0 * off_diag.abs() + 1e-6);
        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            builder.add_entry(i, i, diag).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, off_diag).unwrap();
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, off_diag).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let b = filled_array(n, 1.0);
        let mut x = Array1::zeros([n]);

        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
        let solver = ConjugateGradient::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity)).unwrap();

        assert_spmv_matches(&a, &x, &b, 1e-6);
    }
}

proptest! {
    #[test]
    fn prop_bicgstab_convergence(
        n in 5..15usize,
        scale in 1.0..5.0f64
    ) {
        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            builder.add_entry(i, i, scale * 3.0).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, -scale).unwrap();
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -scale).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let b = filled_array(n, 1.0);
        let mut x = Array1::zeros([n]);

        let config = IterativeSolverConfig::new(1e-7).with_max_iterations(1000);
        let solver = BiCGSTAB::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity)).unwrap();

        for i in 0..n {
            assert!(x[i].is_finite(), "Solution contains non-finite values");
        }

        let residual = residual_norm(&a, &x, &b);
        assert!(residual < 1e-5, "Residual {residual} too large");
    }
}

proptest! {
    #[test]
    fn prop_jacobi_diagonal_scaling(
        n in 5..15usize,
        diag_val in 1.0..10.0f64
    ) {
        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            builder.add_entry(i, i, diag_val).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, 0.1).unwrap();
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, 0.1).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let preconditioner_matrix = basic_preconditioner_matrix(&a);
        let precond = JacobiPreconditioner::new(&preconditioner_matrix).unwrap();
        let r = filled_array(n, diag_val);
        let mut z = Array1::zeros([n]);

        precond.apply_to(&r, &mut z).unwrap();

        for i in 0..n {
            assert_relative_eq!(z[i], 1.0, epsilon = 1e-10);
        }
    }
}

proptest! {
    #[test]
    fn prop_solution_uniqueness(
        n in 5..12usize,
        seed1 in 0.1..5.0f64,
        seed2 in 0.1..5.0f64
    ) {
        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            builder.add_entry(i, i, 2.0).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0).unwrap();
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -1.0).unwrap();
            }
        }
        let a = builder.build().unwrap();

        let b = filled_array(n, 1.0);

        let mut x1 = filled_array(n, seed1);
        let mut x2 = filled_array(n, seed2);

        let config = IterativeSolverConfig::new(1e-9).with_max_iterations(1000);
        let solver = ConjugateGradient::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x1, Some(&identity)).unwrap();
        solver.solve(&a, &b, &mut x2, Some(&identity)).unwrap();

        for i in 0..n {
            assert_relative_eq!(x1[i], x2[i], epsilon = 1e-6);
        }
    }
}
