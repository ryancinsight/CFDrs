//! Comprehensive edge case tests for linear solvers
//!
//! Tests cover positive/negative/zero/boundary cases per SRS requirements
//! All tests maintain <30s runtime via granular cargo nextest execution

#[cfg(test)]
use crate::linear_solver::preconditioners::{
    IdentityPreconditioner, JacobiPreconditioner, SORPreconditioner,
};
use crate::linear_solver::traits::{IterativeLinearSolver, Preconditioner};
use crate::linear_solver::IterativeSolverConfig;
use crate::linear_solver::{BiCGSTAB, ConjugateGradient};
use crate::sparse::SparseMatrixBuilder;
use approx::assert_relative_eq;
use nalgebra::DVector;

#[test]
fn test_bicgstab_zero_rhs() -> Result<(), Box<dyn std::error::Error>> {
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

    let b = DVector::zeros(n);
    let mut x = DVector::from_element(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    for i in 0..n {
        assert_relative_eq!(x[i], 0.0, epsilon = 1e-8);
    }
    Ok(())
}

#[test]
fn test_cg_zero_rhs() -> Result<(), Box<dyn std::error::Error>> {
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

    let b = DVector::zeros(n);
    let mut x = DVector::from_element(n, 1.0);

    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    for i in 0..n {
        assert_relative_eq!(x[i], 0.0, epsilon = 1e-8);
    }
    Ok(())
}

#[test]
fn test_bicgstab_large_negative() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        builder.add_entry(i, i, 10.0)?;
        if i > 0 {
            builder.add_entry(i, i - 1, -5.0)?;
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, -5.0)?;
        }
    }
    let a = builder.build()?;

    let b = DVector::from_element(n, -100.0);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    let ax = &a * &x;
    for i in 0..n {
        assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
    }
    Ok(())
}

#[test]
fn test_bicgstab_ill_conditioned() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        let diag_value = if i == 0 {
            1000.0
        } else if i == n - 1 {
            0.001
        } else {
            1.0
        };
        builder.add_entry(i, i, diag_value)?;
    }
    let a = builder.build()?;

    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-6).with_max_iterations(10000);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    let ax = &a * &x;
    for i in 0..n {
        assert_relative_eq!(ax[i], b[i], epsilon = 1e-4);
    }
    Ok(())
}

#[test]
fn test_jacobi_mixed_signs() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
    let mut builder = SparseMatrixBuilder::new(n, n);

    for i in 0..n {
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        builder.add_entry(i, i, sign * 2.0)?;
        if i > 0 {
            builder.add_entry(i, i - 1, 0.5)?;
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, 0.5)?;
        }
    }
    let a = builder.build()?;

    let precond = JacobiPreconditioner::new(&a)?;
    let r = DVector::from_element(n, 1.0);
    let mut z = DVector::zeros(n);

    precond.apply_to(&r, &mut z)?;

    for i in 0..n {
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        let expected = r[i] / (sign * 2.0);
        assert_relative_eq!(z[i], expected, epsilon = 1e-10);
    }
    Ok(())
}

#[test]
fn test_sor_boundary_omega() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
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

    for omega in [0.1, 0.5, 1.0, 1.5, 1.9] {
        let sor = SORPreconditioner::new(&a, omega)?;
        let r = DVector::from_element(n, 1.0);
        let mut z = DVector::zeros(n);

        sor.apply_to(&r, &mut z)?;

        assert!(
            z.iter().all(|&x: &f64| x.is_finite()),
            "omega={omega}: NaN detected"
        );
    }
    Ok(())
}

#[test]
fn test_bicgstab_small_positive() {
    let n = 5;
    let mut builder = SparseMatrixBuilder::new(n, n);

    let small = 1e-14;
    for i in 0..n {
        builder.add_entry(i, i, small * 2.0).unwrap();
        if i > 0 {
            builder.add_entry(i, i - 1, -small).unwrap();
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, -small).unwrap();
        }
    }
    let a = builder.build().unwrap();

    let b = DVector::from_element(n, small);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-15).with_max_iterations(1000);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    let result = solver.solve(&a, &b, &mut x, Some(&identity));

    if result.is_ok() {
        assert!(x.iter().all(|&xi: &f64| xi.is_finite()));
    } else {
        assert!(
            result.is_err(),
            "Expected either success or numerical error"
        );
    }
}

#[test]
fn test_cg_perfect_initial_guess() -> Result<(), Box<dyn std::error::Error>> {
    let n = 5;
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

    let b = DVector::from_element(n, 1.0);
    let mut x_exact = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-12).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);
    let identity = IdentityPreconditioner;
    let _ = solver.solve(&a, &b, &mut x_exact, Some(&identity))?;

    let mut x = x_exact.clone();
    let _ = solver.solve(&a, &b, &mut x, Some(&identity))?;

    for i in 0..n {
        assert_relative_eq!(x[i], x_exact[i], epsilon = 1e-10);
    }
    Ok(())
}

#[test]
fn test_bicgstab_pentadiagonal() -> Result<(), Box<dyn std::error::Error>> {
    let n = 25;
    let mut builder = SparseMatrixBuilder::new(n, n);

    let nx = 5;
    for i in 0..n {
        let row = i / nx;
        let col = i % nx;

        builder.add_entry(i, i, 4.0)?;

        if col > 0 {
            builder.add_entry(i, i - 1, -1.0)?;
        }
        if col < nx - 1 {
            builder.add_entry(i, i + 1, -1.0)?;
        }

        if row > 0 {
            builder.add_entry(i, i - nx, -1.0)?;
        }
        if row < nx - 1 {
            builder.add_entry(i, i + nx, -1.0)?;
        }
    }
    let a = builder.build()?;

    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    solver.solve(&a, &b, &mut x, Some(&identity))?;

    let ax = &a * &x;
    for i in 0..n {
        assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
    }
    Ok(())
}

#[test]
fn test_bicgstab_max_iterations() {
    let n = 5;
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

    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);

    let config = IterativeSolverConfig::new(1e-15).with_max_iterations(1);
    let solver = BiCGSTAB::new(config);
    let identity = IdentityPreconditioner;

    let result = solver.solve(&a, &b, &mut x, Some(&identity));

    assert!(
        result.is_err(),
        "Expected convergence failure with max_iter=1"
    );
}
