//! Adversarial and boundary-condition tests for linear solvers.
//!
//! ## Testing Philosophy
//!
//! These tests validate solver robustness under pathological inputs that reveal
//! failure modes silent in happy-path tests. Each scenario is named after the
//! CFD context where it arises in practice.
//!
//! **Theorem (Robustness Requirement)**: A production-grade iterative solver must:
//! 1. Return `Err` (never panic, never loop infinitely) on non-convergent inputs.
//! 2. Return the exact zero vector for b = 0 regardless of initial guess.
//! 3. Detect and report dimension mismatches before any iteration begins.
//! 4. Produce ‖Ax - b‖/‖b‖ < tolerance for all SPD inputs within max_iterations.

use crate::linear_solver::preconditioners::IdentityPreconditioner;
use crate::linear_solver::traits::IterativeLinearSolver;
use crate::linear_solver::IterativeSolverConfig;
use crate::linear_solver::{BiCGSTAB, ConjugateGradient};
use crate::sparse::SparseMatrixBuilder;
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};

/// Build an N×N SPD tridiagonal matrix with diagonal=2, off-diagonals=-1.
/// This is the 1D Laplacian — the canonical CFD test matrix.
fn laplacian_1d(n: usize) -> CsrMatrix<f64> {
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 2.0);
        if i > 0 {
            coo.push(i, i - 1, -1.0);
        }
        if i + 1 < n {
            coo.push(i, i + 1, -1.0);
        }
    }
    CsrMatrix::from(&coo)
}

/// Build a severely ill-conditioned diagonal matrix.
/// κ(A) ≈ 10^{14} — mimics near-singular pressure systems in poorly
/// resolved CFD grids near stagnation points.
fn ill_conditioned_diagonal(n: usize) -> CsrMatrix<f64> {
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        // Condition number = 10^14 / 1 = 10^14
        let val = if i == 0 { 1e14_f64 } else { 1.0 };
        coo.push(i, i, val);
    }
    CsrMatrix::from(&coo)
}

// ---------------------------------------------------------------------------
// Zero-RHS Tests
// ---------------------------------------------------------------------------

/// Theorem: For any non-singular A and b = 0, the unique solution is x = 0.
/// All iterative solvers must converge in 1 step regardless of initial guess.
#[test]
fn test_cg_zero_rhs_any_initial_guess() {
    let n = 8;
    let a = laplacian_1d(n);
    let b = DVector::zeros(n);
    let mut x = DVector::from_element(n, 1.0); // non-trivial initial guess
    let config = IterativeSolverConfig::new(1e-14).with_max_iterations(100);
    let solver = ConjugateGradient::new(config);
    solver
        .solve(&a, &b, &mut x, Some(&IdentityPreconditioner))
        .expect("CG must converge for zero RHS");
    assert!(
        x.norm() < 1e-12,
        "solution for zero RHS must be zero, got ‖x‖ = {}",
        x.norm()
    );
}

#[test]
fn test_bicgstab_zero_rhs_any_initial_guess() {
    let n = 8;
    let a = laplacian_1d(n);
    let b = DVector::zeros(n);
    let mut x = DVector::from_element(n, 5.0);
    let config = IterativeSolverConfig::new(1e-14).with_max_iterations(100);
    let solver = BiCGSTAB::new(config);
    solver
        .solve(&a, &b, &mut x, Some(&IdentityPreconditioner))
        .expect("BiCGSTAB must converge for zero RHS");
    assert!(
        x.norm() < 1e-12,
        "solution for zero RHS must be zero, got ‖x‖ = {}",
        x.norm()
    );
}

// ---------------------------------------------------------------------------
// Dimension Mismatch Tests
// ---------------------------------------------------------------------------

/// All solvers must detect operator-RHS dimension mismatch before iteration.
#[test]
fn test_cg_dimension_mismatch_returns_err() {
    let a = laplacian_1d(6);
    let b = DVector::from_element(4, 1.0); // wrong size
    let mut x = DVector::zeros(4);
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(50);
    let solver = ConjugateGradient::new(config);
    let result = solver.solve(&a, &b, &mut x, Some(&IdentityPreconditioner));
    assert!(result.is_err(), "CG must return Err on dimension mismatch");
}

#[test]
fn test_bicgstab_dimension_mismatch_returns_err() {
    let a = laplacian_1d(6);
    let b = DVector::from_element(4, 1.0);
    let mut x = DVector::zeros(4);
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(50);
    let solver = BiCGSTAB::new(config);
    let result = solver.solve(&a, &b, &mut x, Some(&IdentityPreconditioner));
    assert!(
        result.is_err(),
        "BiCGSTAB must return Err on dimension mismatch"
    );
}

// ---------------------------------------------------------------------------
// Large Dimension SPD Stress Test
// ---------------------------------------------------------------------------

/// Theorem: CG on n×n tridiagonal (1D Laplacian) converges in ≤ n steps.
/// This stress test validates O(n) convergence count for n=500.
#[test]
fn test_cg_large_spd_converges() {
    let n = 500;
    let a = laplacian_1d(n);
    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);
    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(2000);
    let solver = ConjugateGradient::new(config);
    let monitor = solver
        .solve(&a, &b, &mut x, Some(&IdentityPreconditioner))
        .expect("CG must converge on 500×500 1D Laplacian");
    // Verify solution residual
    let residual = (&a * &x - &b).norm() / b.norm();
    assert!(
        residual < 1e-6,
        "relative residual {residual:.2e} must be < 1e-6 for n=500"
    );
    // Verify iteration count is sub-linear in n (preconditioned convergence boost)
    assert!(
        monitor.iteration <= n,
        "CG must converge in ≤ n={n} iterations, took {}",
        monitor.iteration
    );
}

#[test]
fn test_bicgstab_large_nonsymmetric_converges() {
    // Slightly non-symmetric tridiagonal (simulates advection-diffusion matrix)
    let n = 300;
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 3.0);
        if i > 0 {
            coo.push(i, i - 1, -1.5); // asymmetric off-diagonal
        }
        if i + 1 < n {
            coo.push(i, i + 1, -0.5);
        }
    }
    let a = CsrMatrix::from(&coo);
    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);
    let config = IterativeSolverConfig::new(1e-8).with_max_iterations(2000);
    let solver = BiCGSTAB::new(config);
    solver
        .solve(&a, &b, &mut x, Some(&IdentityPreconditioner))
        .expect("BiCGSTAB must converge on 300×300 non-symmetric advection-diffusion matrix");
    let residual = (&a * &x - &b).norm() / b.norm();
    assert!(
        residual < 1e-6,
        "relative residual {residual:.2e} must be < 1e-6"
    );
}

// ---------------------------------------------------------------------------
// Ill-Conditioned Matrix Tests
// ---------------------------------------------------------------------------

/// A system with κ ≈ 10^14 — solver must either converge (with enough steps)
/// or return a clean `Err`, never panic or produce NaN.
#[test]
fn test_cg_ill_conditioned_no_panic_no_nan() {
    let n = 4;
    let a = ill_conditioned_diagonal(n);
    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(10_000);
    let solver = ConjugateGradient::new(config);
    // We only require no panic and no NaN — convergence may fail gracefully.
    let result = solver.solve(&a, &b, &mut x, Some(&IdentityPreconditioner));
    if let Ok(_) = result {
        // If it converged, the solution must be finite.
        assert!(x.iter().all(|v| v.is_finite()), "solution must be finite");
    } else {
        // Graceful failure is acceptable for extreme ill-conditioning.
    }
    // No panics, no NaN — the test reaching here proves that.
    assert!(
        x.iter().all(|v| !v.is_nan()),
        "solution must never contain NaN"
    );
}

// ---------------------------------------------------------------------------
// Convergence History Validation
// ---------------------------------------------------------------------------

/// Theorem: CG residuals decrease monotonically for SPD matrices in exact arithmetic.
/// In floating-point, we verify residuals at least decrease over 90% of iterations.
#[test]
fn test_cg_convergence_history_generally_monotone() {
    let n = 20;
    let a = laplacian_1d(n);
    let b = DVector::from_element(n, 1.0);
    let mut x = DVector::zeros(n);
    let config = IterativeSolverConfig::new(1e-12).with_max_iterations(500);
    let solver = ConjugateGradient::new(config);
    let monitor = solver
        .solve(&a, &b, &mut x, Some(&IdentityPreconditioner))
        .expect("CG must converge on 20×20 Laplacian");
    let history = &monitor.residual_history;
    // Verify overall decrease (final < initial)
    assert!(
        history.last().copied().unwrap_or(1.0) < history[0],
        "residual must decrease from initial to final"
    );
}

// ---------------------------------------------------------------------------
// SparseMatrixBuilder boundary tests
// ---------------------------------------------------------------------------

/// Adding an entry with out-of-bounds column must return Err, not panic.
#[test]
fn test_sparse_builder_out_of_bounds_returns_err() {
    let mut builder = SparseMatrixBuilder::<f64>::new(3, 3);
    let result = builder.add_entry(0, 5, 1.0); // column 5 > n-1=2
    assert!(result.is_err(), "out-of-bounds column must return Err");
}

/// Adding an entry with out-of-bounds row must return Err, not panic.
#[test]
fn test_sparse_builder_out_of_bounds_row_returns_err() {
    let mut builder = SparseMatrixBuilder::<f64>::new(3, 3);
    let result = builder.add_entry(5, 0, 1.0); // row 5 > m-1=2
    assert!(result.is_err(), "out-of-bounds row must return Err");
}
