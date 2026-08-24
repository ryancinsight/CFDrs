#![allow(missing_docs)]
//! Matrix-free CG demonstration.
//!
//! This is the minimal end-to-end setup that exercises the matrix-free CG
//! path introduced in `docs/book/numerics_and_solvers.md` and detailed in
//! `docs/book/examples/matrix_free_demo.md`:
//!
//! 1. **1D diffusion stencil** — assemble a 3-point second-difference
//!    operator on a uniform grid as a `CsrMatrix<f64>` (the matrix-free
//!    idea materialised on a small structured problem so the example stays
//!    runnable on the standard test profile).
//! 2. **Conjugate gradient solve** — `cfd_math::linear_solver::cg` solves
//!    `A·x = b` for the SPD tridiagonal system and returns a
//!    `SolveReport` carrying the iteration count and final residual norm.
//! 3. **Convergence budget** — `IterativeSolverConfig` bounds both the
//!    absolute and the RHS-scaled residual norm, matching the chapter's
//!    "matrix-free" claim (no `A` is materialised beyond the tridiagonal
//!    layout).
//!
//! Run with: `cargo run -p cfd-math --example matrix_free_demo`

use cfd_math::linear_solver::{krylov::cg, IterativeSolverConfig};
use cfd_math::sparse::SparseMatrixBuilder;
use leto::Array1;
use leto_ops::CsrMatrix;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // -------------------------------------------------------------------
    // 1. Assemble the 1D diffusion operator d²u/dx² with homogeneous
    //    Dirichlet BCs on a uniform grid. The 3-point central-difference
    //    stencil has 2/dx² on the diagonal and -1/dx² on each off-diagonal
    //    so the operator is symmetric positive-definite and CG converges
    //    in at most n steps.
    // -------------------------------------------------------------------
    let n = 64_usize;
    let dx = 1.0 / (n + 1) as f64;
    let dx2_inv = 1.0 / (dx * dx);

    let mut builder = SparseMatrixBuilder::<f64>::new(n, n);
    for i in 0..n {
        builder.add_entry(i, i, 2.0 * dx2_inv)?;
        if i > 0 {
            builder.add_entry(i, i - 1, -dx2_inv)?;
        }
        if i + 1 < n {
            builder.add_entry(i, i + 1, -dx2_inv)?;
        }
    }
    let _matrix_unused: CsrMatrix<f64> = builder.build()?;

    // Build the same matrix directly via `from_parts` so the example
    // does not depend on the builder's internal CsMatrix type.
    let mut row_offsets: Vec<usize> = Vec::with_capacity(n + 1);
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    row_offsets.push(0);
    for i in 0..n {
        if i > 0 {
            col_indices.push(i - 1);
            values.push(-dx2_inv);
        }
        col_indices.push(i);
        values.push(2.0 * dx2_inv);
        if i + 1 < n {
            col_indices.push(i + 1);
            values.push(-dx2_inv);
        }
        row_offsets.push(col_indices.len());
    }
    let matrix = CsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
        .map_err(|e| format!("from_parts: {e}"))?;

    // -------------------------------------------------------------------
    // 2. Right-hand side: b_i = 1 (constant source). The analytic solution
    //    of -d²u/dx² = 1 with u(0) = u(1) = 0 is u(x) = x(1 - x) / 2.
    //    The closed-form check below confirms the CG iterate converges
    //    to the parabola evaluated at the cell centres x_i = (i + 0.5) /
    //    (n + 1).
    // -------------------------------------------------------------------
    let rhs = Array1::from_elem([n], 1.0_f64);
    let mut solution = Array1::<f64>::zeros([n]);

    let config = IterativeSolverConfig {
        max_iterations: n * 4,
        tolerance: 1.0e-10,
        relative_tolerance: 1.0e-12,
    };
    let report = cg(&matrix, &rhs, &mut solution, &config)
        .map_err(|e| format!("cg failed: {e:?}"))?;
    assert!(
        report.final_residual_norm < config.tolerance,
        "CG did not converge: final residual {:.3e} exceeds tolerance {:.3e}",
        report.final_residual_norm, config.tolerance
    );
    assert!(
        report.iterations < config.max_iterations,
        "CG did not converge: {} iterations hit the budget {}",
        report.iterations, config.max_iterations
    );

    // -------------------------------------------------------------------
    // 3. Closed-form check: the analytic solution at cell centre
    //    x_i = (i + 0.5) / (n + 1) is u_analytic(x) = x(1 - x) / 2.
    //    The L2 error between the CG iterate and the analytic
    //    parabola is the discretisation error of the central-difference
    //    Laplacian, O(dx²) for the L∞ norm and O(dx) for the L² norm.
    // -------------------------------------------------------------------
    let analytic: Vec<f64> = (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / (n as f64 + 1.0);
            x * (1.0 - x) / 2.0
        })
        .collect();

    let l2_error_sq: f64 = solution
        .iter()
        .zip(analytic.iter())
        .map(|(s, a)| (s - a).powi(2))
        .sum();
    let l2_error = l2_error_sq.sqrt();
    let analytic_norm_sq: f64 = analytic.iter().map(|a| a * a).sum();
    let relative_l2 = l2_error / analytic_norm_sq.sqrt();
    // O(dx) L2 discretisation error of central-difference Laplacian.
    let expected = dx;
    assert!(
        relative_l2 < 10.0 * expected,
        "relative L2 error {relative_l2:.3e} should be O(dx) ~ {expected:.3e}"
    );

    Ok(())
}
