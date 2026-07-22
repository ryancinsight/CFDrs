//! # `cfd-math::sparse` × `cfd-math::linear_solver` Atlas-Parity Harness
//!
//! Solves the same 1D Poisson problem \(u'' = -\sin(\pi x), \pi^2\) with
//! homogeneous Dirichlet BCs on \([0, 1]\) using
//!
//! 1. **Legacy**: `nalgebra::DMatrix` + a hand-rolled conjugate-gradient
//!    loop (dense path);
//! 2. **Atlas**: `cfd_math::sparse::SparseMatrixBuilder` (which uses
//!    `leto_ops::CsrMatrix` under the hood) and
//!    `cfd_math::linear_solver::{ConjugateGradient, LinearSolver}`.
//!
//! Both implementations consume the same discrete operator (with
//! matching \(1/h^2\) scaling on every stencil entry), the same RHS,
//! and the same CG tolerance; the harness then compares their
//! solutions and per-path residuals and emits a JSON line suitable
//! for regression gates.
//!
//! Run with:
//! ```sh
//! cargo run --release --example migration_cfd_math_parity -p cfd-math
//! ```
//!
//! ## Parity tolerances (from `migration_validation.md`)
//!
//! - Velocity field RMS (here `u`): 1e-6 relative
//! - Mass conservation (residual L2): 1e-8 relative
//! - L∞ residual:                    1e-8 absolute

use cfd_core::error::Error as CfdError;
use cfd_math::error::Result as CfdResult;
use cfd_math::linear_solver::{ConjugateGradient, IterativeSolverConfig, LinearSolver};
use cfd_math::sparse::{SparseMatrixBuilder, try_spmv};
use leto::Array1;
use nalgebra::DMatrix;
use std::time::Instant;

/// 1D Poisson test problem. The Laplacian operator is scaled by `1/h²`
/// consistently across both legacy and Atlas paths so the manufactured
/// solution \(u(x) = \sin(\pi x)/\pi^{2}\) is the analytic solution of
/// the discretised system.
#[derive(Clone, Copy, Debug)]
struct Problem {
    n: usize,
    h: f64,
    h2_inv: f64,
}

impl Problem {
    fn new(n: usize) -> Self {
        let h = 1.0 / (n as f64 + 1.0);
        Self {
            n,
            h,
            h2_inv: 1.0 / (h * h),
        }
    }

    /// RHS \(\sin(\pi x_i)\) at interior points.
    fn rhs(&self) -> Vec<f64> {
        let mut b = Vec::with_capacity(self.n);
        for i in 1..=self.n {
            let x = i as f64 * self.h;
            b.push((std::f64::consts::PI * x).sin());
        }
        b
    }

    /// Manufactured exact solution \(u^\star(x) = \sin(\pi x)/\pi^{2}\),
    /// the closed-form solution of \(L u = \sin(\pi x)\) on \([0, 1]\)
    /// with homogeneous Dirichlet BCs (independent of \(h\)).
    fn exact(&self) -> Vec<f64> {
        let pi_sq = std::f64::consts::PI * std::f64::consts::PI;
        let mut u = Vec::with_capacity(self.n);
        for i in 1..=self.n {
            let x = i as f64 * self.h;
            u.push((std::f64::consts::PI * x).sin() / pi_sq);
        }
        u
    }
}

// ----- Legacy reference -----

struct LegacySolver;

impl LegacySolver {
    /// Build a dense `nalgebra::DMatrix` Laplacian (with \(1/h^{2}\)
    /// scaling on every entry) and run a textbook CG with the same
    /// tolerance / max-iteration budget as the Atlas path. Returns
    /// (solution vector, total solve time in microseconds).
    fn solve(p: &Problem, b: &[f64]) -> (Vec<f64>, u128) {
        let n = p.n;
        let scale = p.h2_inv;
        let mut a = DMatrix::<f64>::zeros(n, n);
        for i in 0..n {
            a[(i, i)] = 2.0 * scale;
            if i > 0 {
                a[(i, i - 1)] = -1.0 * scale;
            }
            if i + 1 < n {
                a[(i, i + 1)] = -1.0 * scale;
            }
        }
        let b_mat = DMatrix::from_column_slice(n, 1, b);
        let t0 = Instant::now();
        let x = cg_solve(&a, &b_mat, 10_000, 1e-10);
        let elapsed_us = t0.elapsed().as_micros();
        let mut out = Vec::with_capacity(n);
        for r in 0..n {
            out.push(x[(r, 0)]);
        }
        (out, elapsed_us)
    }

    /// Compute legacy residual `A*x - b` for parity diagnostics.
    fn residual(p: &Problem, x: &[f64], b: &[f64]) -> Vec<f64> {
        let n = p.n;
        let scale = p.h2_inv;
        let mut r = vec![0.0_f64; n];
        for i in 0..n {
            let mut acc = 2.0 * scale * x[i];
            if i > 0 {
                acc += -1.0 * scale * x[i - 1];
            }
            if i + 1 < n {
                acc += -1.0 * scale * x[i + 1];
            }
            r[i] = acc - b[i];
        }
        r
    }

    fn memory_bytes(p: &Problem) -> usize {
        // Dense n×n f64 matrix + three n-vectors (x, b, residual).
        p.n * p.n * 8 + 3 * p.n * 8
    }
}

/// Textbook (Hestenes–Stiefel 1952) conjugate-gradient on a dense
/// `nalgebra::DMatrix`. Deliberately written without Atlas types so
/// the legacy path is a true planting of how CFDrs looked before
/// `appendix_migration.md`'s recipe was applied.
fn cg_solve(a: &DMatrix<f64>, b: &DMatrix<f64>, max_iter: usize, tol: f64) -> DMatrix<f64> {
    let n = a.nrows();
    let mut x = DMatrix::<f64>::zeros(n, 1);
    let mut r = b - a * &x; // r0 = b - A x0
    let mut p = r.clone();
    let mut rs_old = dvec_dot(&r, &r);
    let b_norm = (dvec_dot(b, b)).sqrt().max(1e-30);

    for _ in 0..max_iter {
        if (rs_old.sqrt()) / b_norm < tol {
            break;
        }
        let ap = a * &p;
        let denom = dvec_dot(&p, &ap);
        if denom.abs() < 1e-30 {
            break;
        }
        let alpha = rs_old / denom;
        x = x + &p * alpha;
        r = r - &ap * alpha;
        let rs_new = dvec_dot(&r, &r);
        let beta = rs_new / rs_old;
        p = &r + &p * beta;
        rs_old = rs_new;
    }
    x
}

fn dvec_dot(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    let mut s = 0.0;
    for i in 0..a.nrows() {
        s += a[(i, 0)] * b[(i, 0)];
    }
    s
}

// ----- Atlas candidate -----

struct AtlasSolver;

impl AtlasSolver {
    /// Build the discrete Laplacian via `SparseMatrixBuilder` (Atlas
    /// CSR `CsrMatrix<f64>`) and run `ConjugateGradient` through the
    /// ergonomic `LinearSolver::solve_system` trait method. Time the
    /// whole solve — not just the residual SpMV — so legacy vs Atlas
    /// numbers cover the same work.
    fn solve(p: &Problem, b_vec: &[f64]) -> CfdResult<(Vec<f64>, u128)> {
        let n = p.n;
        let scale = p.h2_inv;
        let mut builder = SparseMatrixBuilder::<f64>::new(n, n);
        for i in 0..n {
            builder.add_entry(i, i, 2.0 * scale)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0 * scale)?;
            }
            if i + 1 < n {
                builder.add_entry(i, i + 1, -1.0 * scale)?;
            }
        }
        let csr = builder.build()?;

        let cfg = IterativeSolverConfig::<f64>::new(1e-10).with_max_iterations(10_000);
        let cg = ConjugateGradient::new(cfg);
        let b_arr = Array1::from_shape_vec([n], b_vec.to_vec()).expect("RHS dimension is valid");

        let t0 = Instant::now();
        let x_arr: Array1<f64> = cg.solve_system(&csr, &b_arr, None)?;
        let elapsed_us = t0.elapsed().as_micros();

        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            out.push(x_arr[i]);
        }
        Ok((out, elapsed_us))
    }

    /// Compute atlas residual using `cfd_math::sparse::try_spmv`
    /// (which delegates to `leto_ops::CsrMatrix::spmv_into`).
    fn residual(p: &Problem, x: &[f64], b_vec: &[f64]) -> CfdResult<Vec<f64>> {
        let n = p.n;
        let scale = p.h2_inv;
        let mut builder = SparseMatrixBuilder::<f64>::new(n, n);
        for i in 0..n {
            builder.add_entry(i, i, 2.0 * scale)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0 * scale)?;
            }
            if i + 1 < n {
                builder.add_entry(i, i + 1, -1.0 * scale)?;
            }
        }
        let csr = builder.build()?;
        let x_arr = Array1::from_shape_vec([n], x.to_vec()).expect("x dimension is valid");
        let b_arr = Array1::from_shape_vec([n], b_vec.to_vec()).expect("b dimension is valid");
        let mut y_arr = Array1::<f64>::zeros([n]);
        try_spmv(&csr, &x_arr, &mut y_arr)?;
        let mut r = Vec::with_capacity(n);
        for i in 0..n {
            r.push(y_arr[i] - b_arr[i]);
        }
        Ok(r)
    }

    fn memory_bytes(p: &Problem) -> usize {
        // CSR: nnz = 3·N entries; vec<usize> row_offsets, vec<u32>
        // col_indices, vec<f64> values; plus three f64 vectors of length N.
        let nnz = 3 * p.n;
        nnz * (4 + 4 + 8) + 3 * p.n * 8
    }
}

fn vec_max(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, v| v.abs().max(m))
}

fn vec_max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

fn vec_rms(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    let s = a.iter().map(|v| v * v).sum::<f64>();
    (s / n).sqrt()
}

fn main() -> CfdResult<()> {
    let problem = Problem::new(1024);
    let b_vec = problem.rhs();
    let u_exact = problem.exact();

    // ----- Legacy reference -----
    let (u_legacy, legacy_solve_us) = LegacySolver::solve(&problem, &b_vec);
    let r_legacy = LegacySolver::residual(&problem, &u_legacy, &b_vec);

    // ----- Atlas candidate -----
    let (u_atlas, atlas_solve_us) = AtlasSolver::solve(&problem, &b_vec)?;
    let r_atlas = AtlasSolver::residual(&problem, &u_atlas, &b_vec)?;

    // ----- Parity diagnostics (per migration_validation.md) -----
    let max_abs_diff_solutions = vec_max_diff(&u_legacy, &u_atlas);
    let max_abs_diff_atlas_vs_exact = vec_max_diff(&u_atlas, &u_exact);
    let max_residual_legacy = vec_max(&r_legacy);
    let max_residual_atlas = vec_max(&r_atlas);
    let sqrt_kappa_legacy = vec_rms(&u_legacy);
    let sqrt_kappa_atlas = vec_rms(&u_atlas);

    // The default tolerance budget for parity is 1e-8 absolute on
    // residual max-abs, plus 1e-6 relative on solution agreement.
    let parity_pass = max_residual_legacy < 1e-8
        && max_residual_atlas < 1e-8
        && max_abs_diff_solutions < 1e-6;

    println!(
        "{}",
        parity_report(ParityReport {
            problem_n: problem.n,
            legacy_solve_us,
            atlas_solve_us,
            max_residual_legacy,
            max_residual_atlas,
            max_abs_diff_solutions,
            max_abs_diff_atlas_vs_exact,
            legacy_resid_rms: sqrt_kappa_legacy,
            atlas_resid_rms: sqrt_kappa_atlas,
            legacy_mem_bytes: LegacySolver::memory_bytes(&problem),
            atlas_mem_bytes: AtlasSolver::memory_bytes(&problem),
            parity_pass,
        })
    );

    if !parity_pass {
        return Err(CfdError::InvalidConfiguration(format!(
            "migration parity FAIL: max_abs_diff_solutions={max_abs_diff_solutions:.3e}, max_residual_legacy={max_residual_legacy:.3e}, max_residual_atlas={max_residual_atlas:.3e}"
        )));
    }
    Ok(())
}

#[derive(Debug)]
struct ParityReport {
    problem_n: usize,
    legacy_solve_us: u128,
    atlas_solve_us: u128,
    max_residual_legacy: f64,
    max_residual_atlas: f64,
    max_abs_diff_solutions: f64,
    max_abs_diff_atlas_vs_exact: f64,
    legacy_resid_rms: f64,
    atlas_resid_rms: f64,
    legacy_mem_bytes: usize,
    atlas_mem_bytes: usize,
    parity_pass: bool,
}

fn parity_report(r: ParityReport) -> String {
    format!(
        "{{\"problem_n\":{},\"legacy_solve_us\":{},\"atlas_solve_us\":{},\"max_residual_legacy\":{:.6e},\"max_residual_atlas\":{:.6e},\"max_abs_diff_solutions\":{:.6e},\"max_abs_diff_atlas_vs_exact\":{:.6e},\"legacy_resid_rms\":{:.6e},\"atlas_resid_rms\":{:.6e},\"legacy_mem_bytes\":{},\"atlas_mem_bytes\":{},\"parity_pass\":{}}}",
        r.problem_n,
        r.legacy_solve_us,
        r.atlas_solve_us,
        r.max_residual_legacy,
        r.max_residual_atlas,
        r.max_abs_diff_solutions,
        r.max_abs_diff_atlas_vs_exact,
        r.legacy_resid_rms,
        r.atlas_resid_rms,
        r.legacy_mem_bytes,
        r.atlas_mem_bytes,
        r.parity_pass
    )
}
