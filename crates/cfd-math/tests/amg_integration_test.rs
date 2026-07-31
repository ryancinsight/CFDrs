//! Integration tests for the AMG preconditioner with BiCGSTAB and GMRES.
//!
//! Oracle: the 2-D five-point Poisson operator with a manufactured solution
//! `x*` and `b = A·x*`. Krylov solves preconditioned by
//! [`AlgebraicMultigrid`] must recover `x*` to the solver tolerance, and the
//! stationary preconditioned iteration must contract in the A-norm with
//! factor `ρ < 1` (Ruge & Stüben 1987; theorem restated in `amg.rs`).

use cfd_math::iterative::{IterativeSolverConfig, Preconditioner, BiCGSTAB, GMRES};
use cfd_math::multigrid::{AMGConfig, AlgebraicMultigrid, CycleType};
use leto::Array1;
use leto_ops::{spmv_into, CsrMatrix};

/// Assemble the 2-D five-point Poisson operator on an `n × n` grid
/// (order `n²`): `4` on the diagonal, `-1` toward each grid neighbour.
fn poisson_2d(n: usize) -> CsrMatrix<f64> {
    let size = n * n;
    let mut row_offsets = Vec::with_capacity(size + 1);
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
    row_offsets.push(0);
    for i in 0..size {
        let row = i / n;
        let col = i % n;
        if row > 0 {
            col_indices.push(i - n);
            values.push(-1.0);
        }
        if col > 0 {
            col_indices.push(i - 1);
            values.push(-1.0);
        }
        col_indices.push(i);
        values.push(4.0);
        if col + 1 < n {
            col_indices.push(i + 1);
            values.push(-1.0);
        }
        if row + 1 < n {
            col_indices.push(i + n);
            values.push(-1.0);
        }
        row_offsets.push(col_indices.len());
    }
    CsrMatrix::from_parts(values, col_indices, row_offsets, size, size)
        .expect("five-point stencil assembly is structurally valid")
}

/// Manufactured smooth solution sampled on the grid.
fn manufactured_solution(size: usize) -> Array1<f64> {
    let values = (0..size)
        .map(|k| {
            let t = k as f64 / size as f64;
            (std::f64::consts::PI * t).sin() + 0.5 * t
        })
        .collect::<Vec<_>>();
    Array1::from_shape_vec([size], values).expect("solution length matches size")
}

/// `b = A · x*` so the solve has an exact reference.
fn rhs_for(a: &CsrMatrix<f64>, x_exact: &Array1<f64>) -> Array1<f64> {
    let mut b = Array1::zeros([x_exact.shape()[0]]);
    spmv_into(a, &x_exact.view(), b.as_slice_mut().expect("contiguous rhs"))
        .expect("SpMV over matching shapes");
    b
}

fn max_abs_diff(lhs: &Array1<f64>, rhs: &Array1<f64>) -> f64 {
    assert_eq!(lhs.shape(), rhs.shape(), "compared arrays share shape");
    (0..lhs.shape()[0])
        .map(|idx| (lhs[idx] - rhs[idx]).abs())
        .fold(0.0, f64::max)
}

/// `‖v‖_A = sqrt(vᵀ A v)` — the energy norm the AMG theory contracts in.
fn energy_norm(a: &CsrMatrix<f64>, v: &Array1<f64>) -> f64 {
    let mut av = Array1::zeros([v.shape()[0]]);
    spmv_into(a, &v.view(), av.as_slice_mut().expect("contiguous workspace"))
        .expect("SpMV over matching shapes");
    let dot = (0..v.shape()[0]).map(|idx| v[idx] * av[idx]).sum::<f64>();
    dot.max(0.0).sqrt()
}

fn amg_for(a: &CsrMatrix<f64>, config: AMGConfig) -> AlgebraicMultigrid<f64> {
    AlgebraicMultigrid::new(a, config).expect("AMG setup on SPD Poisson operator")
}

#[test]
fn test_amg_with_bicgstab() {
    let n = 16;
    let a = poisson_2d(n);
    let x_exact = manufactured_solution(n * n);
    let b = rhs_for(&a, &x_exact);
    let amg = amg_for(&a, AMGConfig::default());

    let solver = BiCGSTAB::new(IterativeSolverConfig::new(1e-10).with_max_iterations(200));
    let mut x = Array1::zeros([n * n]);
    let monitor = solver
        .solve_preconditioned(&a, &b, &amg, &mut x)
        .expect("AMG-preconditioned BiCGSTAB converges on SPD Poisson");

    let final_residual = *monitor
        .residual_history
        .last()
        .expect("monitor records at least the initial residual");
    assert!(
        final_residual <= 1e-10 * monitor.initial_residual.max(1.0),
        "relative residual not reached: final {final_residual}, initial {}",
        monitor.initial_residual
    );
    let error = max_abs_diff(&x, &x_exact);
    assert!(error < 1e-6, "‖x - x*‖∞ = {error}");
}

#[test]
fn test_amg_with_gmres() {
    let n = 16;
    let a = poisson_2d(n);
    let x_exact = manufactured_solution(n * n);
    let b = rhs_for(&a, &x_exact);
    let amg = amg_for(&a, AMGConfig::default());

    let solver = GMRES::new(
        IterativeSolverConfig::new(1e-10).with_max_iterations(200),
        30,
    );
    let mut x = Array1::zeros([n * n]);
    solver
        .solve_preconditioned(&a, &b, &amg, &mut x)
        .expect("AMG-preconditioned GMRES converges on SPD Poisson");

    let error = max_abs_diff(&x, &x_exact);
    assert!(error < 1e-6, "‖x - x*‖∞ = {error}");
}

#[test]
fn test_amg_different_cycles() {
    let n = 16;
    let a = poisson_2d(n);
    let x_exact = manufactured_solution(n * n);
    let b = rhs_for(&a, &x_exact);

    for cycle in [CycleType::VCycle, CycleType::WCycle, CycleType::FCycle] {
        let amg = amg_for(
            &a,
            AMGConfig {
                cycle_type: cycle,
                ..AMGConfig::default()
            },
        );
        let solver = BiCGSTAB::new(IterativeSolverConfig::new(1e-10).with_max_iterations(200));
        let mut x = Array1::zeros([n * n]);
        solver
            .solve_preconditioned(&a, &b, &amg, &mut x)
            .unwrap_or_else(|e| panic!("{cycle:?} preconditioned solve failed: {e}"));
        let error = max_abs_diff(&x, &x_exact);
        assert!(error < 1e-6, "{cycle:?}: ‖x - x*‖∞ = {error}");
    }
}

#[test]
fn test_amg_construction_edge_cases() {
    // A matrix already at or below min_coarse_size builds a single-level
    // hierarchy whose application is still a usable preconditioner.
    let n = 4; // order 16 < default min_coarse_size (50)
    let a = poisson_2d(n);
    let amg = amg_for(&a, AMGConfig::default());
    let r = manufactured_solution(n * n);
    let mut z = Array1::zeros([n * n]);
    amg.apply_to(&r, &mut z)
        .expect("single-level AMG application");
    // A preconditioner approximating A⁻¹ must not annihilate a nonzero
    // residual: z = M⁻¹ r ≠ 0.
    assert!(
        (0..z.shape()[0]).any(|idx| z[idx] != 0.0),
        "preconditioner returned the zero vector for a nonzero residual"
    );

    // Dimension mismatch at application is a typed error, not a panic.
    let mut wrong = Array1::zeros([n * n + 1]);
    assert!(
        amg.apply_to(&r, &mut wrong).is_err(),
        "output-length mismatch must be rejected"
    );

    // A capped hierarchy depth (max_levels = 1) is honored.
    let deep = poisson_2d(16);
    let amg_shallow = amg_for(
        &deep,
        AMGConfig {
            max_levels: 1,
            ..AMGConfig::default()
        },
    );
    let r = manufactured_solution(256);
    let mut z = Array1::zeros([256]);
    amg_shallow
        .apply_to(&r, &mut z)
        .expect("depth-1 AMG application");
}

#[test]
fn test_amg_two_grid_convergence_factor() {
    // Stationary preconditioned Richardson iteration
    //   x_{k+1} = x_k + M⁻¹ (b − A x_k)
    // has error propagation e_{k+1} = (I − M⁻¹A) e_k. Ruge–Stüben theory
    // bounds the A-norm contraction factor ρ < 1 independently of h for
    // SPD M-matrices; assert the contraction on two grid sizes.
    for n in [8usize, 16] {
        let size = n * n;
        let a = poisson_2d(n);
        let x_exact = manufactured_solution(size);
        let b = rhs_for(&a, &x_exact);
        let amg = amg_for(&a, AMGConfig::default());

        let mut x = Array1::zeros([size]);
        let mut residual = Array1::zeros([size]);
        let mut z = Array1::zeros([size]);

        let mut error = Array1::zeros([size]);
        for idx in 0..size {
            error[idx] = x_exact[idx] - x[idx];
        }
        let mut previous_energy = energy_norm(&a, &error);
        assert!(previous_energy > 0.0, "nonzero initial error required");

        let mut worst_factor = 0.0_f64;
        for _ in 0..5 {
            spmv_into(
                &a,
                &x.view(),
                residual.as_slice_mut().expect("contiguous residual"),
            )
            .expect("SpMV over matching shapes");
            for idx in 0..size {
                residual[idx] = b[idx] - residual[idx];
            }
            amg.apply_to(&residual, &mut z).expect("AMG application");
            for idx in 0..size {
                x[idx] += z[idx];
            }

            for idx in 0..size {
                error[idx] = x_exact[idx] - x[idx];
            }
            let energy = energy_norm(&a, &error);
            if previous_energy > 0.0 {
                worst_factor = worst_factor.max(energy / previous_energy);
            }
            previous_energy = energy;
            if energy == 0.0 {
                break;
            }
        }
        assert!(
            worst_factor < 1.0,
            "n = {n}: A-norm contraction factor {worst_factor} is not < 1"
        );
    }
}
