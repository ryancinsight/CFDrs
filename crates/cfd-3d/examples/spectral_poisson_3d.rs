//! 3D Spectral Poisson Solver
//!
//! Demonstrates the spectral solver on a 3D Poisson problem with a known
//! analytical solution, verifying exponential convergence in the number of modes.
//!
//! # Problem
//!
//! Solve: -nabla^2 u = f  on [0,1]^3
//! with Dirichlet BCs: u = 0 on all boundaries
//!
//! Source term: f = 3 * pi^2 * sin(pi*x) * sin(pi*y) * sin(pi*z)
//! Exact solution: u = sin(pi*x) * sin(pi*y) * sin(pi*z)
//!
//! # Reference
//!
//! Canuto, C. et al. (2006). "Spectral Methods: Fundamentals in Single Domains".
//! Springer. Chapter 3: Spectral approximation on the interval.

use cfd_3d::spectral::{SpectralConfig, SpectralSolver};
use cfd_3d::spectral::poisson::PoissonBoundaryCondition;
use cfd_3d::spectral::solver::PoissonProblem;
use nalgebra::DVector;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D Spectral Poisson Solver — Convergence Study");
    println!("================================================\n");

    println!("Problem:");
    println!("  -laplacian(u) = 3*pi^2 * sin(pi*x)*sin(pi*y)*sin(pi*z)");
    println!("  Domain: [0,1]^3");
    println!("  BCs: Dirichlet u=0 on all boundaries");
    println!("  Exact: u = sin(pi*x)*sin(pi*y)*sin(pi*z)");

    // Convergence study over increasing mode counts
    let mode_counts = [4, 6, 8, 10, 12];

    println!("\nConvergence Study:");
    println!("  {:>6} {:>12} {:>12} {:>12}", "N", "L2 Error", "Linf Error", "DOFs");
    println!("  {}", "-".repeat(48));

    let mut prev_l2 = f64::MAX;

    for &n in &mode_counts {
        let config = SpectralConfig::<f64>::new(n, n, n)?;
        let mut solver = SpectralSolver::new(config)?;

        // Build source term on Chebyshev-Gauss-Lobatto nodes
        let source = build_source_term(n);
        let problem = PoissonProblem {
            source_term: source,
            bc_x: (
                PoissonBoundaryCondition::Dirichlet(0.0),
                PoissonBoundaryCondition::Dirichlet(0.0),
            ),
            bc_y: (
                PoissonBoundaryCondition::Dirichlet(0.0),
                PoissonBoundaryCondition::Dirichlet(0.0),
            ),
            bc_z: (
                PoissonBoundaryCondition::Dirichlet(0.0),
                PoissonBoundaryCondition::Dirichlet(0.0),
            ),
        };

        let solution = solver.solve(&problem)?;

        // Compute error against exact solution
        let (l2_err, linf_err) = compute_errors(&solution.u, n);
        let dofs = n * n * n;

        println!("  {:>6} {:>12.4e} {:>12.4e} {:>12}", n, l2_err, linf_err, dofs);

        // Verify convergence (error should decrease)
        if n > mode_counts[0] {
            assert!(
                l2_err < prev_l2,
                "L2 error should decrease with more modes: {} >= {}",
                l2_err,
                prev_l2
            );
        }
        prev_l2 = l2_err;
    }

    println!("\nAll convergence checks PASSED.");
    println!("Spectral 3D Poisson example completed.");
    Ok(())
}

/// Build the source term f = 3*pi^2 * sin(pi*x)*sin(pi*y)*sin(pi*z)
/// evaluated on Chebyshev-Gauss-Lobatto nodes mapped to [0,1].
fn build_source_term(n: usize) -> DVector<f64> {
    let mut source = DVector::zeros(n * n * n);

    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);

                let idx = k * n * n + j * n + i;
                source[idx] = 3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
            }
        }
    }

    source
}

/// Chebyshev-Gauss-Lobatto node mapped from [-1,1] to [0,1].
fn chebyshev_node(i: usize, n: usize) -> f64 {
    let theta = PI * i as f64 / (n - 1).max(1) as f64;
    0.5 * (1.0 - theta.cos())
}

/// Compute L2 and L-infinity errors against the exact solution.
fn compute_errors(u_numerical: &DVector<f64>, n: usize) -> (f64, f64) {
    let mut l2_sum = 0.0;
    let mut linf = 0.0_f64;
    let total = n * n * n;

    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);

                let exact = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                let idx = k * n * n + j * n + i;
                let err = (u_numerical[idx] - exact).abs();

                l2_sum += err * err;
                linf = linf.max(err);
            }
        }
    }

    let l2 = (l2_sum / total as f64).sqrt();
    (l2, linf)
}
