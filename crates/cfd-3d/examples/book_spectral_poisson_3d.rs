#![allow(missing_docs)]
//! Book example: 3D spectral Poisson solve with Apollo-backed FFT operators.

use cfd_3d::spectral::poisson::PoissonBoundaryCondition;
use cfd_3d::spectral::solver::PoissonProblem;
use cfd_3d::spectral::{SpectralConfig, SpectralSolver};
use leto::Array1;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = 8usize;
    let config = SpectralConfig::<f64>::new(n, n, n)?;
    let mut solver = SpectralSolver::new(config)?;

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
    let l2_error = l2_error_against_exact(&solution.u, n);
    assert!(l2_error.is_finite(), "L2 error must be finite");
    assert!(l2_error < 0.4, "L2 error expected < 0.4, got {l2_error}");

    Ok(())
}

fn build_source_term(n: usize) -> Array1<f64> {
    let mut source = Array1::zeros([n * n * n]);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);
                let idx = k * n * n + j * n + i;
                source[[idx]] = 3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
            }
        }
    }
    source
}

fn l2_error_against_exact(u_numerical: &Array1<f64>, n: usize) -> f64 {
    let mut sum = 0.0f64;
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);
                let exact = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                let idx = k * n * n + j * n + i;
                let err = u_numerical[[idx]] - exact;
                sum += err * err;
            }
        }
    }
    (sum / (n * n * n) as f64).sqrt()
}

fn chebyshev_node(i: usize, n: usize) -> f64 {
    let theta = PI * i as f64 / (n - 1).max(1) as f64;
    0.5 * (1.0 - theta.cos())
}
