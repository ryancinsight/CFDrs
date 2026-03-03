//! Root-level integration tests that validate CFDrs solvers against known
//! analytical solutions.
//!
//! These replicate the same problems used in the `validation_vs_python_cfd`
//! and `validation_vs_pmocz` examples, but as `#[test]` functions with hard
//! error tolerances suitable for CI.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Test 1: Poiseuille Channel Flow
// ---------------------------------------------------------------------------

use cfd_2d::solvers::poiseuille::{BloodModel as PoisBloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
use cfd_core::physics::fluid::blood::CassonBlood;

#[test]
fn poiseuille_2d_analytical_match() {
    // Zero-yield Casson = Newtonian, mu = 0.001 Pa-s
    let blood = CassonBlood::new(1000.0, 0.0, 0.001, 0.0);
    let config = PoiseuilleConfig {
        height: 100e-6,
        width: 500e-6,
        length: 1e-3,
        ny: 101,
        pressure_gradient: 1000.0,
        tolerance: 1e-6,
        max_iterations: 1000,
        relaxation_factor: 0.7,
    };
    let mut solver = PoiseuilleFlow2D::new(config, PoisBloodModel::Casson(blood));
    solver.solve().unwrap();

    let u_max = solver.max_velocity();
    // Analytical: u_max = H^2 * |dP/dx| / (8 * mu)
    let u_max_analytical = (100e-6_f64).powi(2) * 1000.0 / (8.0 * 0.001);
    let error = ((u_max - u_max_analytical) / u_max_analytical).abs();
    assert!(
        error < 0.01,
        "Poiseuille u_max error {:.2e} exceeds 1%",
        error
    );
}

// ---------------------------------------------------------------------------
// Test 2: Spectral Poisson Convergence
// ---------------------------------------------------------------------------

use cfd_3d::spectral::solver::PoissonProblem;
use cfd_3d::spectral::{PoissonBoundaryCondition, SpectralConfig, SpectralSolver};
use nalgebra::DVector;

#[test]
fn spectral_poisson_l2_error_below_threshold() {
    let n = 8;
    let config = SpectralConfig::<f64>::new(n, n, n).unwrap();
    let mut solver = SpectralSolver::new(config).unwrap();

    // Build source for -Laplacian(u) = 3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)
    // on the Chebyshev Gauss-Lobatto domain [-1, 1]^3
    let mut source = DVector::zeros(n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);
                // Flat index: i * ny * nz + j * nz + k  (row-major, x outermost)
                source[i * n * n + j * n + k] =
                    3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
            }
        }
    }

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

    let solution = solver.solve(&problem).unwrap();

    // Compute L2 error against exact solution sin(pi*x)*sin(pi*y)*sin(pi*z)
    let mut l2_sum = 0.0;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);
                let exact = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                let err = solution.u[i * n * n + j * n + k] - exact;
                l2_sum += err * err;
            }
        }
    }
    let l2 = (l2_sum / (n * n * n) as f64).sqrt();
    assert!(
        l2 < 0.5,
        "Spectral Poisson L2 error = {}, expected < 0.5",
        l2
    );
}

/// Chebyshev Gauss-Lobatto collocation point on [-1, 1].
///
/// x_i = cos(i * pi / (n - 1)),  i = 0 .. n-1
fn chebyshev_node(i: usize, n: usize) -> f64 {
    let theta = PI * i as f64 / (n - 1).max(1) as f64;
    theta.cos()
}

// ---------------------------------------------------------------------------
// Test 3: Ghia Reference Data Available
// ---------------------------------------------------------------------------

use cfd_validation::benchmarks::LidDrivenCavity;

#[test]
fn ghia_reference_data_loads() {
    // Verify the validation crate can load Ghia Re=100 data
    let cavity = LidDrivenCavity::new(1.0, 1.0, 100.0);
    let u_data = cavity.ghia_u_centerline(100.0);
    // Check data has expected number of points (Ghia 1982 Table I has 17 points)
    assert!(
        !u_data.is_empty(),
        "Ghia u-centerline data must be non-empty"
    );
    assert_eq!(
        u_data.len(),
        17,
        "Ghia 1982 Table I has 17 tabulated points"
    );
}
