//! Phase 7 — Property-based tests for the cfd-3d crate.
//!
//! Uses `proptest` to verify invariants across randomised parameter ranges
//! and deterministic tests for spectral Poisson solver correctness.

use proptest::prelude::*;

// ---------------------------------------------------------------------------
// Property test: SpectralConfig::new() with valid sizes always succeeds
// ---------------------------------------------------------------------------

proptest! {
    #[test]
    fn spectral_config_valid_sizes(n in 2usize..=16) {
        use cfd_3d::SpectralConfig;
        let result = SpectralConfig::<f64>::new(n, n, n);
        prop_assert!(result.is_ok(), "SpectralConfig::new({n},{n},{n}) must succeed");
    }
}

// ---------------------------------------------------------------------------
// Property test: SpectralSolver construction succeeds for valid configs
// ---------------------------------------------------------------------------

proptest! {
    #[test]
    fn spectral_solver_creation_valid(n in 2usize..=8) {
        use cfd_3d::{SpectralConfig, SpectralSolver};
        let config = SpectralConfig::<f64>::new(n, n, n)
            .expect("SpectralConfig::new must succeed");
        let result = SpectralSolver::new(config);
        prop_assert!(result.is_ok(), "SpectralSolver::new must succeed for n={n}");
    }
}

// ---------------------------------------------------------------------------
// Deterministic test: Spectral Poisson solver produces finite solutions
// ---------------------------------------------------------------------------

#[test]
fn spectral_poisson_error_is_finite() {
    use cfd_3d::spectral::solver::{PoissonProblem, SpectralConfig, SpectralSolver};
    use cfd_3d::spectral::PoissonBoundaryCondition;
    use nalgebra::DVector;

    let n = 4;
    let config = SpectralConfig::<f64>::new(n, n, n).unwrap();
    let mut solver = SpectralSolver::new(config).unwrap();
    let source = DVector::from_element(n * n * n, 1.0);
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
    for val in solution.u.iter() {
        assert!(val.is_finite(), "Solution must be finite, got {val}");
    }
}

// ---------------------------------------------------------------------------
// Deterministic test: Poisson solution has correct dimensions
// ---------------------------------------------------------------------------

#[test]
fn spectral_poisson_solution_dimensions() {
    use cfd_3d::spectral::solver::{PoissonProblem, SpectralConfig, SpectralSolver};
    use cfd_3d::spectral::PoissonBoundaryCondition;
    use nalgebra::DVector;

    let nx = 4;
    let ny = 5;
    let nz = 3;
    let config = SpectralConfig::<f64>::new(nx, ny, nz).unwrap();
    let mut solver = SpectralSolver::new(config).unwrap();
    let source = DVector::from_element(nx * ny * nz, 0.0);
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

    assert_eq!(solution.nx, nx, "Solution nx must match config");
    assert_eq!(solution.ny, ny, "Solution ny must match config");
    assert_eq!(solution.nz, nz, "Solution nz must match config");
    assert_eq!(
        solution.u.len(),
        nx * ny * nz,
        "Solution vector length must equal nx*ny*nz"
    );
}
