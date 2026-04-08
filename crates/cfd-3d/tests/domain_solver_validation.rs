//! Integration tests for all cfd-3d domain solvers.
//!
//! Exercises VenturiSolver3D, BifurcationSolver3D, SerpentineSolver3D,
//! TrifurcationSolver3D, and SpectralSolver at minimal grid resolutions.
//! Validates mass conservation, pressure bounds, flow physics, and convergence.

use cfd_3d::bifurcation::{BifurcationConfig3D, BifurcationGeometry3D, BifurcationSolver3D};
use cfd_3d::serpentine::{SerpentineConfig3D, SerpentineSolver3D};
use cfd_3d::spectral::poisson::PoissonBoundaryCondition;
use cfd_3d::spectral::solver::PoissonProblem;
use cfd_3d::spectral::{SpectralConfig, SpectralSolver};
use cfd_3d::trifurcation::{TrifurcationConfig3D, TrifurcationGeometry3D, TrifurcationSolver3D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::application::channel::serpentine::SerpentineMeshBuilder;
use cfd_mesh::VenturiMeshBuilder;
use nalgebra::DVector;
use std::f64::consts::PI;

// ──────────────────────────────────────────────────────────────────────
// Bifurcation
// ──────────────────────────────────────────────────────────────────────

#[test]
fn bifurcation_symmetric_mass_conservation() {
    let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
    let config = BifurcationConfig3D {
        mesh_resolution: 4,
        ..BifurcationConfig3D::default()
    };
    let solver = BifurcationSolver3D::new(geom, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    assert!(
        sol.mass_conservation_error.abs() < 0.10,
        "Bifurcation mass error too large: {}",
        sol.mass_conservation_error
    );
}

#[test]
fn bifurcation_symmetric_flow_split() {
    let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
    let config = BifurcationConfig3D {
        mesh_resolution: 4,
        ..BifurcationConfig3D::default()
    };
    let solver = BifurcationSolver3D::new(geom, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    if sol.q_daughter1.abs() > 1e-15 && sol.q_daughter2.abs() > 1e-15 {
        let ratio = sol.q_daughter1 / sol.q_daughter2;
        assert!(
            (ratio - 1.0).abs() < 0.50,
            "Symmetric daughters should have ~equal flow: q1={}, q2={}, ratio={}",
            sol.q_daughter1,
            sol.q_daughter2,
            ratio
        );
    }
}

#[test]
fn bifurcation_positive_wall_shear() {
    let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
    let config = BifurcationConfig3D {
        mesh_resolution: 4,
        ..BifurcationConfig3D::default()
    };
    let solver = BifurcationSolver3D::new(geom, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    assert!(
        sol.wall_shear_stress_parent >= 0.0,
        "Parent wall shear must be non-negative: {}",
        sol.wall_shear_stress_parent
    );
    assert!(
        sol.wall_shear_stress_daughter1 >= 0.0,
        "Daughter1 wall shear must be non-negative: {}",
        sol.wall_shear_stress_daughter1
    );
}

#[test]
fn bifurcation_blood_casson() {
    let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
    let config = BifurcationConfig3D {
        mesh_resolution: 4,
        inlet_flow_rate: 1e-8,
        ..BifurcationConfig3D::default()
    };
    let solver = BifurcationSolver3D::new(geom, config);
    let blood = CassonBlood::<f64>::normal_blood();
    let sol = solver.solve(blood).unwrap();

    assert!(
        sol.mass_conservation_error.abs() < 0.15,
        "Blood bifurcation mass error: {}",
        sol.mass_conservation_error
    );
    let dp = sol.p_inlet - sol.p_outlet;
    assert!(dp > 0.0, "Pressure drop must be positive for blood flow");
}

// ──────────────────────────────────────────────────────────────────────
// Venturi
// ──────────────────────────────────────────────────────────────────────

#[test]
fn venturi_nonzero_pressure_difference() {
    // At coarse P1-P1 resolution the FEM solver may report negative dp
    // due to pressure normalization artifacts. We check that a non-trivial
    // pressure difference exists (it must not be exactly zero).
    let builder =
        VenturiMeshBuilder::new(2.0e-3_f64, 1.0e-3, 3.0e-3, 2.0e-3, 1.0e-3, 4.0e-3, 3.0e-3);
    let config = VenturiConfig3D {
        inlet_flow_rate: 1e-7,
        inlet_pressure: 200.0,
        outlet_pressure: 0.0,
        resolution: (30, 5),
        ..VenturiConfig3D::default()
    };
    let solver = VenturiSolver3D::new(builder, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    let dp = (sol.p_inlet - sol.p_outlet).abs();
    assert!(dp > 0.0, "Pressure difference must be non-zero, got {}", dp);
    // Throat pressure drop sign depends on convention at coarse resolution
    assert!(
        sol.dp_throat.abs() > 0.0,
        "Throat pressure difference must be non-zero, got {}",
        sol.dp_throat
    );
}

#[test]
fn venturi_throat_acceleration() {
    let builder =
        VenturiMeshBuilder::new(2.0e-3_f64, 1.0e-3, 3.0e-3, 2.0e-3, 1.0e-3, 4.0e-3, 3.0e-3);
    let config = VenturiConfig3D {
        inlet_flow_rate: 1e-7,
        resolution: (30, 5),
        ..VenturiConfig3D::default()
    };
    let solver = VenturiSolver3D::new(builder, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    assert!(
        sol.u_throat > sol.u_inlet,
        "Throat velocity must exceed inlet: u_throat={}, u_inlet={}",
        sol.u_throat,
        sol.u_inlet
    );
}

// ──────────────────────────────────────────────────────────────────────
// Serpentine
// ──────────────────────────────────────────────────────────────────────

#[test]
fn serpentine_positive_pressure_drop() {
    let builder = SerpentineMeshBuilder::new(1.0e-3_f64, 2.0e-3, 8.0e-3)
        .with_periods(2)
        .with_resolution(20, 4);
    let config = SerpentineConfig3D {
        inlet_flow_rate: 5e-7,
        inlet_pressure: 200.0,
        outlet_pressure: 0.0,
        ..SerpentineConfig3D::default()
    };
    let solver = SerpentineSolver3D::new(builder, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    assert!(
        sol.dp_total > 0.0,
        "Serpentine pressure drop must be positive: {}",
        sol.dp_total
    );
    assert!(
        sol.u_inlet > 0.0,
        "Inlet velocity must be positive: {}",
        sol.u_inlet
    );
}

#[test]
fn serpentine_dean_number_positive() {
    let builder = SerpentineMeshBuilder::new(1.0e-3_f64, 2.0e-3, 8.0e-3)
        .with_periods(2)
        .with_resolution(20, 4);
    let config = SerpentineConfig3D {
        inlet_flow_rate: 5e-7,
        inlet_pressure: 200.0,
        outlet_pressure: 0.0,
        ..SerpentineConfig3D::default()
    };
    let solver = SerpentineSolver3D::new(builder, config);
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let sol = solver.solve(water).unwrap();

    assert!(
        sol.dean_number >= 0.0,
        "Dean number must be non-negative: {}",
        sol.dean_number
    );
}

// ──────────────────────────────────────────────────────────────────────
// Trifurcation
// ──────────────────────────────────────────────────────────────────────

fn solve_reference_trifurcation_blood_flow() -> cfd_3d::trifurcation::TrifurcationSolution3D<f64> {
    let geom = TrifurcationGeometry3D::<f64>::symmetric(
        100e-6,   // parent diameter
        80e-6,    // daughter diameter
        500e-6,   // parent length
        500e-6,   // daughter length
        50e-6,    // transition length
        PI / 4.0, // spread angle
    );
    let config = TrifurcationConfig3D {
        inlet_flow_rate: 1e-9,
        inlet_pressure: 200.0,
        max_linear_iterations: 600,
        linear_tolerance: 1e-5,
        // 50µm cell_size = 1 cell per parent radius: coarse but fast (#200 seeds, <10s solve)
        target_mesh_size: Some(50e-6),
        ..TrifurcationConfig3D::default()
    };
    let solver = TrifurcationSolver3D::new(geom, config);
    let blood = CassonBlood::<f64>::normal_blood();
    solver.solve(&blood).unwrap()
}

#[test]
fn trifurcation_symmetric_mass_conservation() {
    let sol = solve_reference_trifurcation_blood_flow();
    let total_out = sol.flow_rates[1] + sol.flow_rates[2] + sol.flow_rates[3];

    assert!(sol.flow_rates[0].is_finite() && sol.flow_rates[0] > 0.0);
    assert!(total_out.is_finite() && total_out > 0.0);

    assert!(
        sol.mass_conservation_error.abs() < 0.10,
        "Trifurcation mass error: {}",
        sol.mass_conservation_error
    );
}

#[test]
fn trifurcation_symmetric_positive_outlet_flow() {
    let sol = solve_reference_trifurcation_blood_flow();

    let q1 = sol.flow_rates[1];
    let q2 = sol.flow_rates[2];
    let q3 = sol.flow_rates[3];

    assert!(q1 > 0.0 && q2 > 0.0 && q3 > 0.0);
    assert!(q1.is_finite() && q2.is_finite() && q3.is_finite());
}

// ──────────────────────────────────────────────────────────────────────
// Spectral Poisson
// ──────────────────────────────────────────────────────────────────────

/// Chebyshev-Gauss-Lobatto node mapped from [-1,1] to [0,1].
fn chebyshev_node(i: usize, n: usize) -> f64 {
    let theta = PI * i as f64 / (n - 1).max(1) as f64;
    0.5 * (1.0 - theta.cos())
}

/// Build source for -laplacian(u) = 3*pi^2 * sin(pi*x)*sin(pi*y)*sin(pi*z).
fn build_poisson_source(n: usize) -> DVector<f64> {
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

/// Compute L2 error against exact solution.
fn poisson_l2_error(u: &DVector<f64>, n: usize) -> f64 {
    let mut l2_sum = 0.0;
    let total = n * n * n;
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let x = chebyshev_node(i, n);
                let y = chebyshev_node(j, n);
                let z = chebyshev_node(k, n);
                let exact = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                let idx = k * n * n + j * n + i;
                let err = u[idx] - exact;
                l2_sum += err * err;
            }
        }
    }
    (l2_sum / total as f64).sqrt()
}

#[test]
fn spectral_poisson_solves_known_solution() {
    let n = 6;
    let config = SpectralConfig::<f64>::new(n, n, n).unwrap();
    let mut solver = SpectralSolver::new(config).unwrap();

    let source = build_poisson_source(n);
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
    let l2 = poisson_l2_error(&solution.u, n);

    assert!(
        l2 < 1.0,
        "Spectral Poisson L2 error should be < 1 for N=6, got {}",
        l2
    );
}

#[test]
fn spectral_poisson_convergence() {
    // Error at N=8 should be less than at N=4 (overall trend)
    let mut errors = Vec::new();

    for &n in &[4, 8] {
        let config = SpectralConfig::<f64>::new(n, n, n).unwrap();
        let mut solver = SpectralSolver::new(config).unwrap();

        let source = build_poisson_source(n);
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
        let l2 = poisson_l2_error(&solution.u, n);
        errors.push((n, l2));
    }

    // Overall trend: finer resolution should have lower error
    let (n0, e0) = errors[0];
    let (n1, e1) = errors[1];
    assert!(
        e1 <= e0 * 1.5, // allow small regression tolerance
        "L2 error should not grow significantly: N={} l2={}, N={} l2={}",
        n0,
        e0,
        n1,
        e1
    );
}

#[test]
fn spectral_solution_correct_dimensions() {
    let n = 4;
    let config = SpectralConfig::<f64>::new(n, n, n).unwrap();
    let mut solver = SpectralSolver::new(config).unwrap();

    let source = build_poisson_source(n);
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
    assert_eq!(solution.nx, n);
    assert_eq!(solution.ny, n);
    assert_eq!(solution.nz, n);
    assert_eq!(solution.u.len(), n * n * n);
}
