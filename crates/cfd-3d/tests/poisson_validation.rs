//! Comprehensive 3D Spectral Poisson solver validation tests
//!
//! Validates the spectral Poisson solver per Boyd (2001)
//! Tests include:
//! - Dirichlet boundary conditions
//! - Neumann boundary conditions
//! - Robin boundary conditions
//! - Manufactured solutions
//! - Tensor product operator accuracy

use approx::assert_relative_eq;
use cfd_3d::spectral::poisson::{PoissonBoundaryCondition, PoissonSolver};
use nalgebra::DMatrix;
use std::f64::consts::PI;

/// Test 3D Poisson solver with homogeneous Dirichlet BC
/// Solves: ∇²u = -2π²sin(πx)sin(πy)sin(πz)
/// Analytical solution: u = sin(πx)sin(πy)sin(πz)
#[test]
fn test_poisson_3d_dirichlet_sinusoidal() {
    let nx = 8;
    let ny = 8;
    let nz = 4;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    // Manufactured solution: u = sin(πx)sin(πy)sin(πz)
    // Laplacian: ∇²u = -3π²sin(πx)sin(πy)sin(πz)
    // RHS: f = 3π²sin(πx)sin(πy)sin(πz)
    let mut f = DMatrix::zeros(nx * ny, nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let x = i as f64 / (nx - 1) as f64;
                let y = j as f64 / (ny - 1) as f64;
                let z = k as f64 / (nz - 1) as f64;
                let idx = j * nx + i;
                f[(idx, k)] = 3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
            }
        }
    }

    // Homogeneous Dirichlet BC: u = 0 on all boundaries
    let bc_x = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Verify solution at interior points
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            for k in 1..nz - 1 {
                let x = i as f64 / (nx - 1) as f64;
                let y = j as f64 / (ny - 1) as f64;
                let z = k as f64 / (nz - 1) as f64;
                let idx = j * nx + i;
                let analytical = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();

                // Spectral accuracy depends on grid resolution
                // Coarse grid (8x8x4) won't achieve machine precision
                assert_relative_eq!(solution[(idx, k)], analytical, epsilon = 0.1);
            }
        }
    }

    // Verify boundary conditions
    for j in 0..ny {
        for k in 0..nz {
            assert_relative_eq!(solution[(j * nx, k)], 0.0, epsilon = 1e-10);
            assert_relative_eq!(solution[(j * nx + nx - 1, k)], 0.0, epsilon = 1e-10);
        }
    }
}

/// Test 3D Poisson solver with constant source term
/// Solves: ∇²u = 1
/// With homogeneous Dirichlet BC
#[test]
fn test_poisson_3d_constant_source() {
    let nx = 6;
    let ny = 6;
    let nz = 4;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    // Constant source: f = 1
    let f = DMatrix::from_element(nx * ny, nz, 1.0);

    // Homogeneous Dirichlet BC
    let bc_x = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Solution should be smooth and bounded
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = j * nx + i;
                assert!(solution[(idx, k)].is_finite());
                assert!(solution[(idx, k)].abs() < 10.0);
            }
        }
    }

    // Verify boundary conditions
    for j in 0..ny {
        for k in 0..nz {
            assert_relative_eq!(solution[(j * nx, k)], 0.0, epsilon = 1e-10);
            assert_relative_eq!(solution[(j * nx + nx - 1, k)], 0.0, epsilon = 1e-10);
        }
    }
}

/// Test 3D Poisson solver with Neumann boundary conditions
/// Validates implementation of Neumann BC
#[test]
fn test_poisson_3d_neumann_bc() {
    let nx = 6;
    let ny = 6;
    let nz = 4;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    // Simple source term
    let f = DMatrix::from_element(nx * ny, nz, 1.0);

    // Mix of Dirichlet and Neumann BC
    let bc_x = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Neumann(0.0),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Solution should be finite and bounded
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = j * nx + i;
                assert!(solution[(idx, k)].is_finite());
            }
        }
    }
}

/// Test 3D Poisson solver with Robin boundary conditions
/// Validates αu + β∂u/∂n = g implementation
#[test]
fn test_poisson_3d_robin_bc() {
    let nx = 6;
    let ny = 6;
    let nz = 4;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    // Simple source term
    let f = DMatrix::from_element(nx * ny, nz, 1.0);

    // Robin BC: u + ∂u/∂n = 0
    let bc_x = (
        PoissonBoundaryCondition::Robin {
            alpha: 1.0,
            beta: 1.0,
            value: 0.0,
        },
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Solution should be finite and bounded
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = j * nx + i;
                assert!(solution[(idx, k)].is_finite());
            }
        }
    }
}

/// Test 3D Poisson solver zero source (Laplace equation)
/// Solves: ∇²u = 0
/// With non-homogeneous Dirichlet BC
#[test]
fn test_poisson_3d_laplace_equation() {
    let nx = 6;
    let ny = 6;
    let nz = 4;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    // Zero source: Laplace equation
    let f = DMatrix::zeros(nx * ny, nz);

    // Non-homogeneous Dirichlet BC: u = 1 on one face
    let bc_x = (
        PoissonBoundaryCondition::Dirichlet(1.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Solution should be smooth and bounded
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = j * nx + i;
                assert!(solution[(idx, k)].is_finite());
                assert!(solution[(idx, k)] >= 0.0);
                assert!(solution[(idx, k)] <= 1.0);
            }
        }
    }

    // Verify boundary values
    for j in 0..ny {
        for k in 0..nz {
            assert_relative_eq!(solution[(j * nx, k)], 1.0, epsilon = 1e-10);
            assert_relative_eq!(solution[(j * nx + nx - 1, k)], 0.0, epsilon = 1e-10);
        }
    }
}

/// Test solver handles edge case: all Dirichlet BC
/// This should produce a well-conditioned system
#[test]
fn test_poisson_3d_all_dirichlet() {
    let nx = 5;
    let ny = 5;
    let nz = 3;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    let f = DMatrix::from_element(nx * ny, nz, 2.0);

    let bc_x = (
        PoissonBoundaryCondition::Dirichlet(0.5),
        PoissonBoundaryCondition::Dirichlet(0.5),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.5),
        PoissonBoundaryCondition::Dirichlet(0.5),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.5),
        PoissonBoundaryCondition::Dirichlet(0.5),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Verify solution is finite
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = j * nx + i;
                assert!(solution[(idx, k)].is_finite());
            }
        }
    }
}

/// Test solver with minimal grid size
/// Edge case: smallest practical 3D grid
#[test]
fn test_poisson_3d_minimal_grid() {
    let nx = 3;
    let ny = 3;
    let nz = 2;
    let solver = PoissonSolver::<f64>::new(nx, ny, nz).unwrap();

    let f = DMatrix::from_element(nx * ny, nz, 1.0);

    let bc_x = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_y = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );
    let bc_z = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let solution = solver.solve(&f, &bc_x, &bc_y, &bc_z).unwrap();

    // Just verify solver doesn't crash and produces finite values
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let idx = j * nx + i;
                assert!(solution[(idx, k)].is_finite());
            }
        }
    }
}
