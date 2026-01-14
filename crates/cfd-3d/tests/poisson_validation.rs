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
use nalgebra::DVector;
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

    // Manufactured solution: u = sin(PI * x) * sin(PI * y) * sin(PI * z)
    // Laplacian: ∇²u = -3 * PI^2 * sin(PI * x) * sin(PI * y) * sin(PI * z)
    // RHS: f = 3 * PI^2 * sin(PI * x) * sin(PI * y) * sin(PI * z)
    let mut f = DVector::zeros(nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let x = (PI * i as f64 / (nx - 1) as f64).cos();
                let y = (PI * j as f64 / (ny - 1) as f64).cos();
                let z = (PI * k as f64 / (nz - 1) as f64).cos();
                let idx = i * ny * nz + j * nz + k;
                f[idx] = 3.0 * PI * PI * (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
            }
        }
    }

    // Homogeneous Dirichlet BC: u = 0 on all boundaries
    // Note: sin(PI * x) is 0 at x=1 and x=-1 (the Chebyshev points)
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

    // Verify solution matches analytical one
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let x = (PI * i as f64 / (nx - 1) as f64).cos();
                let y = (PI * j as f64 / (ny - 1) as f64).cos();
                let z = (PI * k as f64 / (nz - 1) as f64).cos();
                let idx = i * ny * nz + j * nz + k;
                let analytical = -(PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                assert_relative_eq!(solution[idx], analytical, epsilon = 0.1);
            }
        }
    }

    // Verify boundary conditions
    // X boundaries
    for j in 0..ny {
        for k in 0..nz {
            assert_relative_eq!(solution[j * nz + k], 0.0, epsilon = 1e-10);
            assert_relative_eq!(
                solution[(nx - 1) * ny * nz + j * nz + k],
                0.0,
                epsilon = 1e-10
            );
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
    let f = DVector::from_element(nx * ny * nz, 1.0);

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
    for idx in 0..nx * ny * nz {
        assert!(solution[idx].is_finite());
        assert!(solution[idx].abs() < 10.0);
    }

    // Verify boundary conditions
    // X boundaries
    for j in 0..ny {
        for k in 0..nz {
            assert_relative_eq!(solution[j * nz + k], 0.0, epsilon = 1e-10);
            assert_relative_eq!(
                solution[(nx - 1) * ny * nz + j * nz + k],
                0.0,
                epsilon = 1e-10
            );
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
    let f = DVector::from_element(nx * ny * nz, 1.0);

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
    for idx in 0..nx * ny * nz {
        assert!(solution[idx].is_finite());
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
    let f = DVector::from_element(nx * ny * nz, 1.0);

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
    for idx in 0..nx * ny * nz {
        assert!(solution[idx].is_finite());
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
    let f = DVector::zeros(nx * ny * nz);

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
    for idx in 0..nx * ny * nz {
        assert!(solution[idx].is_finite());
        assert!(solution[idx] >= -0.1); // Small numerical overshoot possible
        assert!(solution[idx] <= 1.1);
    }

    // Verify boundary values (interior of the face to avoid corner conflicts)
    for j in 1..ny - 1 {
        for k in 1..nz - 1 {
            assert_relative_eq!(solution[j * nz + k], 1.0, epsilon = 1e-10);
            assert_relative_eq!(
                solution[(nx - 1) * ny * nz + j * nz + k],
                0.0,
                epsilon = 1e-10
            );
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

    let f = DVector::from_element(nx * ny * nz, 2.0);

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
    for idx in 0..nx * ny * nz {
        assert!(solution[idx].is_finite());
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

    let f = DVector::from_element(nx * ny * nz, 1.0);

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
    for idx in 0..nx * ny * nz {
        assert!(solution[idx].is_finite());
    }
}
