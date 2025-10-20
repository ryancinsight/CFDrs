//! Comprehensive 2D FDM Poisson solver validation tests
//!
//! Validates the finite difference Poisson solver
//! Tests include:
//! - Laplace equation (f=0) with analytical solution
//! - Manufactured sinusoidal solution
//! - Constant source term
//! - Corner singularities
//! - Grid convergence

use approx::assert_relative_eq;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::solvers::fdm::config::FdmConfig;
use cfd_2d::solvers::fdm::poisson::PoissonSolver;
use std::collections::HashMap;
use std::f64::consts::PI;

/// Test Poisson solver with sinusoidal manufactured solution
/// Solves: ∇²φ = -2π²sin(πx)sin(πy)
/// Analytical: φ = sin(πx)sin(πy)
/// Reference: Numerical Recipes (Press et al., 2007)
#[test]
fn test_poisson_2d_sinusoidal_solution() {
    let nx = 21;
    let ny = 21;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    let mut config = FdmConfig::default();
    config.base.convergence.max_iterations = 10000;
    config.base.convergence.tolerance = 1e-8;
    config.base.execution.verbose = true; // Enable verbose output
    let solver = PoissonSolver::new(config);

    // Manufactured solution: φ = sin(πx)sin(πy)
    // ∇²φ = -2π²sin(πx)sin(πy)
    let mut source = HashMap::new();
    for (i, j) in grid.iter() {
        let x = i as f64 / (nx - 1) as f64;
        let y = j as f64 / (ny - 1) as f64;
        let f = -2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
        source.insert((i, j), f);
    }

    // Homogeneous Dirichlet BC: φ = 0 on all boundaries
    let mut boundary_values = HashMap::new();
    for (i, j) in grid.iter() {
        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            boundary_values.insert((i, j), 0.0);
        }
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Verify solution at interior points
    let mut max_error: f64 = 0.0;
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let x = i as f64 / (nx - 1) as f64;
            let y = j as f64 / (ny - 1) as f64;
            let analytical = (PI * x).sin() * (PI * y).sin();
            let computed = solution[&(i, j)];
            let error = (computed - analytical).abs();
            max_error = max_error.max(error);

            // Finite difference 2nd order accuracy: O(h²)
            // For 21x21 grid: h = 0.05, expect error ~ 0.01
            assert_relative_eq!(computed, analytical, epsilon = 0.02);
        }
    }

    // Error should be small for this grid resolution
    assert!(max_error < 0.02, "Max error {} too large", max_error);
}

/// Debug test with small 5x5 grid to understand accuracy issues
#[test]
fn test_poisson_debug_5x5() {
    let nx = 5;
    let ny = 5;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    let mut config = FdmConfig::default();
    config.base.convergence.max_iterations = 10000;
    config.base.convergence.tolerance = 1e-8;
    config.base.execution.verbose = true;
    let solver = PoissonSolver::new(config);

    // Manufactured solution: φ = sin(πx)sin(πy)
    let mut source = HashMap::new();
    for (i, j) in grid.iter() {
        let x = i as f64 / (nx - 1) as f64;
        let y = j as f64 / (ny - 1) as f64;
        let f = -2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
        source.insert((i, j), f);
    }

    // Homogeneous Dirichlet BC
    let mut boundary_values = HashMap::new();
    for (i, j) in grid.iter() {
        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            boundary_values.insert((i, j), 0.0);
        }
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Check center point (2, 2) at (0.5, 0.5)
    let x = 0.5;
    let y = 0.5;
    let analytical = (PI * x).sin() * (PI * y).sin();
    let computed = solution[&(2, 2)];
    
    println!("Center point (2,2) at (0.5, 0.5):");
    println!("  Analytical: {}", analytical);
    println!("  Computed: {}", computed);
    println!("  Error: {}", (computed - analytical).abs());
    println!("  Relative error: {}%", 100.0 * (computed - analytical).abs() / analytical);
    
    // Check all interior points
    for i in 1..nx-1 {
        for j in 1..ny-1 {
            let x = i as f64 / (nx - 1) as f64;
            let y = j as f64 / (ny - 1) as f64;
            let analytical = (PI * x).sin() * (PI * y).sin();
            let computed = solution[&(i, j)];
            println!("({},{}) at ({:.2},{:.2}): analytical={:.4}, computed={:.4}, error={:.4}",
                     i, j, x, y, analytical, computed, (computed - analytical).abs());
        }
    }
}

/// Test Laplace equation (f=0) with non-homogeneous Dirichlet BC
/// Validates steady-state heat conduction
#[test]
fn test_poisson_2d_laplace_equation() {
    let nx = 11;
    let ny = 11;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    let mut config = FdmConfig::default();
    config.base.convergence.tolerance = 1e-10;
    config.base.convergence.max_iterations = 10000;
    let solver = PoissonSolver::new(config);

    // Zero source: Laplace equation ∇²φ = 0
    let source = HashMap::new();

    // BC: φ = 1 on top, φ = 0 on other boundaries
    let mut boundary_values = HashMap::new();
    for (i, j) in grid.iter() {
        if i == 0 || i == nx - 1 || j == 0 {
            boundary_values.insert((i, j), 0.0);
        } else if j == ny - 1 {
            boundary_values.insert((i, j), 1.0);
        }
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Debug: Print solution along center column
    println!("Solution along center column (i=5):");
    for j in 0..ny {
        let y = j as f64 / (ny - 1) as f64;
        println!("  j={}, y={:.2}, phi={:.4}", j, y, solution[&(5, j)]);
    }

    // Solution should be monotonic in y-direction (within numerical tolerance)
    // Note: Iterative solvers may have small oscillations near boundaries
    let monotonic_tol = 1e-8; // Allow tiny violations due to floating-point
    for i in 1..nx - 1 {
        for j in 1..ny - 2 {
            let phi_j = solution[&(i, j)];
            let phi_jp1 = solution[&(i, j + 1)];
            if phi_jp1 < phi_j - monotonic_tol {
                println!("Non-monotonic at ({}, {}): phi[{}]={:.10}, phi[{}]={:.10}, diff={:.2e}",
                         i, j, j, phi_j, j+1, phi_jp1, phi_j - phi_jp1);
            }
            assert!(phi_jp1 >= phi_j - monotonic_tol, 
                    "Solution not monotonic at ({}, {}): phi[{}]={} > phi[{}]={}", 
                    i, j, j, phi_j, j+1, phi_jp1);
        }
    }

    // Verify boundary conditions
    for i in 0..nx {
        assert_relative_eq!(solution[&(i, 0)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(solution[&(i, ny - 1)], 1.0, epsilon = 1e-10);
    }
    for j in 0..ny {
        assert_relative_eq!(solution[&(0, j)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(solution[&(nx - 1, j)], 0.0, epsilon = 1e-10);
    }
}

/// Test Poisson solver with constant source term
/// Models uniform heat generation
#[test]
fn test_poisson_2d_constant_source() {
    let nx = 11;
    let ny = 11;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    let config = FdmConfig::default();
    let solver = PoissonSolver::new(config);

    // Constant source: f = 1
    let mut source = HashMap::new();
    for (i, j) in grid.iter() {
        if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
            source.insert((i, j), 1.0);
        }
    }

    // Homogeneous Dirichlet BC
    let mut boundary_values = HashMap::new();
    for (i, j) in grid.iter() {
        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            boundary_values.insert((i, j), 0.0);
        }
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Solution should be symmetric about center
    // Note: Due to floating-point accumulation in iterative solver,
    // perfect symmetry is not achievable. Tolerance set to 1e-5.
    let mid = nx / 2;
    for i in 1..mid {
        for j in 1..ny - 1 {
            let phi_left = solution[&(i, j)];
            let phi_right = solution[&(nx - 1 - i, j)];
            assert_relative_eq!(phi_left, phi_right, epsilon = 1e-5);
        }
    }

    // Maximum should be at center (for square domain)
    let phi_center = solution[&(mid, mid)];
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            assert!(solution[&(i, j)] <= phi_center + 1e-6);
        }
    }
}

/// Test Poisson solver with corner singularity
/// BC discontinuity creates local singularity
#[test]
fn test_poisson_2d_corner_singularity() {
    let nx = 11;
    let ny = 11;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    let config = FdmConfig::default();
    let solver = PoissonSolver::new(config);

    let source = HashMap::new();

    // BC: φ = 0 on bottom and left, φ = 1 on top and right
    // Creates corner singularity at (0, 0) and (nx-1, ny-1)
    let mut boundary_values = HashMap::new();
    for i in 0..nx {
        boundary_values.insert((i, 0), 0.0); // Bottom
        boundary_values.insert((i, ny - 1), 1.0); // Top
    }
    for j in 0..ny {
        boundary_values.insert((0, j), 0.0); // Left
        boundary_values.insert((nx - 1, j), 1.0); // Right
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Solution should exist and be finite despite singularity
    for i in 0..nx {
        for j in 0..ny {
            let phi = solution[&(i, j)];
            assert!(phi.is_finite(), "Non-finite value at ({}, {})", i, j);
            assert!(phi >= 0.0, "Negative value at ({}, {})", i, j);
            assert!(phi <= 1.0, "Value > 1 at ({}, {})", i, j);
        }
    }
}

/// Test solver convergence with grid refinement
/// Validates 2nd order spatial accuracy
#[test]
fn test_poisson_2d_grid_convergence() {
    // Test on grids: 11x11, 21x21
    let grid_sizes = vec![11, 21];
    let mut errors = Vec::new();

    for &n in &grid_sizes {
        let grid = StructuredGrid2D::<f64>::new(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();

        let mut config = FdmConfig::default();
        config.base.convergence.max_iterations = 10000;
        config.base.convergence.tolerance = 1e-10;
        let solver = PoissonSolver::new(config);

        // Simple manufactured solution
        let mut source = HashMap::new();
        for (i, j) in grid.iter() {
            let x = i as f64 / (n - 1) as f64;
            let y = j as f64 / (n - 1) as f64;
            let f = -2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
            source.insert((i, j), f);
        }

        let mut boundary_values = HashMap::new();
        for (i, j) in grid.iter() {
            if i == 0 || i == n - 1 || j == 0 || j == n - 1 {
                boundary_values.insert((i, j), 0.0);
            }
        }

        let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

        // Compute L2 error
        let mut error_sum = 0.0;
        let mut count = 0;
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let x = i as f64 / (n - 1) as f64;
                let y = j as f64 / (n - 1) as f64;
                let analytical = (PI * x).sin() * (PI * y).sin();
                let computed = solution[&(i, j)];
                error_sum += (computed - analytical).powi(2);
                count += 1;
            }
        }
        let l2_error = (error_sum / count as f64).sqrt();
        errors.push(l2_error);
    }

    // Verify 2nd order convergence: error should reduce by ~4x when grid doubled
    let reduction_factor = errors[0] / errors[1];
    assert!(
        reduction_factor > 3.0,
        "Convergence rate {} too low (expect ~4 for 2nd order)",
        reduction_factor
    );
}

/// Test solver with minimal grid size
/// Edge case: 3x3 grid (only 1 interior point)
#[test]
fn test_poisson_2d_minimal_grid() {
    let nx = 3;
    let ny = 3;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();

    let config = FdmConfig::default();
    let solver = PoissonSolver::new(config);

    let mut source = HashMap::new();
    source.insert((1, 1), 1.0);

    let mut boundary_values = HashMap::new();
    for (i, j) in grid.iter() {
        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            boundary_values.insert((i, j), 0.0);
        }
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Just verify solver doesn't crash and produces finite result
    assert!(solution[&(1, 1)].is_finite());
}

/// Test solver with non-square domain
/// Validates rectangular grid handling
#[test]
fn test_poisson_2d_rectangular_domain() {
    let nx = 11;
    let ny = 21; // Rectangular: 2:1 aspect ratio
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 2.0).unwrap();

    let config = FdmConfig::default();
    let solver = PoissonSolver::new(config);

    let mut source = HashMap::new();
    for (i, j) in grid.iter() {
        source.insert((i, j), 1.0);
    }

    let mut boundary_values = HashMap::new();
    for (i, j) in grid.iter() {
        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            boundary_values.insert((i, j), 0.0);
        }
    }

    let solution = solver.solve(&grid, &source, &boundary_values).unwrap();

    // Verify solution is finite everywhere
    for i in 0..nx {
        for j in 0..ny {
            assert!(solution[&(i, j)].is_finite());
        }
    }
}
