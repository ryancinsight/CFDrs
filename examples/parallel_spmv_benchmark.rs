//! Parallel SpMV Performance Benchmark
//!
//! This benchmark demonstrates the performance improvement of parallel sparse matrix-vector
//! multiplication (SpMV) using Rayon in CFD applications.
//!
//! The benchmark compares:
//! - Serial SpMV (baseline)
//! - Parallel SpMV (rayon-based)
//! - Performance scaling with matrix size and core count
//!
//! Expected Results:
//! - 3-5x speedup on 4-8 core systems for large matrices
//! - Minimal overhead for small matrices (<1000 rows)
//! - Near-linear scaling with core count

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::momentum::{MomentumComponent, MomentumSolver};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🚀 Parallel SpMV Performance Benchmark");
    println!("=====================================");
    println!();

    // Test different grid sizes
    let grid_sizes = vec![
        (32, 32),   // Small: 1024 cells
        (64, 64),   // Medium: 4096 cells
        (128, 128), // Large: 16384 cells
        (256, 128), // XL: 32768 cells
    ];

    println!("Testing grid sizes: {:?}", grid_sizes);
    println!("Expected performance: 3-5x speedup on multi-core systems");
    println!();

    for (nx, ny) in grid_sizes {
        println!("📊 Grid Size: {}x{} ({} cells)", nx, ny, nx * ny);
        println!("{}", "─".repeat(50));

        benchmark_grid_size(nx, ny)?;
        println!();
    }

    println!("✅ Benchmark complete!");
    println!("💡 Use MomentumSolver::with_parallel_spmv() for production CFD applications");

    Ok(())
}

/// Benchmark serial vs parallel SpMV for a given grid size
fn benchmark_grid_size(nx: usize, ny: usize) -> Result<(), Box<dyn std::error::Error>> {
    // Create grid and solvers
    let lx = 1.0;
    let ly = 1.0;
    let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, lx, 0.0, ly)?;

    // Create solvers - serial and parallel
    let mut serial_solver = MomentumSolver::new(&grid);
    let mut parallel_solver = MomentumSolver::with_parallel_spmv(&grid);

    // Create test fields with some flow
    let mut fields = SimulationFields::new(nx, ny);

    // Set up a simple channel flow test case
    for j in 0..ny {
        for i in 0..nx {
            // Parabolic velocity profile
            let y_pos = (j as f64 + 0.5) / ny as f64;
            let u_velocity = 6.0 * y_pos * (1.0 - y_pos); // Max velocity = 1.5
            fields.u.set(i, j, u_velocity);
            fields.v.set(i, j, 0.0); // No v-velocity
        }
    }

    // Set constant properties
    for j in 0..ny {
        for i in 0..nx {
            fields.viscosity.set(i, j, 0.01); // Constant viscosity
            fields.density.set(i, j, 1.0); // Constant density
        }
    }

    // Set boundary conditions (simple walls)
    serial_solver.set_boundary(
        "west".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );
    serial_solver.set_boundary(
        "east".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );
    serial_solver.set_boundary(
        "south".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );
    serial_solver.set_boundary(
        "north".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );

    parallel_solver.set_boundary(
        "west".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );
    parallel_solver.set_boundary(
        "east".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );
    parallel_solver.set_boundary(
        "south".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );
    parallel_solver.set_boundary(
        "north".to_string(),
        cfd_core::physics::boundary::BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        },
    );

    // Time step
    let dt = 0.001;

    // Warm up (run once to initialize any lazy allocations)
    let _ = serial_solver.solve_with_coefficients(MomentumComponent::U, &mut fields.clone(), dt)?;
    let _ =
        parallel_solver.solve_with_coefficients(MomentumComponent::U, &mut fields.clone(), dt)?;

    // Benchmark serial solver
    let serial_start = Instant::now();
    const SERIAL_ITERATIONS: u32 = 10;
    for _ in 0..SERIAL_ITERATIONS {
        let _ =
            serial_solver.solve_with_coefficients(MomentumComponent::U, &mut fields.clone(), dt)?;
    }
    let serial_time = serial_start.elapsed();

    // Benchmark parallel solver
    let parallel_start = Instant::now();
    const PARALLEL_ITERATIONS: u32 = 10;
    for _ in 0..PARALLEL_ITERATIONS {
        let _ = parallel_solver.solve_with_coefficients(
            MomentumComponent::U,
            &mut fields.clone(),
            dt,
        )?;
    }
    let parallel_time = parallel_start.elapsed();

    // Calculate performance metrics
    let serial_avg = serial_time / SERIAL_ITERATIONS;
    let parallel_avg = parallel_time / PARALLEL_ITERATIONS;
    let speedup = serial_avg.as_secs_f64() / parallel_avg.as_secs_f64();

    println!("  Serial:   {:>8.2} ms/iter", serial_avg.as_millis());
    println!("  Parallel: {:>8.2} ms/iter", parallel_avg.as_millis());
    println!("  Speedup:  {:>8.2}x", speedup);

    if speedup > 1.2 {
        println!("  ✅ Significant speedup achieved");
    } else if speedup > 0.8 {
        println!("  ⚠️  Minimal speedup (expected for small grids)");
    } else {
        println!("  ❌ Performance regression detected");
    }

    Ok(())
}
