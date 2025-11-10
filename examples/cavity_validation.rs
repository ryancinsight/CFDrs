//! Lid-Driven Cavity Validation Against Ghia et al. (1982)
//!
//! This example validates the vorticity-stream solver against established benchmarks.

use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::vorticity_stream::{VorticityStreamConfig, VorticityStreamSolver};
use cfd_core::error::Result;

/// Ghia et al. (1982) benchmark data for Re=100
/// Vertical centerline u-velocity at selected points
const GHIA_RE100_U: &[(f64, f64)] = &[
    (1.0000, 1.00000), // Lid
    (0.9531, 0.68717),
    (0.5000, -0.20581), // Center
    (0.0547, -0.03717),
    (0.0000, 0.00000), // Bottom
];

fn main() -> Result<()> {
    println!("=== LID-DRIVEN CAVITY VALIDATION ===\n");

    // Create 41x41 grid (standard for Re=100)
    let nx = 41;
    let ny = 41;
    let reynolds = 100.0;

    println!("Configuration:");
    println!("  Grid: {}x{}", nx, ny);
    println!("  Reynolds: {}", reynolds);

    // Create unit square grid
    let grid = StructuredGrid2D::<f64>::unit_square(nx, ny)?;

    // Configure solver
    let mut config = VorticityStreamConfig::<f64>::default();
    config.base.convergence.max_iterations = 10000;
    config.base.convergence.tolerance = 1e-6;
    config.time_step = 0.001;

    // Create and initialize solver
    let mut solver = VorticityStreamSolver::new(config, &grid, reynolds);
    solver.initialize_lid_driven_cavity(1.0)?;

    println!("\nRunning simulation...");

    // Time stepping
    let mut iteration = 0;
    let max_iterations = 1000; // Reduced for quick validation

    while iteration < max_iterations {
        solver.step()?;

        if iteration % 200 == 0 {
            println!("  Iteration {}", iteration);
        }

        iteration += 1;
    }

    println!("\n=== VALIDATION RESULTS ===\n");

    // Extract centerline velocity
    let i_center = nx / 2;

    println!("Vertical centerline u-velocity comparison:");
    println!("  y/L     u_computed  u_ghia     error");
    println!("  ----    ----------  -------    -----");

    let velocity_field = solver.velocity_field();
    let mut total_error = 0.0;
    for &(y_pos, u_ghia) in GHIA_RE100_U.iter() {
        let j = ((y_pos * (ny - 1) as f64).round() as usize).min(ny - 1);
        let u_computed = velocity_field[i_center][j].x;
        let error = (u_computed - u_ghia).abs();
        total_error += error * error;

        println!(
            "  {:.4}  {:+.6}  {:+.6}  {:.4}",
            y_pos, u_computed, u_ghia, error
        );
    }

    let rms_error = (total_error / GHIA_RE100_U.len() as f64).sqrt();

    println!("\nRMS Error: {:.4}", rms_error);

    // Check for physics validity
    let stream_function = solver.stream_function();
    let mut has_circulation = false;
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            if stream_function[i][j].abs() > 1e-3 {
                has_circulation = true;
                break;
            }
        }
    }

    println!("\n=== PHYSICS CHECK ===");
    if has_circulation {
        println!("✅ Flow circulation detected");
    } else {
        println!("❌ No flow circulation - solver may not be working");
    }

    if rms_error < 0.1 {
        println!("✅ Validation PASSED - Error within acceptable range");
    } else if rms_error < 0.5 {
        println!("⚠️  Validation MARGINAL - Solver needs tuning");
    } else {
        println!("❌ Validation FAILED - Physics incorrect");
    }

    Ok(())
}
