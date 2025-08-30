//! Lid-Driven Cavity - Minimum Viable Product
//! 
//! This example implements a complete, validated lid-driven cavity solver
//! and compares results against Ghia et al. (1982) benchmark data.
//!
//! Reference: Ghia, U., Ghia, K.N., Shin, C.T., 1982. 
//! "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method"
//! Journal of Computational Physics, 48, 387-411.

extern crate nalgebra;
extern crate num_traits;
extern crate cfd_2d;
extern crate cfd_core;

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Ghia et al. (1982) benchmark data for Re=100
/// Format: (y/L, u/U) along vertical centerline x=0.5
const GHIA_RE100_U: &[(f64, f64)] = &[
    (1.0000, 1.00000),  // Lid
    (0.9766, 0.84123),
    (0.9688, 0.78871),
    (0.9609, 0.73722),
    (0.9531, 0.68717),
    (0.8516, 0.23151),
    (0.7344, 0.00332),
    (0.6172, -0.13641),
    (0.5000, -0.20581),
    (0.4531, -0.21090),
    (0.2813, -0.15662),
    (0.1719, -0.10150),
    (0.1016, -0.06434),
    (0.0703, -0.04775),
    (0.0625, -0.04192),
    (0.0547, -0.03717),
    (0.0000, 0.00000),  // Bottom wall
];

/// Ghia et al. (1982) benchmark data for Re=100
/// Format: (x/L, v/U) along horizontal centerline y=0.5
const GHIA_RE100_V: &[(f64, f64)] = &[
    (1.0000, 0.00000),  // Right wall
    (0.9688, -0.05906),
    (0.9609, -0.07391),
    (0.9531, -0.08864),
    (0.9453, -0.10313),
    (0.9063, -0.16914),
    (0.8594, -0.22445),
    (0.8047, -0.24533),
    (0.5000, 0.05454),
    (0.2344, 0.17527),
    (0.2266, 0.17507),
    (0.1563, 0.16077),
    (0.0938, 0.12317),
    (0.0781, 0.10890),
    (0.0703, 0.10091),
    (0.0625, 0.09233),
    (0.0000, 0.00000),  // Left wall
];

fn run_lid_driven_cavity<T: RealField + FromPrimitive + Copy + Send + Sync>(
    nx: usize,
    ny: usize,
    reynolds: T,
    max_iterations: usize,
) -> cfd_core::error::Result<cfd_2d::physics::vorticity_stream::VorticityStreamSolver<T>> {
    println!("=== Lid-Driven Cavity Solver ===");
    println!("Grid: {}x{}", nx, ny);
    println!("Reynolds number: {:.0}", reynolds.to_subset().unwrap_or(0.0));
    println!("Max iterations: {}\n", max_iterations);

    // Create unit square grid
    let grid = cfd_2d::grid::StructuredGrid2D::<T>::unit_square(nx, ny)?;
    
    // Configure solver
    let mut config = cfd_2d::physics::vorticity_stream::VorticityStreamConfig::<T>::default();
    config.base.max_iterations = max_iterations;
    config.base.tolerance = T::from_f64(1e-6).unwrap_or_else(T::zero);
    config.stream_tolerance = T::from_f64(1e-8).unwrap_or_else(T::zero);
    config.vorticity_tolerance = T::from_f64(1e-8).unwrap_or_else(T::zero);
    config.time_step = T::from_f64(0.001).unwrap_or_else(T::zero);
    
    // Create and initialize solver
    let mut solver = cfd_2d::physics::vorticity_stream::VorticityStreamSolver::new(config, &grid, reynolds);
    solver.initialize_lid_driven_cavity(T::one())?;
    
    // Time stepping
    println!("Starting time integration...");
    let mut converged = false;
    let mut iteration = 0;
    
    while iteration < max_iterations && !converged {
        let old_psi = solver.psi.clone(); // Store for convergence check
        
        // Single time step
        solver.step()?;
        
        // Check convergence every 100 iterations
        if iteration % 100 == 0 {
            let mut max_change = T::zero();
            for i in 0..nx {
                for j in 0..ny {
                    let change = (solver.psi[i][j] - old_psi[i][j]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                }
            }
            
            if iteration % 1000 == 0 {
                println!("  Iteration {}: max change = {:.2e}", 
                    iteration, 
                    max_change.to_subset().unwrap_or(0.0));
            }
            
            if max_change < config.base.tolerance {
                converged = true;
                println!("Converged at iteration {}", iteration);
            }
        }
        
        iteration += 1;
    }
    
    if !converged {
        println!("Warning: Did not converge in {} iterations", max_iterations);
    }
    
    Ok(solver)
}

fn validate_against_ghia<T: RealField + Copy>(
    solver: &cfd_2d::physics::vorticity_stream::VorticityStreamSolver<T>,
) -> f64 {
    println!("\n=== Validation Against Ghia et al. (1982) ===");
    
    let nx = solver.nx;
    let ny = solver.ny;
    
    // Extract u-velocity along vertical centerline (x = 0.5)
    println!("\nU-velocity along vertical centerline:");
    println!("  y/L     u_computed  u_ghia     error");
    
    let mut total_error = 0.0;
    let mut n_points = 0;
    
    for &(y_pos, u_ghia) in GHIA_RE100_U.iter() {
        // Find nearest grid point
        let j = ((y_pos * (ny - 1) as f64).round() as usize).min(ny - 1);
        let i = nx / 2; // Centerline
        
        let u_computed = solver.u[i][j].x.to_subset().unwrap_or(0.0);
        let error = (u_computed - u_ghia).abs();
        total_error += error * error;
        n_points += 1;
        
        if y_pos == 1.0 || y_pos == 0.5 || y_pos == 0.0 {
            println!("  {:.4}  {:+.6}  {:+.6}  {:.2e}", 
                y_pos, u_computed, u_ghia, error);
        }
    }
    
    // Extract v-velocity along horizontal centerline (y = 0.5)
    println!("\nV-velocity along horizontal centerline:");
    println!("  x/L     v_computed  v_ghia     error");
    
    for &(x_pos, v_ghia) in GHIA_RE100_V.iter() {
        // Find nearest grid point
        let i = ((x_pos * (nx - 1) as f64).round() as usize).min(nx - 1);
        let j = ny / 2; // Centerline
        
        let v_computed = solver.u[i][j].y.to_subset().unwrap_or(0.0);
        let error = (v_computed - v_ghia).abs();
        total_error += error * error;
        n_points += 1;
        
        if x_pos == 1.0 || x_pos == 0.5 || x_pos == 0.0 {
            println!("  {:.4}  {:+.6}  {:+.6}  {:.2e}", 
                x_pos, v_computed, v_ghia, error);
        }
    }
    
    let rms_error = (total_error / n_points as f64).sqrt();
    println!("\nRMS Error: {:.2e}", rms_error);
    
    // Success criteria: RMS error < 5% for this grid resolution
    if rms_error < 0.05 {
        println!("✅ VALIDATION PASSED");
    } else {
        println!("❌ VALIDATION FAILED - Error too large");
    }
    
    rms_error
}

fn main() -> cfd_core::error::Result<()> {
    println!("╔════════════════════════════════════════╗");
    println!("║   LID-DRIVEN CAVITY - MVP              ║");
    println!("║   Validation against Ghia et al. 1982  ║");
    println!("╚════════════════════════════════════════╝\n");
    
    // Run simulation for Re=100 (matches Ghia benchmark)
    let solver = run_lid_driven_cavity::<f64>(
        41,     // Grid size (41x41 is reasonable for Re=100)
        41,     
        100.0,  // Reynolds number
        50000,  // Max iterations
    )?;
    
    // Validate against benchmark
    let error = validate_against_ghia(&solver);
    
    // Output flow field statistics
    println!("\n=== Flow Field Statistics ===");
    
    // Find maximum vorticity
    let mut max_vorticity = 0.0;
    let mut min_vorticity = 0.0;
    for i in 0..solver.nx {
        for j in 0..solver.ny {
            let omega = solver.omega[i][j];
            if omega > max_vorticity {
                max_vorticity = omega;
            }
            if omega < min_vorticity {
                min_vorticity = omega;
            }
        }
    }
    
    println!("Max vorticity: {:.3}", max_vorticity);
    println!("Min vorticity: {:.3}", min_vorticity);
    
    // Find primary vortex center
    let mut max_psi = f64::NEG_INFINITY;
    let mut vortex_i = 0;
    let mut vortex_j = 0;
    
    for i in 1..solver.nx-1 {
        for j in 1..solver.ny-1 {
            if solver.psi[i][j] > max_psi {
                max_psi = solver.psi[i][j];
                vortex_i = i;
                vortex_j = j;
            }
        }
    }
    
    let vortex_x = vortex_i as f64 / (solver.nx - 1) as f64;
    let vortex_y = vortex_j as f64 / (solver.ny - 1) as f64;
    
    println!("\nPrimary vortex center:");
    println!("  Position: ({:.3}, {:.3})", vortex_x, vortex_y);
    println!("  Stream function: {:.6}", max_psi);
    
    // Expected vortex center for Re=100 from literature: approximately (0.62, 0.74)
    let expected_x = 0.62;
    let expected_y = 0.74;
    let position_error = ((vortex_x - expected_x).powi(2) + 
                          (vortex_y - expected_y).powi(2)).sqrt();
    
    println!("  Expected position: ({:.3}, {:.3})", expected_x, expected_y);
    println!("  Position error: {:.3}", position_error);
    
    if position_error < 0.1 {
        println!("  ✅ Vortex position matches literature");
    } else {
        println!("  ⚠️ Vortex position differs from literature");
    }
    
    println!("\n=== MVP Summary ===");
    if error < 0.05 && position_error < 0.1 {
        println!("✅ PHYSICS VALIDATION SUCCESSFUL");
        println!("The lid-driven cavity solver produces physically correct results");
        println!("that match established benchmarks from literature.");
    } else {
        println!("⚠️ PHYSICS VALIDATION INCOMPLETE");
        println!("The solver runs but needs refinement to match benchmarks.");
    }
    
    Ok(())
}