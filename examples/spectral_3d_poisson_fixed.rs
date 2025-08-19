//! Example demonstrating spectral methods for solving 3D Poisson equation
//!
//! This example solves the Poisson equation ∇²u = f using Chebyshev spectral methods

use cfd_3d::spectral::{SpectralConfig, SpectralSolver};
use nalgebra::Vector3;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("===========================================");
    println!("Spectral Method: 3D Poisson Equation");
    println!("===========================================\n");
    
    // Create configuration for spectral solver
    // Using small grid for demonstration (increase for higher accuracy)
    // Spectral methods achieve exponential convergence for smooth solutions
    
    println!("Setting up spectral configuration...");
    
    let base = cfd_core::SolverConfig::<f64>::builder()
        .tolerance(1e-8)
        .max_iterations(100)
        .verbosity(2) // verbose = true means verbosity level 2
        .build_base();

    let config = SpectralConfig {
        base,
        nx_modes: 8,
        ny_modes: 8,
        nz_modes: 8,
        dt: None,
    };
    
    println!("Spectral configuration:");
    println!("  Grid: {}×{}×{} modes", config.nx_modes, config.ny_modes, config.nz_modes);
    println!("  Tolerance: {:.2e}", config.tolerance());
    println!("  Max iterations: {}", config.max_iterations());
    
    // Create spectral solver
    let mut solver = SpectralSolver::new(config.clone());
    
    println!("\nInitializing test problem...");
    
    // Set up a test problem: u = sin(πx)sin(πy)sin(πz)
    // This gives f = -3π²sin(πx)sin(πy)sin(πz)
    let nx = config.nx_modes;
    let ny = config.ny_modes;
    let nz = config.nz_modes;
    
    // Initialize source term f
    let mut source = vec![vec![vec![0.0; nz]; ny]; nx];
    let mut exact_solution = vec![vec![vec![0.0; nz]; ny]; nx];
    
    // Generate Chebyshev collocation points
    let x_points: Vec<f64> = (0..nx).map(|i| {
        -((PI * i as f64) / (nx as f64 - 1.0)).cos()
    }).collect();
    
    let y_points: Vec<f64> = (0..ny).map(|j| {
        -((PI * j as f64) / (ny as f64 - 1.0)).cos()
    }).collect();
    
    let z_points: Vec<f64> = (0..nz).map(|k| {
        -((PI * k as f64) / (nz as f64 - 1.0)).cos()
    }).collect();
    
    println!("Setting up source term and exact solution...");
    
    // Map Chebyshev points from [-1, 1] to [0, 1] for our test function
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // Map from [-1, 1] to [0, 1]
                let x = (x_points[i] + 1.0) / 2.0;
                let y = (y_points[j] + 1.0) / 2.0;
                let z = (z_points[k] + 1.0) / 2.0;
                
                // Exact solution: u = sin(πx)sin(πy)sin(πz)
                exact_solution[i][j][k] = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                
                // Source term: f = -3π²u
                source[i][j][k] = -3.0 * PI * PI * exact_solution[i][j][k];
            }
        }
    }
    
    println!("Source term statistics:");
    let max_source = source.iter()
        .flat_map(|plane| plane.iter())
        .flat_map(|row| row.iter())
        .fold(0.0f64, |max, &val| max.max(val.abs()));
    println!("  Max |f|: {:.6}", max_source);
    
    // Solve Poisson equation
    println!("\nSolving Poisson equation ∇²u = f...");
    
    let solution = solver.solve_poisson(&source)?;
    
    println!("Solution obtained!");
    
    // Compute error
    println!("\nComputing error metrics...");
    
    let mut max_error = 0.0;
    let mut l2_error = 0.0;
    let mut l2_norm = 0.0;
    
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let error = (solution[i][j][k] - exact_solution[i][j][k]).abs();
                max_error = max_error.max(error);
                l2_error += error * error;
                l2_norm += exact_solution[i][j][k] * exact_solution[i][j][k];
            }
        }
    }
    
    l2_error = (l2_error / (nx * ny * nz) as f64).sqrt();
    l2_norm = (l2_norm / (nx * ny * nz) as f64).sqrt();
    let relative_l2_error = if l2_norm > 1e-10 { l2_error / l2_norm } else { l2_error };
    
    println!("\n===========================================");
    println!("Results:");
    println!("===========================================");
    println!("  Max error: {:.6e}", max_error);
    println!("  L2 error: {:.6e}", l2_error);
    println!("  Relative L2 error: {:.6e}", relative_l2_error);
    
    // Verify spectral convergence
    if max_error < 1e-6 {
        println!("\n✓ Solution converged to expected accuracy!");
        println!("  Spectral methods achieve exponential convergence");
        println!("  for smooth solutions like this test case.");
    } else {
        println!("\n⚠ Higher error than expected.");
        println!("  Try increasing the number of modes for better accuracy.");
    }
    
    // Sample the solution at a few points
    println!("\nSolution samples:");
    let sample_points = vec![
        (nx/4, ny/4, nz/4),
        (nx/2, ny/2, nz/2),
        (3*nx/4, 3*ny/4, 3*nz/4),
    ];
    
    for (i, j, k) in sample_points {
        let x = (x_points[i] + 1.0) / 2.0;
        let y = (y_points[j] + 1.0) / 2.0;
        let z = (z_points[k] + 1.0) / 2.0;
        
        println!("  u({:.2}, {:.2}, {:.2}):", x, y, z);
        println!("    Computed: {:.6}", solution[i][j][k]);
        println!("    Exact:    {:.6}", exact_solution[i][j][k]);
        println!("    Error:    {:.6e}", (solution[i][j][k] - exact_solution[i][j][k]).abs());
    }
    
    println!("\n===========================================");
    println!("Spectral method demonstration complete!");
    println!("===========================================");
    
    Ok(())
}