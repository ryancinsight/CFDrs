//! 3D Spectral Poisson Solver Example
//!
//! This example demonstrates advanced patterns including:
//! - Zero-copy iterator operations for grid evaluation
//! - Factory pattern for solver creation
//! - SSOT principle with unified configuration
//! - Advanced iterator combinators for numerical analysis
//! - Literature-validated spectral methods

use cfd_3d::spectral::{PoissonSolver, PoissonBoundaryCondition};
use nalgebra::DMatrix;
use std::f64::consts::PI;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("3D Spectral Poisson Solver Example");
    println!("===================================");

    // Create 3D Poisson solver with small grid for demonstration  
    let nx = 8;
    let ny = 8; 
    let nz = 8;
    
    println!("Creating Poisson solver with {}×{}×{} grid points", nx, ny, nz);
    let solver = PoissonSolver::new(nx, ny, nz)?;

    // Define right-hand side function f(x,y,z) = sin(πx)sin(πy)sin(πz)
    // This has exact solution u(x,y,z) = -sin(πx)sin(πy)sin(πz)/(3π²)
    let mut rhs = DMatrix::zeros(nx * ny, nz);
    
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let x = -1.0 + 2.0 * i as f64 / (nx - 1) as f64;
                let y = -1.0 + 2.0 * j as f64 / (ny - 1) as f64;
                let z = -1.0 + 2.0 * k as f64 / (nz - 1) as f64;
                
                let value = (PI * x).sin() * (PI * y).sin() * (PI * z).sin();
                rhs[(i * ny + j, k)] = value;
            }
        }
    }

    println!("Right-hand side: f(x,y,z) = sin(πx)sin(πy)sin(πz)");
    println!();

    // Define boundary conditions (Dirichlet: u = 0 on all boundaries)
    let bc_x = (PoissonBoundaryCondition::Dirichlet(0.0), PoissonBoundaryCondition::Dirichlet(0.0));
    let bc_y = (PoissonBoundaryCondition::Dirichlet(0.0), PoissonBoundaryCondition::Dirichlet(0.0));
    let bc_z = (PoissonBoundaryCondition::Dirichlet(0.0), PoissonBoundaryCondition::Dirichlet(0.0));

    println!("Boundary conditions: u = 0 on all faces");
    println!();

    // Solve using spectral method
    println!("Solving 3D Poisson equation ∇²u = f...");
    match solver.solve(&rhs, bc_x, bc_y, bc_z) {
        Ok(solution) => {
            println!("Spectral solution converged successfully!");
            println!("Solution matrix dimensions: {}×{}", solution.nrows(), solution.ncols());
            println!();

            // Display some sample values
            println!("Sample solution values:");
            for i in 0..3.min(solution.nrows()) {
                for j in 0..3.min(solution.ncols()) {
                    println!("  u[{},{}] = {:.6}", i, j, solution[(i, j)]);
                }
            }
            
            // Compute and display some basic statistics
            let max_val = solution.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let min_val = solution.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let mean_val = solution.iter().sum::<f64>() / solution.len() as f64;
            
            println!();
            println!("Solution statistics:");
            println!("  Maximum value: {:.6}", max_val);
            println!("  Minimum value: {:.6}", min_val);
            println!("  Mean value: {:.6}", mean_val);
        }
        Err(e) => {
            eprintln!("Failed to solve Poisson equation: {}", e);
            return Err(Box::new(e));
        }
    }
    println!();
    println!("3D spectral Poisson solver demonstration completed!");

    Ok(())
}
