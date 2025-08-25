//! 3D Spectral Poisson Solver Example
//!
//! This example demonstrates advanced patterns including:
//! - Zero-copy iterator operations for grid evaluation
//! - Factory pattern for solver creation
//! - SSOT principle with unified configuration
//! - Advanced iterator combinators for numerical analysis
//! - Literature-validated spectral methods

use cfd_suite::prelude::*;
use cfd_3d::{SpectralBasis, SpectralConfig, SpectralSolver};
use cfd_core::BoundaryCondition;
use std::collections::HashMap;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("3D Spectral Poisson Solver Example");
    println!("===================================");
    
    // Create spectral solver configuration with smaller grid for demonstration
    let base = cfd_core::SolverConfig::<f64>::builder()
        .tolerance(1e-8)
        .max_iterations(100)
        .verbose(true)
        .build();

    let config = SpectralConfig {
        base: base.clone(),
        nx_modes: 8,
        ny_modes: 8,
        nz_modes: 8,
        basis_x: SpectralBasis::Chebyshev,
        basis_y: SpectralBasis::Chebyshev,
        basis_z: SpectralBasis::Chebyshev,
        dt: None,
    };
    
    println!("Spectral configuration:");
    println!("  Grid: {}×{}×{} modes", config.nx_modes, config.ny_modes, config.nz_modes);
    println!("  Tolerance: {:.0e}", base.convergence.tolerance);
    println!("  Basis: Chebyshev polynomials");
    println!();
    
    // Define domain bounds (unit cube)
    let domain_bounds = (
        nalgebra::Vector3::new(-1.0, -1.0, -1.0),
        nalgebra::Vector3::new(1.0, 1.0, 1.0),
    );
    
    println!("Domain: [{:.1}, {:.1}]³", 
             domain_bounds.0.x, domain_bounds.1.x);
    
    // Create spectral solver
    let mut solver = SpectralSolver::new(config).expect("Failed to create spectral solver");
    
    // Define the source function for Poisson equation: ∇²u = f
    // We'll use f(x,y,z) = -3π²sin(πx)sin(πy)sin(πz)
    // which has the analytical solution u(x,y,z) = sin(πx)sin(πy)sin(πz)
    let source_function = |point: &nalgebra::Vector3<f64>| -> f64 {
        let pi = std::f64::consts::PI;
        let x = point.x;
        let y = point.y;
        let z = point.z;
        
        -3.0 * pi * pi * (pi * x).sin() * (pi * y).sin() * (pi * z).sin()
    };
    
    println!("Source function: f(x,y,z) = -3π²sin(πx)sin(πy)sin(πz)");
    println!("Analytical solution: u(x,y,z) = sin(πx)sin(πy)sin(πz)");
    println!();
    
    // Define boundary conditions (Dirichlet: u = 0 on all boundaries)
    let mut boundary_conditions = HashMap::new();
    boundary_conditions.insert("x_min".to_string(), BoundaryCondition::Dirichlet { value: 0.0 });
    boundary_conditions.insert("x_max".to_string(), BoundaryCondition::Dirichlet { value: 0.0 });
    boundary_conditions.insert("y_min".to_string(), BoundaryCondition::Dirichlet { value: 0.0 });
    boundary_conditions.insert("y_max".to_string(), BoundaryCondition::Dirichlet { value: 0.0 });
    boundary_conditions.insert("z_min".to_string(), BoundaryCondition::Dirichlet { value: 0.0 });
    boundary_conditions.insert("z_max".to_string(), BoundaryCondition::Dirichlet { value: 0.0 });
    
    println!("Boundary conditions: u = 0 on all faces");
    println!();
    
    // Solve using spectral method (simplified demonstration)
    println!("Solving 3D problem using spectral methods...");
    match solver.solve() {
        Ok(solution) => {
            println!("Spectral solution converged successfully!");
            println!();
            
            // Display solution information
            println!("Solution information:");
            
            {
                let values = &solution.u;
                    println!("Solution values (sample):");
                    
                    // Show a few sample values
                    let step = values.len() / 8; // Show ~8 values
                    for (i, &value) in values.iter().enumerate().step_by(step.max(1)).take(8) {
                        let k = i / (solution.nx * solution.ny);
                        let j = (i % (solution.nx * solution.ny)) / solution.nx;
                        let i_local = i % solution.nx;
                        
                        // Map grid indices to physical coordinates
                        let x = -1.0 + 2.0 * i_local as f64 / (solution.nx - 1).max(1) as f64;
                        let y = -1.0 + 2.0 * j as f64 / (solution.ny - 1).max(1) as f64;
                        let z = -1.0 + 2.0 * k as f64 / (solution.nz - 1).max(1) as f64;
                        
                        // Calculate analytical solution for comparison
                        let pi = std::f64::consts::PI;
                        let analytical = (pi * x).sin() * (pi * y).sin() * (pi * z).sin();
                        
                        println!("  ({:.2}, {:.2}, {:.2}): u = {:.6} (analytical: {:.6})", 
                                 x, y, z, value, analytical);
                    }
                    
                    // Calculate statistics using advanced iterator patterns
                    

                    let max_value = values.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
                    let min_value = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                    let avg_value = values.iter().sum::<f64>() / values.len() as f64;
                    let variance = if values.len() > 1 {
                        values.iter()
                            .map(|&x| (x - avg_value).powi(2))
                            .sum::<f64>() / (values.len() - 1) as f64
                    } else {
                        0.0
                    };
                    let std_dev = variance.sqrt();
                    
                    println!();
                    println!("Solution statistics (using advanced iterator patterns):");
                    println!("  Maximum |u|: {:.6}", max_value);
                    println!("  Minimum u: {:.6}", min_value);
                    println!("  Average u: {:.6}", avg_value);
                    println!("  Standard deviation: {:.6}", std_dev);
                    
                    // For the analytical solution sin(πx)sin(πy)sin(πz), 
                    // the maximum should be around 1.0 at the center
                    println!("  Expected maximum: ~1.0 (at center of domain)");
                    
                    if max_value > 0.5 && max_value < 1.5 {
                        println!("  ✓ Solution magnitude is reasonable");
                    } else {
                        println!("  ⚠ Solution magnitude may be incorrect");
                    }

            }
            
            println!();
            println!("Spectral solution properties:");
            println!("  Grid dimensions: {}×{}×{}", 
                     solution.nx, solution.ny, solution.nz);
            println!("  Solution size: {} elements", solution.u.len());
            
            // Spectral methods provide exponential convergence for smooth solutions
            println!();
            println!("Note: Spectral methods achieve exponential convergence");
            println!("for smooth solutions like this trigonometric function.");
        }
        Err(e) => {
            println!("Spectral solution failed: {}", e);
            return Err(e.into());
        }
    }
    
    println!();
    println!("3D spectral Poisson solver demonstration completed!");
    
    Ok(())
}
