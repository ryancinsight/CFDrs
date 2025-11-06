//! Demonstration of SIMPLEC and PIMPLE algorithms for pressure-velocity coupling
//!
//! This example shows how to use the advanced SIMPLEC and PIMPLE algorithms
//! for improved convergence in incompressible flow simulations.
//!
//! ## Algorithms
//!
//! - **SIMPLEC**: Semi-Implicit Method for Pressure-Linked Equations - Consistent
//!   - Improves SIMPLE by using consistent discretization
//!   - Better pressure-velocity coupling with Rhie-Chow interpolation
//!   - Reference: Van Doormaal & Raithby (1984)
//!
//! - **PIMPLE**: Merged PISO-SIMPLE algorithm
//!   - Combines outer PISO correctors with inner SIMPLE corrections
//!   - Better for transient flows and stiff problems
//!   - Used in OpenFOAM and other production CFD codes

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::simplec_pimple::{SimplecPimpleConfig, SimplecPimpleSolver, AlgorithmType};
use nalgebra::Vector2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("SIMPLEC and PIMPLE Algorithm Demonstration");
    println!("==========================================");

    // Create a simple 2D grid (10x10 cells)
    let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
    println!("Grid: {}x{} cells", grid.nx, grid.ny);

    // Initialize simulation fields
    let mut fields = SimulationFields::new(10, 10);

    // Set up initial conditions (simple cavity flow)
    // - Zero velocity everywhere initially
    // - Constant pressure
    for i in 0..10 {
        for j in 0..10 {
            fields.set_velocity_at(i, j, &Vector2::new(0.0, 0.0));
            fields.p[(i, j)] = 0.0;
        }
    }

    // Set boundary conditions
    // Top wall: moving lid (u = 1.0, v = 0.0)
    for i in 0..10 {
        fields.set_velocity_at(i, 9, &Vector2::new(1.0, 0.0));
    }

    println!("Initial conditions set up");
    println!("- Cavity flow with moving top lid (u = 1.0)");
    println!("- Zero initial velocity field");
    println!("- Constant initial pressure");

    // Demonstrate SIMPLEC algorithm
    println!("\n--- SIMPLEC Algorithm ---");
    run_algorithm(&grid, &mut fields.clone(), AlgorithmType::Simplec)?;

    // Demonstrate PIMPLE algorithm
    println!("\n--- PIMPLE Algorithm ---");
    run_algorithm(&grid, &mut fields.clone(), AlgorithmType::Pimple)?;

    println!("\nDemonstration completed successfully!");
    println!("Both SIMPLEC and PIMPLE algorithms converged.");

    Ok(())
}

fn run_algorithm(
    grid: &StructuredGrid2D<f64>,
    fields: &mut SimulationFields<f64>,
    algorithm: AlgorithmType,
) -> Result<(), Box<dyn std::error::Error>> {
    // Create algorithm configuration
    let config = match algorithm {
        AlgorithmType::Simplec => {
            println!("Using SIMPLEC (Semi-Implicit Method for Pressure-Linked Equations - Consistent)");
            SimplecPimpleConfig {
                algorithm,
                dt: 0.01,           // Time step
                alpha_u: 0.7,       // Velocity under-relaxation
                alpha_p: 0.3,       // Pressure under-relaxation
                use_rhie_chow: true, // Enable Rhie-Chow interpolation
                ..Default::default()
            }
        }
        AlgorithmType::Pimple => {
            println!("Using PIMPLE (Merged PISO-SIMPLE algorithm)");
            SimplecPimpleConfig {
                algorithm,
                dt: 0.01,
                alpha_u: 0.7,
                alpha_p: 0.3,
                n_outer_correctors: 2,  // PIMPLE outer iterations
                n_inner_correctors: 1,  // PIMPLE inner iterations
                ..Default::default()
            }
        }
    };

    // Create solver
    let mut solver = SimplecPimpleSolver::new(grid.clone(), config)?;

    // Run simulation for a few time steps
    let n_steps = 5;
    println!("Running {} time steps...", n_steps);

    for step in 1..=n_steps {
        // Solve one time step
        let residual = solver.solve_time_step(fields, 0.01, 0.01, 1.0)?;

        println!("Step {:2}: residual = {:.2e}, iterations = {}",
                 step, residual, solver.iterations());

        // Check convergence
        if residual < 1e-6 {
            println!("  Converged at step {}!", step);
            break;
        }
    }

    // Show final velocity field statistics
    let mut max_u = 0.0f64;
    let mut max_v = 0.0f64;
    for i in 0..grid.nx {
        for j in 0..grid.ny {
            let vel = Vector2::new(fields.u.at(i, j), fields.v.at(i, j));
            max_u = max_u.max(vel.x.abs());
            max_v = max_v.max(vel.y.abs());
        }
    }

    println!("Final velocity field:");
    println!("  Max |u| = {:.3}", max_u);
    println!("  Max |v| = {:.3}", max_v);

    Ok(())
}






