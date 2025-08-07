//! 2D Heat Diffusion Example
//!
//! This example demonstrates solving the 2D heat equation using finite differences.
//! We solve: ∇²T = 0 with specified boundary conditions.

use cfd_suite::prelude::*;
use cfd_2d::{GridEdge, AdvectionDiffusionSolver};
use std::collections::HashMap;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== 2D Heat Diffusion Example ===");
    
    // Create a 2D grid for a rectangular domain
    let mut grid = StructuredGrid2D::<f64>::new(
        20,    // nx: 20 cells in x-direction
        15,    // ny: 15 cells in y-direction
        0.0,   // x_min
        2.0,   // x_max
        0.0,   // y_min
        1.5,   // y_max
    )?;
    
    println!("Created {}x{} grid", grid.nx(), grid.ny());
    println!("Grid spacing: {:?}", grid.spacing());
    
    // Set boundary conditions for heat diffusion
    // Left wall: T = 100°C (hot)
    // Right wall: T = 0°C (cold)
    // Top and bottom walls: insulated (∂T/∂n = 0, approximated as T = neighbor)
    
    grid.set_edge_boundary(GridEdge::Left, BoundaryType::Wall);
    grid.set_edge_boundary(GridEdge::Right, BoundaryType::Wall);
    grid.set_edge_boundary(GridEdge::Bottom, BoundaryType::Wall);
    grid.set_edge_boundary(GridEdge::Top, BoundaryType::Wall);
    
    // Set boundary values
    let mut boundary_values = HashMap::new();
    
    // Hot left wall
    for j in 0..grid.ny() {
        boundary_values.insert((0, j), 100.0);
    }
    
    // Cold right wall
    for j in 0..grid.ny() {
        boundary_values.insert((grid.nx() - 1, j), 0.0);
    }
    
    // Insulated top and bottom walls (simplified as fixed temperature)
    // Note: This is a simplified approach. For true Neumann boundary conditions
    // (zero gradient), we would need to modify the finite difference stencil
    // to enforce ∂T/∂n = 0 at the boundary.
    for i in 1..grid.nx()-1 {
        boundary_values.insert((i, 0), 50.0);           // Bottom wall (simplified)
        boundary_values.insert((i, grid.ny() - 1), 50.0); // Top wall (simplified)
    }
    
    // No source term for steady-state heat conduction
    let source = HashMap::new();
    
    // Configure the Poisson solver
    let config = FdmConfig {
        tolerance: 1e-8,
        max_iterations: 2000,
        relaxation_factor: 1.0,
        verbose: true,
    };
    
    // Solve the heat equation
    println!("\nSolving 2D heat equation...");
    let solver = PoissonSolver::new(config);
    let solution = solver.solve(&grid, &source, &boundary_values)?;
    
    // Analyze results
    println!("\n=== Results Analysis ===");
    
    // Find temperature statistics
    let mut min_temp = f64::INFINITY;
    let mut max_temp = f64::NEG_INFINITY;
    let mut sum_temp = 0.0;
    let mut count = 0;
    
    for (i, j) in grid.iter() {
        if let Some(&temp) = solution.get(&(i, j)) {
            min_temp = min_temp.min(temp);
            max_temp = max_temp.max(temp);
            sum_temp += temp;
            count += 1;
        }
    }
    
    let avg_temp = sum_temp / count as f64;
    
    println!("Temperature range: {:.2}°C to {:.2}°C", min_temp, max_temp);
    println!("Average temperature: {:.2}°C", avg_temp);
    
    // Print temperature distribution along centerline
    println!("\n=== Temperature along horizontal centerline ===");
    let center_j = grid.ny() / 2;
    println!("Position (x)\tTemperature (°C)");
    
    for i in 0..grid.nx() {
        if let Some(&temp) = solution.get(&(i, center_j)) {
            let center = grid.cell_center(i, center_j)?;
            println!("{:.3}\t\t{:.2}", center.x, temp);
        }
    }
    
    // Calculate heat flux at boundaries (approximate)
    println!("\n=== Heat Flux Analysis ===");
    let (dx, _dy) = grid.spacing();
    
    // Heat flux at left boundary (hot wall)
    let mut total_heat_flux_left = 0.0;
    for j in 1..grid.ny()-1 {
        if let (Some(&t_boundary), Some(&t_interior)) = 
            (solution.get(&(0, j)), solution.get(&(1, j))) {
            let heat_flux = -(t_interior - t_boundary) / dx; // q = -k * dT/dx
            total_heat_flux_left += heat_flux;
        }
    }
    
    println!("Total heat flux from left wall: {:.2} W/m²", total_heat_flux_left);
    
    // Demonstrate advection-diffusion solver with flow
    println!("\n=== Advection-Diffusion Example ===");
    
    // Create velocity field (simple uniform flow from left to right)
    let mut velocity_x = HashMap::new();
    let mut velocity_y = HashMap::new();
    
    for (i, j) in grid.iter() {
        velocity_x.insert((i, j), 1.0); // 1 m/s in x-direction
        velocity_y.insert((i, j), 0.0); // No y-velocity
    }
    
    // Solve advection-diffusion equation
    let diffusivity = 0.1; // Thermal diffusivity
    
    let advdiff_config = FdmConfig {
        tolerance: 1e-6,
        max_iterations: 1000,
        relaxation_factor: 0.8,
        verbose: false,
    };
    
    let advdiff_solver = AdvectionDiffusionSolver::new(advdiff_config);
    let advdiff_solution = advdiff_solver.solve_steady(
        &grid,
        &velocity_x,
        &velocity_y,
        diffusivity,
        &source,
        &boundary_values,
    )?;
    
    // Compare solutions
    println!("\n=== Comparison: Pure Diffusion vs Advection-Diffusion ===");
    println!("Position (x)\tDiffusion\tAdv-Diff\tDifference");
    
    for i in (0..grid.nx()).step_by(4) {
        if let (Some(&temp_diff), Some(&temp_advdiff)) = 
            (solution.get(&(i, center_j)), advdiff_solution.get(&(i, center_j))) {
            let center = grid.cell_center(i, center_j)?;
            let difference = temp_advdiff - temp_diff;
            println!("{:.3}\t\t{:.2}\t\t{:.2}\t\t{:.2}", 
                     center.x, temp_diff, temp_advdiff, difference);
        }
    }
    
    println!("\n=== Example Complete ===");
    println!("This example demonstrated:");
    println!("1. Creating structured 2D grids");
    println!("2. Setting boundary conditions");
    println!("3. Solving Poisson equation (heat diffusion)");
    println!("4. Solving advection-diffusion equation");
    println!("5. Analyzing and comparing results");
    
    Ok(())
}
