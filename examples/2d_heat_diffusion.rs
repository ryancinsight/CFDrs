//! 2D Heat Diffusion Example
//!
//! This example demonstrates solving the 2D heat equation using finite differences.
//! We solve: ∇²T = 0 with specified boundary conditions.

// Correct imports from specific modules
use cfd_2d::grid::{BoundaryType, StructuredGrid2D};
use cfd_2d::solvers::fdm::{AdvectionDiffusionSolver, FdmConfig, PoissonSolver};
use cfd_core::prelude::SolverConfig;
use std::collections::HashMap;

type GridScalarField = HashMap<(usize, usize), f64>;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("=== 2D Heat Diffusion Example ===");

    // Create a 2D grid for a rectangular domain
    let mut grid = StructuredGrid2D::<f64>::new(
        20,  // nx: 20 cells in x-direction
        15,  // ny: 15 cells in y-direction
        0.0, // x_min
        2.0, // x_max
        0.0, // y_min
        1.5, // y_max
    )?;

    println!("Created {}x{} grid", grid.nx, grid.ny);
    println!("Grid spacing: {:?}", grid.spacing());

    // Set boundary conditions for heat diffusion
    // Left wall: T = 100°C (hot)
    // Right wall: T = 0°C (cold)
    // Top and bottom walls: insulated (∂T/∂n = 0, approximated as T = neighbor)

    // Set boundaries for all edges at once
    grid.set_edge_boundaries(BoundaryType::Wall);

    // Set boundary values using iterator patterns for zero-copy operations
    let boundary_values: HashMap<(usize, usize), f64> = [
        // Hot left wall - using iterator combinators
        (0..grid.ny).map(|j| ((0, j), 100.0)).collect::<Vec<_>>(),
        // Cold right wall - functional approach
        (0..grid.ny)
            .map(|j| ((grid.nx - 1, j), 0.0))
            .collect::<Vec<_>>(),
        // Insulated top and bottom walls (simplified as fixed temperature)
        // Note: This demonstrates iterator chaining for boundary condition setup
        (1..grid.nx - 1)
            .flat_map(|i| {
                [
                    ((i, 0), 50.0),           // Bottom wall (simplified)
                    ((i, grid.ny - 1), 50.0), // Top wall (simplified)
                ]
            })
            .collect::<Vec<_>>(),
    ]
    .into_iter()
    .flatten()
    .collect();

    // No source term for steady-state heat conduction
    let source = HashMap::new();

    // Configure the Poisson solver using builder pattern
    let base = SolverConfig::<f64>::builder()
        .tolerance(1e-8)
        .max_iterations(2000)
        .relaxation_factor(1.0)
        .verbose(true)
        .build();

    let config = FdmConfig { base };

    // Solve the heat equation
    println!("\nSolving 2D heat equation...");
    let solver = PoissonSolver::new(config);
    let solution = solver.solve(&grid, &source, &boundary_values)?;

    // Analyze results
    println!("\n=== Results Analysis ===");

    // Find temperature statistics using iterator combinators for zero-copy analysis

    let temperatures: Vec<f64> = grid
        .iter()
        .filter_map(|(i, j)| solution.get(&(i, j)).copied())
        .collect();

    let (min_temp, max_temp, avg_temp) = temperatures.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY, 0.0),
        |(min, max, sum), &temp| (min.min(temp), max.max(temp), sum + temp),
    );

    let avg_temp = avg_temp / temperatures.len() as f64;

    println!("Range: {:.2}°C to {:.2}°C", min_temp, max_temp);
    println!("Average temperature: {:.2}°C", avg_temp);

    // Print temperature distribution along centerline using iterator patterns
    println!("\n=== Horizontal centerline distribution ===");
    let center_j = grid.ny / 2;
    println!("Position (x)\tTemperature (°C)");

    let (x_min, x_max, _, _) = grid.bounds;
    let dx = (x_max - x_min) / grid.nx as f64;
    (0..grid.nx)
        .filter_map(|i| {
            solution.get(&(i, center_j)).map(|&temp| {
                let x = x_min + (i as f64 + 0.5) * dx;
                (x, temp)
            })
        })
        .for_each(|(x, temp)| println!("{:.3}\t\t{:.2}", x, temp));

    // Calculate heat flux at boundaries (approximate)
    println!("\n=== Heat Flux Analysis ===");
    let (dx, _dy) = grid.spacing();

    // Heat flux at left boundary using iterator combinators for zero-copy calculation
    let total_heat_flux_left: f64 = (1..grid.ny - 1)
        .filter_map(|j| {
            solution
                .get(&(0, j))
                .zip(solution.get(&(1, j)))
                .map(|(&t_boundary, &t_interior)| -(t_interior - t_boundary) / dx)
        })
        .sum();

    println!(
        "Total heat flux from left wall: {:.2} W/m²",
        total_heat_flux_left
    );

    // Demonstrate advection-diffusion solver with flow
    println!("\n=== Advection-Diffusion Example ===");

    // Create velocity field using iterator patterns for zero-copy field initialization
    let (velocity_x, velocity_y): (GridScalarField, GridScalarField) = grid
        .iter()
        .map(|(i, j)| [((i, j), 1.0), ((i, j), 0.0)]) // [x_velocity, y_velocity]
        .fold(
            (HashMap::new(), HashMap::new()),
            |(mut vx, mut vy), [x_entry, y_entry]| {
                vx.insert(x_entry.0, x_entry.1);
                vy.insert(y_entry.0, y_entry.1);
                (vx, vy)
            },
        );

    // Solve advection-diffusion equation
    let diffusivity = 0.1; // Thermal diffusivity

    let advdiff_base = SolverConfig::<f64>::builder()
        .tolerance(1e-6)
        .max_iterations(1000)
        .relaxation_factor(0.8)
        .verbose(false)
        .build();

    let advdiff_config = FdmConfig { base: advdiff_base };

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

    for i in (0..grid.nx).step_by(4) {
        if let (Some(&temp_diff), Some(&temp_advdiff)) = (
            solution.get(&(i, center_j)),
            advdiff_solution.get(&(i, center_j)),
        ) {
            let x = x_min + (i as f64 + 0.5) * dx;
            let difference = temp_advdiff - temp_diff;
            println!(
                "{:.3}\t\t{:.2}\t\t{:.2}\t\t{:.2}",
                x, temp_diff, temp_advdiff, difference
            );
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
