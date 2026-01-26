//! Comprehensive convergence monitoring example
//!
//! Demonstrates proper convergence detection and diagnostic techniques for CFD solvers.

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::momentum::{BoundarySetup, MomentumComponent, MomentumSolver};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_validation::convergence::{ConvergenceMonitor, ConvergenceStatus};

/// Analytical solution for Poiseuille flow
/// u(y) = (1/2μ) * (dp/dx) * y * (H - y)
fn poiseuille_analytical(y: f64, height: f64, pressure_gradient: f64, viscosity: f64) -> f64 {
    -0.5 / viscosity * pressure_gradient * y * (height - y)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================");
    println!("  Convergence Monitoring Demonstration");
    println!("=================================================\n");

    // Run multiple test cases
    test_simple_diffusion()?;
    test_poiseuille_with_monitoring()?;
    test_oscillating_solution()?;

    println!("\n=================================================");
    println!("  All convergence monitoring tests complete");
    println!("=================================================");

    Ok(())
}

/// Test Case 1: Simple diffusion equation (well-behaved convergence)
fn test_simple_diffusion() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n1. SIMPLE DIFFUSION - Well-Behaved Convergence");
    println!("   -------------------------------------------");

    let nx = 50;
    let ny = 50;
    let alpha = 0.1; // Thermal diffusivity
    let dx = 1.0 / (nx - 1) as f64;
    let dt = 0.2 * dx * dx / alpha; // CFL condition

    let mut temperature = vec![vec![0.0; ny]; nx];
    let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-4, 1000);

    // Initialize with hot center
    temperature[nx / 2][ny / 2] = 100.0;

    println!("   Grid: {}x{}", nx, ny);
    println!("   Time step: {:.6e}", dt);

    for iteration in 0..1000 {
        let mut temp_new = temperature.clone();
        let mut max_change: f64 = 0.0;

        // Apply explicit diffusion
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let laplacian = (temperature[i + 1][j]
                    + temperature[i - 1][j]
                    + temperature[i][j + 1]
                    + temperature[i][j - 1]
                    - 4.0 * temperature[i][j])
                    / (dx * dx);

                temp_new[i][j] = temperature[i][j] + dt * alpha * laplacian;
                max_change = max_change.max((temp_new[i][j] - temperature[i][j]).abs());
            }
        }

        temperature = temp_new;
        monitor.update(max_change);

        if iteration % 100 == 0 {
            println!(
                "   Iteration {}: max change = {:.2e}",
                iteration, max_change
            );
        }

        match monitor.check_status() {
            ConvergenceStatus::Converged {
                final_error,
                iterations,
                criterion,
            } => {
                println!("\n   ✓ Converged after {} iterations", iterations);
                println!("     Final error: {:.2e}", final_error);
                println!("     Criterion: {:?}", criterion);
                break;
            }
            ConvergenceStatus::Diverging { growth_rate, .. } => {
                println!(
                    "\n   ✗ Divergence detected! Growth rate: {:.2}",
                    growth_rate
                );
                break;
            }
            ConvergenceStatus::Stalled {
                stall_error,
                stall_iterations,
            } => {
                println!(
                    "\n   ⚠ Convergence stalled at iteration {}",
                    stall_iterations
                );
                println!("     Stall error: {:.2e}", stall_error);
                break;
            }
            _ => {}
        }
    }

    Ok(())
}

/// Test Case 2: Poiseuille flow with detailed monitoring
fn test_poiseuille_with_monitoring() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n2. POISEUILLE FLOW - Convergence Diagnostics");
    println!("   -----------------------------------------");

    let channel_height = 1.0;
    let channel_length = 4.0;
    let viscosity = 1e-3;
    let pressure_gradient = -1.0;

    let nx = 41;
    let ny = 21;
    let dt = 1e10; // Steady-state

    let grid = StructuredGrid2D::new(nx, ny, 0.0, channel_length, 0.0, channel_height)?;
    let mut solver = MomentumSolver::new(&grid);

    // Use the sophisticated boundary condition setup framework
    BoundarySetup::new()
        .south(BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        })
        .north(BoundaryCondition::Dirichlet {
            value: 0.0,
            component_values: None,
        })
        .west(BoundaryCondition::Neumann { gradient: 0.0 })
        .east(BoundaryCondition::Neumann { gradient: 0.0 })
        .apply(&mut solver);

    // Validate boundary conditions
    solver.validate_boundary_conditions()?;

    let mut fields = SimulationFields::new(nx, ny);

    // Set pressure gradient
    let dx = channel_length / (nx - 1) as f64;
    for i in 0..nx {
        for j in 0..ny {
            let x = i as f64 * dx;
            fields.p.set(i, j, pressure_gradient * x);
        }
    }

    let mut velocity_monitor = ConvergenceMonitor::<f64>::new(2e-4, 1e-3, 10000);
    let mut residual_monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 10000);

    println!("   Grid: {}x{}", nx, ny);
    println!("   Viscosity: {:.2e}", viscosity);
    println!("   Pressure gradient: {:.2}", pressure_gradient);

    for iteration in 0..1000 {
        let u_old = fields.u.clone();

        solver.solve(MomentumComponent::U, &mut fields, dt)?;

        // Calculate velocity change and residual
        let mut max_change: f64 = 0.0;
        let mut residual: f64 = 0.0;
        let mut count = 0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let change = (fields.u.at(i, j) - u_old.at(i, j)).abs();
                max_change = max_change.max(change);
                residual += change * change;
                count += 1;
            }
        }

        residual = (residual / count as f64).sqrt();

        velocity_monitor.update(max_change);
        residual_monitor.update(residual);

        if iteration % 10 == 0 {
            let u_center = fields.u.at(nx / 2, ny / 2);
            println!(
                "   Iteration {}: max_change = {:.2e}, L2_residual = {:.2e}, u_center = {:.6}",
                iteration, max_change, residual, u_center
            );
        }

        // Check both monitors
        let vel_status = velocity_monitor.check_status();
        let res_status = residual_monitor.check_status();

        if vel_status.is_converged() && res_status.is_converged() {
            println!("\n   ✓ Converged after {} iterations", iteration);

            // Validate against analytical solution
            let dy = channel_height / (ny - 1) as f64;
            let mut max_error: f64 = 0.0;

            for j in 0..ny {
                let y = j as f64 * dy;
                let u_numerical = fields.u.at(nx / 2, j);
                let u_analytical =
                    poiseuille_analytical(y, channel_height, pressure_gradient, viscosity);
                let error = (u_numerical - u_analytical).abs();
                max_error = max_error.max(error);
            }

            println!("     Max error vs analytical: {:.2e}", max_error);

            if max_error < 1.0 {
                println!("     ✓ Solution within acceptable error");
            } else {
                println!("     ⚠ Solution error exceeds tolerance");
            }

            break;
        }

        if matches!(vel_status, ConvergenceStatus::Diverging { .. })
            || matches!(res_status, ConvergenceStatus::Diverging { .. })
        {
            println!("\n   ✗ Divergence detected!");
            println!("     Velocity status: {:?}", vel_status);
            println!("     Residual status: {:?}", res_status);
            break;
        }
    }

    Ok(())
}

/// Test Case 3: Oscillating solution (stall detection)
fn test_oscillating_solution() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n3. OSCILLATING SOLUTION - Stall Detection");
    println!("   --------------------------------------");

    let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-4, 100);

    // Simulate oscillating convergence
    let mut error = 1.0;
    for iteration in 0..50 {
        // Add oscillation that prevents convergence
        let oscillation = 0.1 * (iteration as f64 * 0.5).sin();
        let current_error = error * (1.0 + oscillation);

        monitor.update(current_error);
        error *= 0.98; // Slow reduction

        if iteration % 10 == 0 {
            println!("   Iteration {}: error = {:.2e}", iteration, current_error);
        }

        match monitor.check_status() {
            ConvergenceStatus::Stalled {
                stall_error,
                stall_iterations,
            } => {
                println!("\n   ⚠ Stall detected at iteration {}", stall_iterations);
                println!("     Stall error: {:.2e}", stall_error);
                println!("     Recommendation: Reduce time step or adjust relaxation");
                break;
            }
            ConvergenceStatus::Converged { .. } => {
                println!("\n   ✓ Converged despite oscillations");
                break;
            }
            _ => {}
        }
    }

    Ok(())
}
