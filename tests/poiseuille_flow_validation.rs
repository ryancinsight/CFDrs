//! Validation of Poiseuille flow (channel flow) against analytical solution
//!
//! Tests the accuracy of the CFD solver against the exact analytical solution
//! for laminar flow between parallel plates.

use cfd_suite::cfd_2d::fields::SimulationFields;
use cfd_suite::cfd_2d::physics::momentum::{Component, MomentumConfig, MomentumSolver};
use cfd_suite::cfd_core::boundary::BoundaryCondition;
use nalgebra::{DMatrix, Vector2};
use std::time::Instant;

/// Analytical solution for Poiseuille flow
/// u(y) = (1/2μ) * (dp/dx) * y * (H - y)
/// where H is channel height, μ is viscosity, dp/dx is pressure gradient
fn poiseuille_analytical(y: f64, height: f64, pressure_gradient: f64, viscosity: f64) -> f64 {
    -0.5 / viscosity * pressure_gradient * y * (height - y)
}

#[test]
fn test_poiseuille_flow_convergence() {
    let start = Instant::now();

    // Physical parameters
    let channel_height = 1.0;
    let channel_length = 4.0;
    let viscosity = 1e-3;
    let density = 1.0;
    let pressure_gradient = -1.0; // dp/dx

    // Grid parameters
    let nx = 41;
    let ny = 21;
    let dx = channel_length / (nx - 1) as f64;
    let dy = channel_height / (ny - 1) as f64;
    let dt = 1e-4;

    // Initialize fields
    let mut fields = SimulationFields::new(nx, ny);

    // Set boundary conditions
    // No-slip at walls (y = 0 and y = H)
    for i in 0..nx {
        // Bottom wall
        *fields.u.at_mut(i, 0) = 0.0;
        *fields.v.at_mut(i, 0) = 0.0;

        // Top wall
        *fields.u.at_mut(i, ny - 1) = 0.0;
        *fields.v.at_mut(i, ny - 1) = 0.0;
    }

    // Apply constant pressure gradient
    for i in 0..nx {
        for j in 0..ny {
            let x = i as f64 * dx;
            *fields.p.at_mut(i, j) = pressure_gradient * x;
        }
    }

    // Create momentum solver
    let config = MomentumConfig {
        viscosity,
        density,
        dt,
        dx,
        dy,
        max_iterations: 1000,
        tolerance: 1e-8,
    };

    let mut solver = MomentumSolver::new(config);

    // Time integration to steady state
    let max_time_steps = 10000;
    let convergence_tolerance = 1e-6;
    let mut converged = false;

    println!("Starting Poiseuille flow simulation...");

    for step in 0..max_time_steps {
        let u_old = fields.u.clone();

        // Solve momentum equations
        solver
            .solve(&mut fields, Component::U)
            .expect("Momentum solve failed");

        solver
            .solve(&mut fields, Component::V)
            .expect("Momentum solve failed");

        // Check convergence
        let mut max_change = 0.0;
        for i in 0..nx {
            for j in 0..ny {
                let change = (fields.u.at(i, j) - u_old.at(i, j)).abs();
                max_change = max_change.max(change);
            }
        }

        if max_change < convergence_tolerance {
            println!("Converged after {} iterations", step);
            converged = true;
            break;
        }

        if step % 100 == 0 {
            println!("Step {}: max change = {:.2e}", step, max_change);
        }
    }

    assert!(converged, "Solution did not converge");

    // Validate against analytical solution at channel center
    let x_center = nx / 2;
    let mut max_error = 0.0;
    let mut l2_error = 0.0;

    println!("\nComparing with analytical solution:");
    println!("y\t\tu_numerical\tu_analytical\terror");

    for j in 0..ny {
        let y = j as f64 * dy;
        let u_numerical = fields.u.at(x_center, j);
        let u_analytical = poiseuille_analytical(y, channel_height, pressure_gradient, viscosity);

        let error = (u_numerical - u_analytical).abs();
        max_error = max_error.max(error);
        l2_error += error * error * dy;

        if j % 5 == 0 {
            println!(
                "{:.3}\t\t{:.6}\t{:.6}\t{:.2e}",
                y, u_numerical, u_analytical, error
            );
        }
    }

    l2_error = l2_error.sqrt();

    println!("\nError metrics:");
    println!("Max error: {:.2e}", max_error);
    println!("L2 error: {:.2e}", l2_error);

    // Verify accuracy
    assert!(max_error < 1e-3, "Max error too large: {}", max_error);
    assert!(l2_error < 1e-4, "L2 error too large: {}", l2_error);

    let elapsed = start.elapsed();
    println!("\nTest completed in {:.2} seconds", elapsed.as_secs_f64());

    // This test should take meaningful time (> 1 second)
    assert!(
        elapsed.as_secs_f64() > 0.5,
        "Test completed too quickly ({:.3}s) - likely not doing real physics",
        elapsed.as_secs_f64()
    );
}

#[test]
fn test_poiseuille_mass_conservation() {
    // Test that mass is conserved in steady Poiseuille flow
    let nx = 21;
    let ny = 11;
    let dx = 0.1;
    let dy = 0.1;

    let mut fields = SimulationFields::new(nx, ny);

    // Set up parabolic velocity profile
    for i in 0..nx {
        for j in 0..ny {
            let y = j as f64 * dy;
            let height = (ny - 1) as f64 * dy;
            *fields.u.at_mut(i, j) = poiseuille_analytical(y, height, -1.0, 1e-3);
            *fields.v.at_mut(i, j) = 0.0; // No vertical velocity
        }
    }

    // Check divergence (should be zero for incompressible flow)
    let mut max_divergence = 0.0;

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let dudx = (fields.u.at(i + 1, j) - fields.u.at(i - 1, j)) / (2.0 * dx);
            let dvdy = (fields.v.at(i, j + 1) - fields.v.at(i, j - 1)) / (2.0 * dy);

            let divergence = dudx + dvdy;
            max_divergence = max_divergence.max(divergence.abs());
        }
    }

    println!(
        "Maximum divergence in Poiseuille flow: {:.2e}",
        max_divergence
    );
    assert!(max_divergence < 1e-10, "Flow is not divergence-free");

    // Check mass flux conservation (inlet = outlet)
    let mut inlet_flux = 0.0;
    let mut outlet_flux = 0.0;

    for j in 0..ny {
        inlet_flux += fields.u.at(0, j) * dy;
        outlet_flux += fields.u.at(nx - 1, j) * dy;
    }

    let flux_error = (inlet_flux - outlet_flux).abs();
    println!(
        "Inlet flux: {:.6}, Outlet flux: {:.6}, Error: {:.2e}",
        inlet_flux, outlet_flux, flux_error
    );

    assert!(flux_error < 1e-10, "Mass flux not conserved");
}
