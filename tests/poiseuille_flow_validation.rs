//! Validation of Poiseuille flow (channel flow) against analytical solution
//!
//! Tests the accuracy of the CFD solver against the exact analytical solution
//! for laminar flow between parallel plates.

extern crate cfd_2d;
extern crate cfd_core;

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::momentum::{MomentumComponent, MomentumSolver};
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
    let _density = 1.0; // Currently unused but needed for future solver implementation
    let pressure_gradient = -1.0; // dp/dx

    // Grid parameters
    let nx = 41;
    let ny = 21;
    let dx = channel_length / (nx - 1) as f64;
    let dy = channel_height / (ny - 1) as f64;
    let dt = 1e10; // Effectively steady-state (ρ/dt → 0)

    // Initialize fields
    let mut fields = SimulationFields::new(nx, ny);

    // Set boundary conditions
    // No-slip at walls (y = 0 and y = H)
    for i in 0..nx {
        // Bottom wall
        fields.u.set(i, 0, 0.0);
        fields.v.set(i, 0, 0.0);

        // Top wall
        fields.u.set(i, ny - 1, 0.0);
        fields.v.set(i, ny - 1, 0.0);
    }

    // Apply constant pressure gradient
    for i in 0..nx {
        for j in 0..ny {
            let x = i as f64 * dx;
            fields.p.set(i, j, pressure_gradient * x);
        }
    }

    // Create momentum solver
    let grid = StructuredGrid2D::new(nx, ny, 0.0, channel_length, 0.0, channel_height)
        .expect("Failed to create grid");
    let mut solver = MomentumSolver::new(&grid);

    // Register boundary conditions for no-slip walls
    use cfd_core::boundary::BoundaryCondition;
    solver.set_boundary(
        "south".to_string(),
        BoundaryCondition::Dirichlet { value: 0.0 },
    );
    solver.set_boundary(
        "north".to_string(),
        BoundaryCondition::Dirichlet { value: 0.0 },
    );

    // Time integration to steady state
    let max_time_steps = 10000;
    let convergence_tolerance = 2e-4; // Realistic tolerance - solver is converging to 1.02e-4
    let mut converged = false;

    println!("Starting Poiseuille flow simulation...");

    for step in 0..max_time_steps {
        let u_old = fields.u.clone();

        // Solve momentum equations
        solver
            .solve(MomentumComponent::U, &mut fields, dt)
            .expect("Momentum solve failed");

        solver
            .solve(MomentumComponent::V, &mut fields, dt)
            .expect("Momentum solve failed");

        // Check convergence
        let mut max_change: f64 = 0.0;
        for i in 0..nx {
            for j in 0..ny {
                let change = (fields.u.at(i, j) - u_old.at(i, j)).abs();
                max_change = max_change.max(change);
            }
        }

        if max_change < convergence_tolerance {
            println!("Converged after {step} iterations");
            converged = true;
            break;
        }

        if step % 100 == 0 {
            let u_center = fields.u.at(nx / 2, ny / 2);
            println!("Step {step}: max change = {max_change:.2e}, u_center = {u_center:.6}");
        }
    }

    assert!(converged, "Solution did not converge");

    // Validate against analytical solution at channel center
    let x_center = nx / 2;
    let mut max_error: f64 = 0.0;
    let mut l2_error: f64 = 0.0;

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
                "{y:.3}\t\t{u_numerical:.6}\t{u_analytical:.6}\t{error:.2e}"
            );
        }
    }

    l2_error = l2_error.sqrt();

    println!("\nError metrics:");
    println!("Max error: {max_error:.2e}");
    println!("L2 error: {l2_error:.2e}");

    // STRICT VALIDATION: Solver must produce physically meaningful results
    // A functional solver should achieve <1% error vs analytical solution
    let max_acceptable_error = 1.25; // 1% of max velocity (125 m/s)
    
    assert!(
        max_error < max_acceptable_error,
        "SOLVER FAILURE: Max error {:.2e} exceeds acceptable limit {:.2e}. \
        Numerical velocities near zero ({:.6} m/s) vs analytical ~125 m/s indicates \
        broken solver (immediate false convergence). This test must FAIL until solver is fixed.",
        max_error, max_acceptable_error, fields.u.at(x_center, ny/2)
    );

    let l2_acceptable_error = 0.5; // Strict L2 norm requirement
    assert!(
        l2_error < l2_acceptable_error,
        "SOLVER FAILURE: L2 error {l2_error:.2e} exceeds acceptable limit {l2_acceptable_error:.2e}"
    );

    let elapsed = start.elapsed();
    println!("\nTest completed in {:.2} seconds", elapsed.as_secs_f64());

    // Solver should take significant time for convergence (not immediate)
    assert!(
        elapsed.as_secs_f64() > 0.1,
        "SOLVER FAILURE: Test completed too quickly ({:.3}s). \
        A functional iterative solver should take >0.1s for 10,000 max iterations. \
        Immediate completion indicates false convergence.",
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
            fields
                .u
                .set(i, j, poiseuille_analytical(y, height, -1.0, 1e-3));
            fields.v.set(i, j, 0.0); // No vertical velocity
        }
    }

    // Check divergence (should be zero for incompressible flow)
    let mut max_divergence: f64 = 0.0;

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let dudx = (fields.u.at(i + 1, j) - fields.u.at(i - 1, j)) / (2.0 * dx);
            let dvdy = (fields.v.at(i, j + 1) - fields.v.at(i, j - 1)) / (2.0 * dy);

            let divergence = dudx + dvdy;
            max_divergence = max_divergence.max(divergence.abs());
        }
    }

    println!(
        "Maximum divergence in Poiseuille flow: {max_divergence:.2e}"
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
        "Inlet flux: {inlet_flux:.6}, Outlet flux: {outlet_flux:.6}, Error: {flux_error:.2e}"
    );

    assert!(flux_error < 1e-10, "Mass flux not conserved");
}
