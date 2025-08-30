//! Comprehensive validation suite using manufactured solutions and analytical benchmarks
//!
//! This example demonstrates validation of numerical methods against known solutions.

use cfd_suite::d2::fields::Field2D;
use cfd_suite::d2::grid::StructuredGrid2D;
use cfd_suite::math::differentiation::FiniteDifference;
use cfd_suite::validation::convergence::{GridConvergenceIndex, RichardsonExtrapolation};
use cfd_suite::validation::error_metrics::ErrorMetric;
use nalgebra::DMatrix;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("==============================================");
    println!("    CFD Validation Suite - Known Solutions    ");
    println!("==============================================\n");

    // Run different validation tests
    validate_diffusion()?;
    validate_advection()?;
    validate_taylor_green()?;
    validate_poiseuille()?;
    validate_grid_convergence()?;

    println!("\n==============================================");
    println!("    All Validation Tests Passed!             ");
    println!("==============================================");

    Ok(())
}

/// Validate 2D diffusion solver using manufactured solution
fn validate_diffusion() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n1. DIFFUSION EQUATION VALIDATION");
    println!("   Method: Manufactured Solution");
    println!("   --------------------------------");

    let alpha = 0.1; // Thermal diffusivity
    let solution = ManufacturedDiffusion::new(alpha);

    // Create grids of different resolutions for convergence study
    let resolutions = vec![16, 32, 64];
    let mut errors = Vec::new();

    for n in &resolutions {
        let grid = StructuredGrid2D::<f64>::new(*n, *n, 0.0, 1.0, 0.0, 1.0)?;
        let dx = 1.0 / (*n as f64);
        let dt = 0.5 * dx * dx / alpha; // CFL condition for diffusion

        // Initialize field with manufactured initial condition
        let mut field = Field2D::new(*n, *n, 0.0);
        for i in 0..*n {
            for j in 0..*n {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                field.set(i, j, solution.initial_condition(x, y, 0.0));
            }
        }

        // Evolve for a short time using explicit finite difference
        let t_final = 0.1;
        let n_steps = (t_final / dt).ceil() as usize;

        for _ in 0..n_steps {
            let mut field_new = field.clone();

            // Apply finite difference diffusion operator
            for i in 1..(*n - 1) {
                for j in 1..(*n - 1) {
                    let laplacian = (field.at(i + 1, j)
                        + field.at(i - 1, j)
                        + field.at(i, j + 1)
                        + field.at(i, j - 1)
                        - 4.0 * field.at(i, j))
                        / (dx * dx);

                    let x = i as f64 * dx;
                    let y = j as f64 * dx;
                    let source = solution.source_term(x, y, 0.0, t_final);

                    field_new.set(i, j, field.at(i, j) + dt * (alpha * laplacian + source));
                }
            }

            field = field_new;
        }

        // Calculate error
        let mut l2_error = 0.0;
        let mut max_error = 0.0;

        for i in 0..*n {
            for j in 0..*n {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                let exact = solution.exact_solution(x, y, 0.0, t_final);
                let numerical = field.at(i, j);
                let error = (numerical - exact).abs();

                l2_error += error * error;
                max_error = max_error.max(error);
            }
        }

        l2_error = (l2_error / (*n * n) as f64).sqrt();
        errors.push((dx, l2_error, max_error));

        println!("   Grid: {}x{}", n, n);
        println!("   L2 Error: {:.6e}", l2_error);
        println!("   Max Error: {:.6e}", max_error);
    }

    // Check convergence rate
    if errors.len() >= 2 {
        let rate = ((errors[1].1 / errors[0].1).ln()) / ((errors[1].0 / errors[0].0).ln());
        println!("   Convergence Rate: {:.2}", rate);
        println!("   Expected: ~2.0 (second-order)");

        // Verify second-order convergence
        assert!(
            rate > 1.8 && rate < 2.2,
            "Convergence rate {:.2} outside expected range [1.8, 2.2]",
            rate
        );
    }

    Ok(())
}

/// Validate advection solver using manufactured solution
fn validate_advection() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n2. ADVECTION EQUATION VALIDATION");
    println!("   Method: Manufactured Solution");
    println!("   --------------------------------");

    let vx = 1.0;
    let vy = 0.5;
    let solution = ManufacturedAdvection::new(vx, vy);

    let n = 64;
    let grid = StructuredGrid2D::<f64>::new(n, n, 0.0, 1.0, 0.0, 1.0)?;
    let dx = 1.0 / (n as f64);
    let dt = 0.5 * dx / (vx.abs().max(vy.abs())); // CFL condition

    // Initialize field
    let mut field = Field2D::new(n, n, 0.0);
    for i in 0..n {
        for j in 0..n {
            let x = i as f64 * dx;
            let y = j as f64 * dx;
            field.set(i, j, solution.initial_condition(x, y, 0.0));
        }
    }

    // Evolve using upwind scheme
    let t_final = 0.1;
    let n_steps = (t_final / dt).ceil() as usize;
    let mut t = 0.0;

    for _ in 0..n_steps {
        let mut field_new = field.clone();

        for i in 1..(n - 1) {
            for j in 1..(n - 1) {
                // Upwind differencing
                let dudx = if vx > 0.0 {
                    (field.at(i, j) - field.at(i - 1, j)) / dx
                } else {
                    (field.at(i + 1, j) - field.at(i, j)) / dx
                };

                let dudy = if vy > 0.0 {
                    (field.at(i, j) - field.at(i, j - 1)) / dx
                } else {
                    (field.at(i, j + 1) - field.at(i, j)) / dx
                };

                field_new.set(i, j, field.at(i, j) - dt * (vx * dudx + vy * dudy));
            }
        }

        field = field_new;
        t += dt;
    }

    // Calculate error
    let mut l2_error = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = i as f64 * dx;
            let y = j as f64 * dx;
            let exact = solution.exact_solution(x, y, 0.0, t_final);
            let numerical = field.at(i, j);
            l2_error += (numerical - exact).powi(2);
        }
    }
    l2_error = (l2_error / (n * n) as f64).sqrt();

    println!("   Grid: {}x{}", n, n);
    println!("   L2 Error: {:.6e}", l2_error);
    println!("   Time Steps: {}", n_steps);

    // Verify accuracy
    assert!(
        l2_error < 0.01,
        "L2 error {:.6e} exceeds tolerance",
        l2_error
    );

    Ok(())
}

/// Validate using Taylor-Green vortex decay
fn validate_taylor_green() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n3. TAYLOR-GREEN VORTEX VALIDATION");
    println!("   Method: Analytical Solution");
    println!("   --------------------------------");

    let nu = 0.01;
    let tg = TaylorGreenVortex::new(nu, 1.0, 1.0);

    // Check kinetic energy decay
    let times = vec![0.0, 0.5, 1.0];
    println!("   Kinetic Energy Decay:");

    for &t in &times {
        let ke_exact = 0.25 * (-4.0 * nu * PI * PI * t).exp();
        println!("   t = {:.1}: KE = {:.6}", t, ke_exact);
    }

    // Verify vorticity evolution at a point
    let x = 0.5;
    let y = 0.5;
    let t = 1.0;

    let velocity = tg.evaluate(x, y, 0.0, t);
    println!(
        "\n   Velocity at ({}, {}, t={}): [{:.6}, {:.6}, {:.6}]",
        x, y, t, velocity.x, velocity.y, velocity.z
    );

    // Check that solution decays as expected
    let v0 = tg.evaluate(x, y, 0.0, 0.0);
    let v1 = tg.evaluate(x, y, 0.0, 1.0);
    let decay_ratio = v1.norm() / v0.norm();
    let expected_decay = (-2.0 * nu * PI * PI).exp();

    println!("   Decay Ratio: {:.6}", decay_ratio);
    println!("   Expected: {:.6}", expected_decay);

    assert!(
        (decay_ratio - expected_decay).abs() < 1e-6,
        "Decay ratio mismatch"
    );

    Ok(())
}

/// Validate Poiseuille flow between parallel plates
fn validate_poiseuille() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n4. POISEUILLE FLOW VALIDATION");
    println!("   Method: Analytical Solution");
    println!("   --------------------------------");

    let h = 1.0; // Channel half-height
    let dp_dx = -1.0; // Pressure gradient
    let mu = 0.01; // Dynamic viscosity

    let poiseuille = PoiseuilleFlow::new(h, dp_dx, mu);

    // Sample velocity profile at different heights
    let n_points = 11;
    println!("   Velocity Profile:");

    for i in 0..n_points {
        let y = -h + 2.0 * h * (i as f64) / ((n_points - 1) as f64);
        let velocity = poiseuille.evaluate(0.0, y, 0.0, 0.0);
        let u_exact = -dp_dx / (2.0 * mu) * (h * h - y * y);

        println!(
            "   y = {:5.2}: u = {:.6} (exact: {:.6})",
            y, velocity.x, u_exact
        );

        assert!(
            (velocity.x - u_exact).abs() < 1e-10,
            "Velocity mismatch at y = {}",
            y
        );
    }

    // Check maximum velocity at centerline
    let u_max = poiseuille.evaluate(0.0, 0.0, 0.0, 0.0).x;
    let u_max_exact = -dp_dx * h * h / (2.0 * mu);

    println!("\n   Maximum Velocity: {:.6}", u_max);
    println!("   Expected: {:.6}", u_max_exact);

    assert!(
        (u_max - u_max_exact).abs() < 1e-10,
        "Maximum velocity mismatch"
    );

    Ok(())
}

/// Validate grid convergence using Richardson extrapolation
fn validate_grid_convergence() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n5. GRID CONVERGENCE STUDY");
    println!("   Method: Richardson Extrapolation");
    println!("   --------------------------------");

    // Use a simple test function: sin(2πx) * cos(2πy)
    let test_fn = |x: f64, y: f64| (2.0 * PI * x).sin() * (2.0 * PI * y).cos();

    // Compute Laplacian on different grid sizes
    let grids = vec![16, 32, 64];
    let mut results = Vec::new();

    for n in &grids {
        let dx = 1.0 / (*n as f64);

        // Evaluate Laplacian at center point (0.5, 0.5)
        let x = 0.5;
        let y = 0.5;

        // Use 5-point stencil
        let f_center = test_fn(x, y);
        let f_right = test_fn(x + dx, y);
        let f_left = test_fn(x - dx, y);
        let f_top = test_fn(x, y + dx);
        let f_bottom = test_fn(x, y - dx);

        let laplacian = (f_right + f_left + f_top + f_bottom - 4.0 * f_center) / (dx * dx);
        results.push((*n, laplacian));

        println!("   Grid {}x{}: Laplacian = {:.6}", n, n, laplacian);
    }

    // Apply Richardson extrapolation
    if results.len() >= 3 {
        let f1 = results[0].1;
        let f2 = results[1].1;
        let f3 = results[2].1;

        let r: f64 = 2.0; // Grid refinement ratio
        let p = ((f3 - f2) / (f2 - f1)).ln() / r.ln();

        println!("\n   Order of Accuracy: {:.2}", p);
        println!("   Expected: ~2.0 (second-order)");

        // Richardson extrapolated value
        let f_exact_estimate = f3 + (f3 - f2) / (r.powf(p) - 1.0);
        println!("   Richardson Extrapolation: {:.6}", f_exact_estimate);

        // True analytical value: -8π² * sin(2πx) * cos(2πy)
        let exact = -8.0 * PI * PI * test_fn(0.5, 0.5);
        println!("   Analytical Value: {:.6}", exact);

        let error = (f_exact_estimate - exact).abs();
        println!("   Extrapolation Error: {:.6e}", error);

        assert!(
            p > 1.9 && p < 2.1,
            "Order of accuracy {:.2} outside expected range",
            p
        );
    }

    Ok(())
}
