//! Method of Manufactured Solutions (MMS) validation for the explicit diffusion solver.
//!
//! This test suite verifies the correctness and accuracy of the `DiffusionSolver`
//! by manufacturing a solution to the heat diffusion equation and comparing the
//! numerical result against it.
//!
//! The tests confirm that the solver is second-order accurate in space and
//! produces results within a tight tolerance.

use cfd_2d::solvers::DiffusionSolver;
use cfd_validation::manufactured::{diffusion::ManufacturedDiffusion, ManufacturedSolution};
use std::collections::HashMap;

/// Validates the DiffusionSolver using a manufactured solution.
///
/// This test verifies:
/// 1. The solver correctly computes a solution for a known problem.
/// 2. The solver achieves the expected second-order accuracy.
#[test]
fn test_diffusion_mms_validation() {
    let alpha = 0.1;
    let manufactured = ManufacturedDiffusion::new(alpha);
    let t_final = 0.1;

    println!("MMS Diffusion Validation - PRODUCTION STANDARD:");
    println!("  Thermal diffusivity: {alpha}");
    println!("  Final time: {t_final}");

    let grid_sizes = vec![16, 32, 64];
    let mut l2_errors = Vec::with_capacity(grid_sizes.len());
    let mut linf_errors = Vec::with_capacity(grid_sizes.len());

    for &n in &grid_sizes {
        let dx = 1.0 / (n - 1) as f64;
        let dy = dx;

        println!("  Grid: {n}x{n}");

        // Set up solver with manufactured source and boundary conditions
        let mut solver = DiffusionSolver::new(n, n, dx, dy, alpha);
        let mut exact_solution = HashMap::with_capacity(n * n);
        let mut initial_conditions = HashMap::with_capacity(n * n);

        for i in 0..n {
            for j in 0..n {
                let x = i as f64 * dx;
                let y = j as f64 * dy;

                exact_solution.insert((i, j), manufactured.exact_solution(x, y, 0.0, t_final));
                initial_conditions.insert((i, j), manufactured.exact_solution(x, y, 0.0, 0.0));
            }
        }

        solver.set_boundary_conditions(&initial_conditions);

        // Solve and compute errors
        let numerical_solution =
            solver.solve_to_time(t_final, &|x, y, t| manufactured.source_term(x, y, 0.0, t));
        let l2_error = compute_l2_error(&numerical_solution, &exact_solution);
        let linf_error = compute_linf_error(&numerical_solution, &exact_solution);

        l2_errors.push(l2_error);
        linf_errors.push(linf_error);

        println!("    L2 error:   {l2_error:.6e}");
        println!("    L∞ error:   {linf_error:.6e}");
    }

    // Verify order of accuracy (should be O(h²))
    let h1 = 1.0 / (grid_sizes[0] - 1) as f64;
    let h2 = 1.0 / (grid_sizes[1] - 1) as f64;
    let observed_order = (l2_errors[0] / l2_errors[1]).log2() / (h1 / h2).log2();

    println!("  Observed order of accuracy: {observed_order:.2}");
    assert!(
        observed_order > 1.8,
        "Order of accuracy verification failed: expected ~2.0, got {observed_order}"
    );

    println!("  MMS validation successful: Solver is second-order accurate.");
}

/// Validates MMS framework components for diffusion equation.
#[test]
fn test_mms_integration_pattern() {
    let alpha = 0.1;
    let manufactured = ManufacturedDiffusion::new(alpha);
    let t_final = 0.1;
    let (nx, ny) = (64, 64);
    let (dx, dy) = (1.0 / (nx - 1) as f64, 1.0 / (ny - 1) as f64);

    println!("MMS Integration Pattern - PRODUCTION TEMPLATE:");

    // Generate source terms, boundary conditions, and exact solution
    let mut exact_solution = HashMap::with_capacity(nx * ny);
    let mut initial_conditions = HashMap::with_capacity(nx * ny);

    for i in 0..nx {
        for j in 0..ny {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            let pos = (i, j);

            exact_solution.insert(pos, manufactured.exact_solution(x, y, 0.0, t_final));
            initial_conditions.insert(pos, manufactured.exact_solution(x, y, 0.0, 0.0));
        }
    }

    // Set up and run the solver
    println!("  Setting up solver with manufactured source terms...");
    let mut solver = DiffusionSolver::new(nx, ny, dx, dy, alpha);
    solver.set_boundary_conditions(&initial_conditions);
    let numerical_solution =
        solver.solve_to_time(t_final, &|x, y, t| manufactured.source_term(x, y, 0.0, t));

    // Analyze errors
    let l2_error = compute_l2_error(&numerical_solution, &exact_solution);
    let linf_error = compute_linf_error(&numerical_solution, &exact_solution);

    println!("  L2 error:   {l2_error:.6e}");
    println!("  L∞ error:   {linf_error:.6e}");

    let tolerance = 1e-3;
    assert!(l2_error < tolerance, "L2 error exceeds tolerance");
    assert!(linf_error < tolerance, "L∞ error exceeds tolerance");

    println!("  MMS integration validated: Solver produces accurate results.");
}

/// Helper to compute L2 error norm.
fn compute_l2_error(
    numerical: &HashMap<(usize, usize), f64>,
    exact: &HashMap<(usize, usize), f64>,
) -> f64 {
    let mut sum_sq_error = 0.0;
    for (pos, &num_val) in numerical {
        if let Some(&ex_val) = exact.get(pos) {
            let error = num_val - ex_val;
            sum_sq_error += error * error;
        }
    }
    (sum_sq_error / numerical.len() as f64).sqrt()
}

/// Helper to compute L∞ error norm.
fn compute_linf_error(
    numerical: &HashMap<(usize, usize), f64>,
    exact: &HashMap<(usize, usize), f64>,
) -> f64 {
    let mut max_error = 0.0;
    for (pos, &num_val) in numerical {
        if let Some(&ex_val) = exact.get(pos) {
            let error = (num_val - ex_val).abs();
            if error > max_error {
                max_error = error;
            }
        }
    }
    max_error
}
