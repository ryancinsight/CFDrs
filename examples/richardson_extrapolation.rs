//! Richardson Extrapolation convergence study example
//!
//! Demonstrates systematic grid convergence analysis using Richardson extrapolation
//! following Roache (1998) methodology and ASME V&V 20-2009 standards.

use cfd_2d::fields::Field2D;
use cfd_2d::grid::StructuredGrid2D;
use cfd_validation::convergence::{ConvergenceStudy, GridConvergenceIndex};
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=========================================================");
    println!("  Richardson Extrapolation - Grid Convergence Study");
    println!("=========================================================\n");

    println!("Following Roache (1998) and ASME V&V 20-2009, this example");
    println!("demonstrates systematic grid convergence analysis using");
    println!("Richardson extrapolation and Grid Convergence Index (GCI).\n");

    // Run convergence studies
    richardson_heat_equation()?;
    richardson_poisson_equation()?;
    gci_uncertainty_quantification()?;

    println!("\n=========================================================");
    println!("  Grid Convergence Study Complete");
    println!("=========================================================");

    Ok(())
}

/// Richardson extrapolation for 2D heat equation
///
/// Uses systematic grid refinement to estimate discretization error
fn richardson_heat_equation() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n1. HEAT EQUATION - Richardson Extrapolation");
    println!("   ----------------------------------------");

    let alpha = 0.1; // Thermal diffusivity
    let t_final = 0.01;

    // Three grid levels with refinement ratio r = 2
    let resolutions = vec![32, 64, 128];
    let mut grid_sizes = Vec::new();
    let mut temperatures_center = Vec::new();

    println!("   Problem: 2D heat diffusion with point source");
    println!("   Refinement ratio: r = 2");
    println!("   Time: t = {:.3}", t_final);

    for &n in &resolutions {
        let _grid = StructuredGrid2D::<f64>::new(n, n, 0.0, 1.0, 0.0, 1.0)?;
        let dx = 1.0 / (n - 1) as f64;
        grid_sizes.push(dx);

        // CFL condition
        let dt = 0.2 * dx * dx / alpha;
        let n_steps = (t_final / dt).ceil() as usize;

        // Initialize with point source at center
        let mut field = Field2D::new(n, n, 0.0);
        field.set(n / 2, n / 2, 100.0);

        // Time evolution
        for _ in 0..n_steps {
            let mut field_new = field.clone();

            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    let laplacian = (field.at(i + 1, j)
                        + field.at(i - 1, j)
                        + field.at(i, j + 1)
                        + field.at(i, j - 1)
                        - 4.0 * field.at(i, j))
                        / (dx * dx);

                    field_new.set(i, j, field.at(i, j) + dt * alpha * laplacian);
                }
            }

            field = field_new;
        }

        // Extract solution at center
        let t_center = field.at(n / 2, n / 2);
        temperatures_center.push(t_center);

        println!("   Grid {}×{}: T_center = {:.6}", n, n, t_center);
    }

    // Richardson extrapolation with r = 2
    let r: f64 = 2.0;
    let f_fine = temperatures_center[2]; // 128×128
    let f_medium = temperatures_center[1]; // 64×64
    let f_coarse = temperatures_center[0]; // 32×32

    // Estimate order of accuracy
    let epsilon_21 = f_medium - f_fine;
    let epsilon_32 = f_coarse - f_medium;
    let p = (epsilon_32 / epsilon_21).abs().log2();

    println!("\n   Richardson Analysis:");
    println!("     ε₂₁ = {:.6e} (64→128 change)", epsilon_21);
    println!("     ε₃₂ = {:.6e} (32→64 change)", epsilon_32);
    println!("     Observed order p = {:.2}", p);

    // Richardson extrapolation to infinite resolution
    let f_extrapolated = f_fine + epsilon_21 / (r.powf(p) - 1.0);
    println!("     T_∞ (extrapolated) = {:.6}", f_extrapolated);

    // Relative error on finest grid
    let relative_error = ((f_fine - f_extrapolated) / f_extrapolated).abs() * 100.0;
    println!("     Relative error (128×128): {:.3}%", relative_error);

    if relative_error < 1.0 {
        println!("     ✓ Solution within 1% of extrapolated value");
    } else {
        println!("     ⚠ Consider finer grids for better accuracy");
    }

    Ok(())
}

/// Richardson extrapolation for Poisson equation
///
/// ∇²u = f(x,y) with manufactured right-hand side
fn richardson_poisson_equation() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n2. POISSON EQUATION - Grid Convergence Study");
    println!("   -----------------------------------------");

    // Manufactured solution: u(x,y) = sin(πx)sin(πy)
    // Source term: f = -2π² sin(πx)sin(πy)

    let resolutions = vec![16, 32, 64, 128];
    let mut grid_sizes = Vec::new();
    let mut l2_errors = Vec::new();
    let mut max_errors = Vec::new();

    println!("   Manufactured solution: u = sin(πx)sin(πy)");
    println!("   Source: f = -2π²sin(πx)sin(πy)");

    for &n in &resolutions {
        let dx = 1.0 / (n - 1) as f64;
        grid_sizes.push(dx);

        // Solve Poisson equation using Jacobi iteration
        let max_iterations = 10000;
        let tolerance = 1e-10;

        let mut u = vec![vec![0.0; n]; n];

        for _iter in 0..max_iterations {
            let mut u_new = u.clone();
            let mut max_change: f64 = 0.0;

            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    let x = i as f64 * dx;
                    let y = j as f64 * dx;

                    let f = -2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();

                    let laplacian_u = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]
                        - 4.0 * u[i][j])
                        / (dx * dx);

                    u_new[i][j] = u[i][j] + 0.25 * (f - laplacian_u) * dx * dx;

                    max_change = max_change.max((u_new[i][j] - u[i][j]).abs());
                }
            }

            u = u_new;

            if max_change < tolerance {
                break;
            }
        }

        // Compute errors
        let mut l2_error = 0.0;
        let mut max_error: f64 = 0.0;
        let mut count = 0;

        let interior = n - 2;
        for (i, row) in u.iter().enumerate().skip(1).take(interior) {
            let x = i as f64 * dx;
            for (j, &u_ij) in row.iter().enumerate().skip(1).take(interior) {
                let y = j as f64 * dx;
                let exact = (PI * x).sin() * (PI * y).sin();
                let error = (u_ij - exact).abs();

                l2_error += error * error;
                max_error = max_error.max(error);
                count += 1;
            }
        }

        l2_error = (l2_error / count as f64).sqrt();
        l2_errors.push(l2_error);
        max_errors.push(max_error);

        println!(
            "   Grid {}×{}: L2={:.4e}, L∞={:.4e}",
            n, n, l2_error, max_error
        );
    }

    // Convergence analysis
    let study = ConvergenceStudy::new(grid_sizes.clone(), l2_errors.clone())?;

    println!("\n   Convergence Study:");
    println!("     Observed order (L2): {:.2}", study.convergence_rate);
    println!("     R²: {:.6}", study.r_squared);
    println!("     Asymptotic range: {}", study.is_asymptotic());

    // Predict error for even finer grid
    // NOTE: Using unwrap() is acceptable for example code demonstration.
    let h_predict = grid_sizes.last().unwrap() / 2.0;
    let error_predict = study.predict_error(h_predict);
    println!(
        "     Predicted error (h={:.4e}): {:.4e}",
        h_predict, error_predict
    );

    // Grid size for target accuracy
    let target_error = 1e-6;
    let h_target = study.grid_size_for_error(target_error)?;
    let n_target = (1.0 / h_target).ceil() as usize;
    println!(
        "     Grid for ε<{:.0e}: {}×{} (h={:.4e})",
        target_error, n_target, n_target, h_target
    );

    Ok(())
}

/// Grid Convergence Index (GCI) for uncertainty quantification
///
/// Following Roache (1998) methodology for reporting numerical uncertainty
fn gci_uncertainty_quantification() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n3. GCI UNCERTAINTY QUANTIFICATION");
    println!("   -------------------------------");

    // Use data from Poisson equation example
    let resolutions = vec![32, 64, 128];
    let refinement_ratio: f64 = 2.0;
    let observed_order: f64 = 2.0;

    // Recompute for demonstration
    let mut solutions = Vec::new();

    for &n in &resolutions {
        let dx = 1.0 / (n - 1) as f64;

        // Quick solve (reduced iterations for speed)
        let max_iterations = 1000;
        let mut u = vec![vec![0.0; n]; n];

        for _ in 0..max_iterations {
            let mut u_new = u.clone();

            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    let x = i as f64 * dx;
                    let y = j as f64 * dx;

                    let f = -2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
                    let laplacian_u = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]
                        - 4.0 * u[i][j])
                        / (dx * dx);

                    u_new[i][j] = u[i][j] + 0.25 * (f - laplacian_u) * dx * dx;
                }
            }

            u = u_new;
        }

        // Extract solution at center
        solutions.push(u[n / 2][n / 2]);
    }

    println!("   Grid solutions:");
    for (i, &sol) in solutions.iter().enumerate() {
        println!("     {}×{}: φ = {:.8}", resolutions[i], resolutions[i], sol);
    }

    // Compute GCI
    let gci = GridConvergenceIndex::<f64>::new(3, observed_order, refinement_ratio);

    let gci_fine = gci.compute_fine(solutions[2], solutions[1]);
    let gci_medium = gci.compute_fine(solutions[1], solutions[0]);

    println!("\n   GCI Analysis:");
    println!("     GCI_fine (128×128): {:.4}%", gci_fine * 100.0);
    println!("     GCI_medium (64×64): {:.4}%", gci_medium * 100.0);

    // Check asymptotic range
    let asymptotic = gci.is_asymptotic(gci_fine, gci_medium);
    println!("     Asymptotic range: {}", asymptotic);

    if asymptotic {
        println!("     ✓ Solutions in asymptotic convergence range");
    } else {
        println!("     ⚠ Not yet in asymptotic range, consider finer grids");
    }

    // Uncertainty bounds
    let (lower, upper) = gci.uncertainty_band(solutions[2], gci_fine);
    println!("\n   Uncertainty Band (128×128):");
    println!("     φ_exact ∈ [{:.8}, {:.8}]", lower, upper);
    println!("     Width: {:.4e}", upper - lower);

    Ok(())
}
