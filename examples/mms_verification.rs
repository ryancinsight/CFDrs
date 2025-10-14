//! Method of Manufactured Solutions (MMS) validation example
//!
//! Demonstrates code verification using MMS following Roache (1998) and Salari & Knupp (2000).
//! MMS enables systematic verification of discretization order by constructing exact solutions
//! and corresponding source terms.

use cfd_2d::fields::Field2D;
use cfd_2d::grid::StructuredGrid2D;
use cfd_validation::convergence::ConvergenceStudy;
use cfd_validation::manufactured::{ManufacturedAdvection, ManufacturedDiffusion, ManufacturedSolution};
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("====================================================");
    println!("  Method of Manufactured Solutions (MMS) Validation");
    println!("====================================================\n");

    println!("Following Roache (1998) and Salari & Knupp (2000),");
    println!("MMS enables code verification by constructing exact");
    println!("solutions with known source terms.\n");

    // Run different MMS tests
    mms_diffusion_order_verification()?;
    mms_advection_order_verification()?;
    mms_laplace_order_verification()?;

    println!("\n====================================================");
    println!("  All MMS Verification Tests Complete");
    println!("====================================================");

    Ok(())
}

/// Verify spatial discretization order for diffusion equation
///
/// ∂u/∂t = α∇²u + S(x,y,t)
///
/// where S is manufactured to match chosen exact solution
fn mms_diffusion_order_verification() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n1. DIFFUSION EQUATION - Order Verification");
    println!("   ----------------------------------------");

    let alpha = 0.1; // Thermal diffusivity
    let solution = ManufacturedDiffusion::new(alpha);

    // Manufactured solution: u(x,y,t) = sin(πx)sin(πy)exp(-t)
    println!("   Manufactured solution: u(x,y,t) = sin(πx)sin(πy)exp(-t)");
    println!("   Expected order: 2nd order (central differences)");

    // Test on multiple grid resolutions
    let resolutions = vec![16, 32, 64, 128];
    let mut grid_sizes = Vec::new();
    let mut errors = Vec::new();

    for &n in &resolutions {
        let _grid = StructuredGrid2D::<f64>::new(n, n, 0.0, 1.0, 0.0, 1.0)?;
        let dx = 1.0 / (n - 1) as f64;
        grid_sizes.push(dx);

        // CFL condition for stability
        let dt = 0.2 * dx * dx / alpha;
        let t_final = 0.01;
        let n_steps = (t_final / dt).ceil() as usize;

        // Initialize with manufactured solution at t=0
        let mut field = Field2D::new(n, n, 0.0);
        for i in 0..n {
            for j in 0..n {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                field.set(i, j, solution.initial_condition(x, y, 0.0));
            }
        }

        // Time integration with explicit Euler
        let mut t = 0.0;
        for _ in 0..n_steps {
            let mut field_new = field.clone();

            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    let laplacian = (field.at(i + 1, j) + field.at(i - 1, j)
                        + field.at(i, j + 1)
                        + field.at(i, j - 1)
                        - 4.0 * field.at(i, j))
                        / (dx * dx);

                    let x = i as f64 * dx;
                    let y = j as f64 * dx;
                    let source = solution.source_term(x, y, 0.0, t);

                    field_new.set(i, j, field.at(i, j) + dt * (alpha * laplacian + source));
                }
            }

            field = field_new;
            t += dt;
        }

        // Compute L2 error
        let mut l2_error = 0.0;
        let mut count = 0;

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                let exact = solution.exact_solution(x, y, 0.0, t);
                let error = field.at(i, j) - exact;
                l2_error += error * error;
                count += 1;
            }
        }

        l2_error = (l2_error / count as f64).sqrt();
        errors.push(l2_error);

        println!("   Grid {}x{}: dx={:.4e}, L2 error={:.4e}", n, n, dx, l2_error);
    }

    // Compute convergence rate
    let study = ConvergenceStudy::new(grid_sizes, errors)?;
    
    println!("\n   Convergence Analysis:");
    println!("     Observed order: {:.2}", study.convergence_rate);
    println!("     R²: {:.6}", study.r_squared);
    println!("     Asymptotic: {}", study.is_asymptotic());

    // Verify 2nd order convergence
    if (study.convergence_rate - 2.0).abs() < 0.3 {
        println!("     ✓ Spatial discretization verified as 2nd order");
    } else {
        println!("     ⚠ Observed order deviates from expected 2nd order");
        println!("       (may be limited by temporal discretization)");
    }

    Ok(())
}

/// Verify spatial discretization order for advection equation
///
/// ∂u/∂t + c·∇u = S(x,y,t)
fn mms_advection_order_verification() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n2. ADVECTION EQUATION - Order Verification");
    println!("   ---------------------------------------");

    let vx: f64 = 1.0; // Advection velocity x-component
    let vy: f64 = 0.5; // Advection velocity y-component
    let solution = ManufacturedAdvection::new(vx, vy);

    println!("   Manufactured solution: u(x,y,t) = sin(x-cₓt)cos(y-cyt)");
    println!("   Velocity: ({:.1}, {:.1})", vx, vy);
    println!("   Expected order: 1st order (upwind scheme)");

    let resolutions = vec![32, 64, 128];
    let mut grid_sizes = Vec::new();
    let mut errors = Vec::new();

    for &n in &resolutions {
        let _grid = StructuredGrid2D::<f64>::new(n, n, 0.0, 1.0, 0.0, 1.0)?;
        let dx = 1.0 / (n - 1) as f64;
        grid_sizes.push(dx);

        // CFL condition for advection
        let c_max = vx.abs().max(vy.abs());
        let dt = 0.5 * dx / c_max;
        let t_final = 0.1;
        let n_steps = (t_final / dt).ceil() as usize;

        // Initialize
        let mut field = Field2D::new(n, n, 0.0);
        for i in 0..n {
            for j in 0..n {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                field.set(i, j, solution.initial_condition(x, y, 0.0));
            }
        }

        // Time integration with upwind scheme
        // Note: Boundaries must be updated to exact solution at each timestep
        // to avoid accumulating errors from stale boundary values
        let mut t = 0.0;
        for _ in 0..n_steps {
            let mut field_new = field.clone();

            // Update interior points with upwind discretization
            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    // Upwind derivatives
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

                    let x = i as f64 * dx;
                    let y = j as f64 * dx;
                    let source = solution.source_term(x, y, 0.0, t);

                    let advection = vx * dudx + vy * dudy;
                    field_new.set(i, j, field.at(i, j) + dt * (source - advection));
                }
            }

            // Update boundary conditions to exact solution at new time
            t += dt;
            for i in 0..n {
                let x = i as f64 * dx;
                // South boundary (j=0)
                let y = 0.0;
                field_new.set(i, 0, solution.exact_solution(x, y, 0.0, t));
                // North boundary (j=n-1)
                let y = (n - 1) as f64 * dx;
                field_new.set(i, n - 1, solution.exact_solution(x, y, 0.0, t));
            }
            for j in 0..n {
                let y = j as f64 * dx;
                // West boundary (i=0)
                let x = 0.0;
                field_new.set(0, j, solution.exact_solution(x, y, 0.0, t));
                // East boundary (i=n-1)
                let x = (n - 1) as f64 * dx;
                field_new.set(n - 1, j, solution.exact_solution(x, y, 0.0, t));
            }

            field = field_new;
        }

        // Compute L2 error
        let mut l2_error = 0.0;
        let mut count = 0;

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                let exact = solution.exact_solution(x, y, 0.0, t);
                let error = field.at(i, j) - exact;
                l2_error += error * error;
                count += 1;
            }
        }

        l2_error = (l2_error / count as f64).sqrt();
        errors.push(l2_error);

        println!("   Grid {}x{}: dx={:.4e}, L2 error={:.4e}", n, n, dx, l2_error);
    }

    // Compute convergence rate
    let study = ConvergenceStudy::new(grid_sizes, errors)?;
    
    println!("\n   Convergence Analysis:");
    println!("     Observed order: {:.2}", study.convergence_rate);
    println!("     R²: {:.6}", study.r_squared);

    if (study.convergence_rate - 1.0).abs() < 0.3 {
        println!("     ✓ Upwind scheme verified as 1st order");
    } else {
        println!("     ⚠ Observed order deviates from expected 1st order");
    }

    Ok(())
}

/// Verify Laplace equation solver with manufactured solution
///
/// ∇²u = S(x,y)
fn mms_laplace_order_verification() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n3. LAPLACE EQUATION - Order Verification");
    println!("   -------------------------------------");

    println!("   Manufactured solution: u(x,y) = sin(2πx)sin(2πy)");
    println!("   Source term: S = -8π²sin(2πx)sin(2πy)");
    println!("   Expected order: 2nd order (5-point stencil)");

    let resolutions = vec![16, 32, 64, 128];
    let mut grid_sizes = Vec::new();
    let mut errors = Vec::new();

    for &n in &resolutions {
        let dx = 1.0 / (n - 1) as f64;
        grid_sizes.push(dx);

        // Build linear system Au = b for discrete Laplacian
        let n_interior = (n - 2) * (n - 2);
        let mut a_matrix = vec![vec![0.0; n_interior]; n_interior];
        let mut b_vector = vec![0.0; n_interior];

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let idx = (i - 1) * (n - 2) + (j - 1);
                
                // 5-point stencil
                a_matrix[idx][idx] = -4.0;
                
                if j > 1 {
                    let idx_left = (i - 1) * (n - 2) + (j - 2);
                    a_matrix[idx][idx_left] = 1.0;
                }
                if j < n - 2 {
                    let idx_right = (i - 1) * (n - 2) + j;
                    a_matrix[idx][idx_right] = 1.0;
                }
                if i > 1 {
                    let idx_down = (i - 2) * (n - 2) + (j - 1);
                    a_matrix[idx][idx_down] = 1.0;
                }
                if i < n - 2 {
                    let idx_up = i * (n - 2) + (j - 1);
                    a_matrix[idx][idx_up] = 1.0;
                }

                // Source term: -8π² sin(2πx)sin(2πy)
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                let source = -8.0 * PI * PI * (2.0 * PI * x).sin() * (2.0 * PI * y).sin();
                b_vector[idx] = source * dx * dx;
            }
        }

        // Simple iterative solver (Jacobi) for demonstration
        let mut u = vec![0.0; n_interior];
        for _ in 0..10000 {
            let mut u_new = u.clone();
            for i in 0..n_interior {
                let mut sum = b_vector[i];
                for j in 0..n_interior {
                    if i != j {
                        sum -= a_matrix[i][j] * u[j];
                    }
                }
                u_new[i] = sum / a_matrix[i][i];
            }
            
            // Check convergence
            let max_change: f64 = u.iter()
                .zip(u_new.iter())
                .map(|(old, new)| (new - old).abs())
                .fold(0.0, f64::max);
            
            u = u_new;
            
            if max_change < 1e-8 {
                break;
            }
        }

        // Compute L2 error against exact solution
        let mut l2_error = 0.0;
        let mut count = 0;

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let idx = (i - 1) * (n - 2) + (j - 1);
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                let exact = (2.0 * PI * x).sin() * (2.0 * PI * y).sin();
                let error = u[idx] - exact;
                l2_error += error * error;
                count += 1;
            }
        }

        l2_error = (l2_error / count as f64).sqrt();
        errors.push(l2_error);

        println!("   Grid {}x{}: dx={:.4e}, L2 error={:.4e}", n, n, dx, l2_error);
    }

    // Compute convergence rate
    let study = ConvergenceStudy::new(grid_sizes, errors)?;
    
    println!("\n   Convergence Analysis:");
    println!("     Observed order: {:.2}", study.convergence_rate);
    println!("     R²: {:.6}", study.r_squared);
    println!("     Asymptotic: {}", study.is_asymptotic());

    if (study.convergence_rate - 2.0).abs() < 0.2 {
        println!("     ✓ Discrete Laplacian verified as 2nd order");
    } else {
        println!("     ⚠ Observed order deviates from expected 2nd order");
    }

    Ok(())
}
