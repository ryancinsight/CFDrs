//! External Package Validation: Python_CFD Comparison
//!
//! This example replicates key problems from the external/Python_CFD/ notebooks
//! (Lorena Barba's "12 Steps to Navier-Stokes") and solves them with CFDrs
//! solvers, comparing results against analytical solutions.
//!
//! Problems validated:
//! 1. 2D Poiseuille Channel Flow (Python_CFD Notebook 16)
//! 2. 2D Poisson Equation (Python_CFD Notebook 14)
//! 3. Ghia Lid-Driven Cavity Reference Data (Python_CFD Ghia-1982.txt)

use std::collections::HashMap;

// --- Problem 1: Poiseuille flow ---
use cfd_2d::solvers::poiseuille::{BloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
use cfd_core::physics::fluid::blood::CassonBlood;

// --- Problem 2: Poisson equation ---
use cfd_2d::grid::{Grid2D, StructuredGrid2D};
use cfd_2d::solvers::fdm::{FdmConfig, PoissonSolver};

// --- Problem 3: Ghia reference data ---
use cfd_validation::benchmarks::LidDrivenCavity;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("External Package Validation: Python_CFD Comparison");
    println!("====================================================");
    println!();

    let mut pass_count = 0u32;
    let mut total_count = 0u32;

    // ====================================================================
    // Problem 1: 2D Poiseuille Channel Flow
    // ====================================================================
    total_count += 1;
    match run_poiseuille_validation() {
        Ok(passed) => {
            if passed {
                pass_count += 1;
            }
        }
        Err(e) => {
            println!("  Status: FAIL (solver error: {})", e);
        }
    }
    println!();

    // ====================================================================
    // Problem 2: 2D Poisson Equation
    // ====================================================================
    total_count += 1;
    match run_poisson_validation() {
        Ok(passed) => {
            if passed {
                pass_count += 1;
            }
        }
        Err(e) => {
            println!("  Status: FAIL (solver error: {})", e);
        }
    }
    println!();

    // ====================================================================
    // Problem 3: Ghia Lid-Driven Cavity Reference Data
    // ====================================================================
    total_count += 1;
    if run_ghia_reference_data() {
        pass_count += 1;
    }
    println!();

    // ====================================================================
    // Summary
    // ====================================================================
    println!("====================================================");
    println!("Summary: {}/{} problems PASS", pass_count, total_count);
    println!();

    if pass_count == total_count {
        println!("All validations passed.");
    } else {
        println!("WARNING: Some validations failed. Review results above.");
    }

    Ok(())
}

// ========================================================================
// Problem 1: 2D Poiseuille Channel Flow
// ========================================================================

/// Validate the Poiseuille flow solver against the analytical solution.
///
/// Physical setup (matching Python_CFD Notebook 16):
///   - Parallel plates, channel height H = 100 um
///   - Pressure gradient |dP/dx| = 1000 Pa/m
///   - Fluid: water-like (mu = 0.001 Pa.s, rho = 1000 kg/m^3)
///
/// Analytical solution (Newtonian):
///   u(y) = (1 / 2mu) * (dP/dx) * y * (H - y)
///   u_max = H^2 * |dP/dx| / (8 * mu)
///   Q/W   = H^3 * |dP/dx| / (12 * mu)
///
/// To emulate a Newtonian fluid with the non-Newtonian Poiseuille solver,
/// we use a Casson blood model with zero yield stress. When tau_y = 0 the
/// Casson apparent viscosity reduces to:
///   mu_app = (sqrt(0) / sqrt(gamma_dot) + sqrt(mu_inf))^2 = mu_inf
/// which is constant -- i.e., Newtonian.
fn run_poiseuille_validation() -> Result<bool, Box<dyn std::error::Error>> {
    println!("Problem 1: 2D Poiseuille Channel Flow");
    println!("  Reference: Python_CFD Notebook 16, Lorena Barba (2014)");

    let height: f64 = 100e-6; // 100 um
    let mu: f64 = 0.001; // Pa.s (water viscosity)
    let dp_dx: f64 = 1000.0; // Pa/m
    let rho: f64 = 1000.0; // kg/m^3

    println!(
        "  Setup: H={:.0}um, dP/dx={} Pa/m, mu={} Pa.s, rho={} kg/m^3",
        height * 1e6,
        dp_dx,
        mu,
        rho
    );

    // Create a Newtonian-equivalent Casson model (zero yield stress)
    let water_as_casson = CassonBlood::new(
        rho, // density
        0.0, // yield_stress = 0 => Newtonian
        mu,  // infinite_shear_viscosity = water viscosity
        0.0, // hematocrit = 0 (not blood)
    );

    let config = PoiseuilleConfig {
        height,
        width: 500e-6,
        length: 1e-3,
        ny: 101,
        pressure_gradient: dp_dx,
        tolerance: 1e-8,
        max_iterations: 2000,
        relaxation_factor: 0.7,
    };

    let mut solver = PoiseuilleFlow2D::new(config, BloodModel::Casson(water_as_casson));
    let iterations = solver.solve()?;

    // --- Compare peak velocity ---
    let u_max_numerical = solver.max_velocity();
    let u_max_analytical = height.powi(2) * dp_dx / (8.0 * mu);
    let error_umax = ((u_max_numerical - u_max_analytical) / u_max_analytical).abs();

    println!("  Analytical u_max = {:.6e} m/s", u_max_analytical);
    println!("  Numerical u_max  = {:.6e} m/s", u_max_numerical);
    println!("  Relative error   = {:.2e}", error_umax);
    println!("  Iterations       = {}", iterations);

    // --- Compare full velocity profile ---
    let analytical_profile = solver.analytical_solution(mu);
    let velocity_profile = solver.velocity_profile();

    let mut max_profile_error: f64 = 0.0;
    for i in 1..velocity_profile.len() - 1 {
        let ana = analytical_profile[i];
        if ana.abs() > 1e-20 {
            let err = ((velocity_profile[i] - ana) / ana).abs();
            if err > max_profile_error {
                max_profile_error = err;
            }
        }
    }
    println!("  Profile max relative error = {:.2e}", max_profile_error);

    // --- Compare volumetric flow rate ---
    let q_per_w_numerical = solver.flow_rate_per_width();
    let q_per_w_analytical = height.powi(3) * dp_dx / (12.0 * mu);
    let error_q = ((q_per_w_numerical - q_per_w_analytical) / q_per_w_analytical).abs();

    println!(
        "  Flow rate Q/W: numerical={:.6e}, analytical={:.6e}, error={:.2e}",
        q_per_w_numerical, q_per_w_analytical, error_q
    );

    // --- Compare wall shear stress ---
    let tau_w_numerical = solver.wall_shear_stress();
    // Analytical: tau_w = mu * du/dy |_{y=0} = mu * (H * dP/dx) / (2 * mu) = H * dP/dx / 2
    let tau_w_analytical = height * dp_dx / 2.0;
    let error_tau = ((tau_w_numerical - tau_w_analytical) / tau_w_analytical).abs();

    println!(
        "  Wall shear stress: numerical={:.6e}, analytical={:.6e}, error={:.2e}",
        tau_w_numerical, tau_w_analytical, error_tau
    );

    // --- Reynolds number for reference ---
    let re = rho * u_max_analytical * height / mu;
    println!("  Reynolds number (analytical) = {:.4}", re);

    // --- Pass/fail ---
    let passed = error_umax < 0.01; // 1% tolerance
    if passed {
        println!("  Status: PASS (u_max error < 1%)");
    } else {
        println!(
            "  Status: FAIL (u_max error = {:.2}%, expected < 1%)",
            error_umax * 100.0
        );
    }

    Ok(passed)
}

// ========================================================================
// Problem 2: 2D Poisson Equation
// ========================================================================

/// Validate the FDM Poisson solver against a known manufactured solution.
///
/// Physical setup (analogous to Python_CFD Notebook 14):
///   - Domain: unit square [0, 1]^2
///   - Manufactured solution: phi(x,y) = sin(pi*x) * sin(pi*y)
///   - Source term: f = -2 * pi^2 * sin(pi*x) * sin(pi*y)
///     (since nabla^2 phi = -pi^2 sin(pi*x) sin(pi*y) - pi^2 sin(pi*x) sin(pi*y))
///   - Dirichlet BCs: phi = sin(pi*x) * sin(pi*y) on all boundaries
///
/// The Poisson equation is:  nabla^2 phi = f
/// with 5-point finite difference stencil on a structured grid.
///
/// NOTE: The StructuredGrid2D uses cell-centered coordinates, so boundary
/// cell centers are offset slightly from the domain edges. The Dirichlet
/// boundary values are evaluated at these cell centers, introducing a small
/// (order dx) offset from the true boundary. This is standard for
/// cell-centered FDM/FVM discretizations.
fn run_poisson_validation() -> Result<bool, Box<dyn std::error::Error>> {
    println!("Problem 2: 2D Poisson Equation");
    println!("  Reference: Python_CFD Notebook 14, Lorena Barba (2014)");
    println!(
        "  Setup: Unit square, phi=sin(pi*x)*sin(pi*y), source=-2*pi^2*sin(pi*x)*sin(pi*y)"
    );

    let nx: usize = 41;
    let ny: usize = 41;
    let grid = StructuredGrid2D::<f64>::unit_square(nx, ny)?;

    let pi = std::f64::consts::PI;

    // Build source term and Dirichlet boundary values
    let mut source: HashMap<(usize, usize), f64> = HashMap::new();
    let mut boundary_values: HashMap<(usize, usize), f64> = HashMap::new();

    for (i, j) in grid.iter() {
        let center = grid.cell_center(i, j)?;
        let x = center.x;
        let y = center.y;

        if grid.is_boundary(i, j) {
            // Dirichlet BC evaluated at cell center
            boundary_values.insert((i, j), (pi * x).sin() * (pi * y).sin());
        } else {
            // Source term: nabla^2 phi = -2*pi^2 * sin(pi*x) * sin(pi*y)
            source.insert((i, j), -2.0 * pi * pi * (pi * x).sin() * (pi * y).sin());
        }
    }

    // Configure solver: increase max iterations and use SOR for faster convergence
    let mut fdm_config = FdmConfig::<f64>::default();
    fdm_config.base.convergence.max_iterations = 10_000;
    fdm_config.base.convergence.tolerance = 1e-8;
    fdm_config.base.numerical.relaxation = 1.5; // SOR over-relaxation

    let poisson_solver = PoissonSolver::new(fdm_config);
    let result = poisson_solver.solve(&grid, &source, &boundary_values)?;

    // Compute error metrics over interior points
    let mut max_err: f64 = 0.0;
    let mut l2_err_sum: f64 = 0.0;
    let mut n_interior: usize = 0;

    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let center = grid.cell_center(i, j)?;
            let x = center.x;
            let y = center.y;
            let phi_exact = (pi * x).sin() * (pi * y).sin();
            let phi_num = result.get(&(i, j)).copied().unwrap_or(0.0);
            let err = (phi_num - phi_exact).abs();
            if err > max_err {
                max_err = err;
            }
            l2_err_sum += err * err;
            n_interior += 1;
        }
    }

    let l2_err = (l2_err_sum / n_interior as f64).sqrt();

    println!("  Grid: {}x{} ({} interior points)", nx, ny, n_interior);
    println!("  L-infinity error = {:.4e}", max_err);
    println!("  L2 error         = {:.4e}", l2_err);

    // Accept error < 5% of the peak value (peak of sin*sin is 1.0)
    // Cell-centered grids introduce O(dx) boundary offset, plus iterative solver tolerance
    let passed = max_err < 0.05;
    if passed {
        println!("  Status: PASS (L-inf error < 5%)");
    } else {
        println!(
            "  Status: FAIL (L-inf error = {:.4e}, expected < 0.05)",
            max_err
        );
    }

    Ok(passed)
}

// ========================================================================
// Problem 3: Ghia Lid-Driven Cavity Reference Data
// ========================================================================

/// Display and validate Ghia et al. (1982) reference data for the
/// lid-driven cavity benchmark at Re = 100.
///
/// This is the standard CFD validation case. The Python_CFD notebooks
/// compare against the same tabulated data from:
///   Ghia, Ghia, and Shin (1982) "High-Re solutions for incompressible
///   flow using the Navier-Stokes equations and a multigrid method",
///   J. Comp. Phys. 48, pp. 387-411.
///
/// We verify that the reference data is correctly loaded from the
/// cfd-validation crate and display it for cross-checking against the
/// Python_CFD/Ghia-1982.txt file.
fn run_ghia_reference_data() -> bool {
    println!("Problem 3: Ghia Lid-Driven Cavity Reference Data (Re=100)");
    println!("  Reference: Ghia et al. (1982), Python_CFD/Ghia-1982.txt");
    println!("  Setup: Unit square, lid velocity U=1, Re=100");
    println!();

    let cavity = LidDrivenCavity::<f64>::new(1.0, 1.0, 100.0);

    let ghia_u = cavity.ghia_u_centerline(100.0);
    let ghia_v = cavity.ghia_v_centerline(100.0);

    // Display u-velocity along vertical centerline
    println!("  Ghia u-velocity along vertical centerline (x=0.5):");
    println!("    {:>8}  {:>10}", "y", "u/U");
    println!("    {:>8}  {:>10}", "--------", "----------");
    for (y_val, u_val) in &ghia_u {
        println!("    {:>8.4}  {:>10.5}", y_val, u_val);
    }

    println!();

    // Display v-velocity along horizontal centerline
    println!("  Ghia v-velocity along horizontal centerline (y=0.5):");
    println!("    {:>8}  {:>10}", "x", "v/U");
    println!("    {:>8}  {:>10}", "--------", "----------");
    for (x_val, v_val) in &ghia_v {
        println!("    {:>8.4}  {:>10.5}", x_val, v_val);
    }

    println!();

    // Validate that reference data was loaded correctly
    let u_ok = !ghia_u.is_empty();
    let v_ok = !ghia_v.is_empty();

    // Spot-check known values from the original Ghia (1982) paper
    let mut spot_checks_pass = true;

    if let Some(&(y, u)) = ghia_u.first() {
        // First entry: y=1.0, u=1.0 (lid velocity)
        if (y - 1.0).abs() > 1e-10 || (u - 1.0).abs() > 1e-10 {
            println!(
                "  WARNING: Ghia u-data first point mismatch: y={}, u={} (expected 1.0, 1.0)",
                y, u
            );
            spot_checks_pass = false;
        }
    }

    if let Some(&(y, u)) = ghia_u.last() {
        // Last entry: y=0.0, u=0.0 (no-slip wall)
        if (y).abs() > 1e-10 || (u).abs() > 1e-10 {
            println!(
                "  WARNING: Ghia u-data last point mismatch: y={}, u={} (expected 0.0, 0.0)",
                y, u
            );
            spot_checks_pass = false;
        }
    }

    // Check that the minimum u-velocity is around -0.21 (known for Re=100)
    if let Some(u_min) = ghia_u.iter().map(|(_, u)| *u).reduce(f64::min) {
        if (u_min - (-0.21090)).abs() > 0.01 {
            println!(
                "  WARNING: Ghia u_min = {} (expected ~-0.21090 for Re=100)",
                u_min
            );
            spot_checks_pass = false;
        }
    }

    let passed = u_ok && v_ok && spot_checks_pass;
    if passed {
        println!(
            "  Status: PASS (reference data loaded: {} u-points, {} v-points, spot-checks OK)",
            ghia_u.len(),
            ghia_v.len()
        );
    } else {
        println!("  Status: FAIL (reference data incomplete or spot-checks failed)");
    }

    passed
}
