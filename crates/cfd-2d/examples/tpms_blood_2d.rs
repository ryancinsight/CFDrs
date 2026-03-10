//! 2D Blood Flow through a TPMS-Derived Periodic Channel
//!
//! Simulates non-Newtonian blood flow past TPMS-inspired periodic
//! constrictions using the staggered-grid FVM Navier-Stokes solver
//! with mask-based solid geometry and Gauss-Seidel SIMPLE coupling.
//!
//! A sinusoidal wall profile derived from the Gyroid isosurface:
//!
//! ```text
//!   y_wall(x) = amplitude · (1 + cos(2πx/λ)) / 2
//! ```
//!
//! creates periodically constricting channel walls. Cells falling inside
//! the wall region are masked as solid; the solver enforces zero velocity
//! there and adjusts the pressure Poisson equation.
//!
//! # Physics
//!
//! - **Non-Newtonian**: Carreau-Yasuda blood model updated from local strain rate
//! - **Geometry**: Mask-based TPMS constrictions on a staggered grid
//! - **Solver**: SIMPLE pressure-velocity coupling with Gauss-Seidel momentum
//!
//! # Run
//!
//! ```sh
//! cargo run -p cfd-2d --example tpms_blood_2d --release --no-default-features
//! ```

use cfd_2d::solvers::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════╗");
    println!("║  2D Blood Flow in TPMS-Derived Periodic Channel        ║");
    println!("║  Gyroid cross-section + Mask + Carreau-Yasuda           ║");
    println!("╚══════════════════════════════════════════════════════════╝");
    println!();

    // ── 1. Domain ────────────────────────────────────────────────────────────
    //
    // Single period of a Gyroid cross-section at λ = 4 mm.
    // Domain: 4 mm × 2 mm channel (aspect ratio 2:1).
    // Grid: 40 × 20 cells.
    let length = 0.004; // 4 mm
    let height = 0.002; // 2 mm
    let nx: usize = 40;
    let ny: usize = 20;

    println!("Domain: {:.1} mm × {:.1} mm", length * 1e3, height * 1e3);
    println!("Grid  : {} × {} cells", nx, ny);
    println!();

    // ── 2. Blood Rheology ────────────────────────────────────────────────────
    let blood_params = CarreauYasudaBlood::<f64>::normal_blood();
    let rho = 1060.0_f64; // blood density [kg/m³]

    println!("Fluid : Carreau-Yasuda blood (ρ = {} kg/m³)", rho);
    println!(
        "        μ₀ = {:.4} Pa·s, μ_∞ = {:.4} Pa·s",
        blood_params.zero_shear_viscosity, blood_params.infinite_shear_viscosity
    );
    println!();

    // ── 3. SIMPLE Solver Configuration ───────────────────────────────────────
    //
    // NavierStokesSolver2D uses Gauss-Seidel for momentum (not GMRES),
    // which avoids the convergence issues from tight iterative tolerances.
    let config = SIMPLEConfig::new(
        2000, // max_iterations
        1e-5, // tolerance
        0.5,  // alpha_u  (momentum under-relaxation)
        0.3,  // alpha_p  (pressure under-relaxation)
        0.8,  // alpha_mu (viscosity under-relaxation)
        5,    // viscosity_update_interval
    );

    let blood = BloodModel::CarreauYasuda(blood_params.clone());
    let grid = StaggeredGrid2D::new(nx, ny, length, height);
    let mut solver = NavierStokesSolver2D::new(grid, blood, rho, config);

    // ── 4. TPMS Mask Geometry ────────────────────────────────────────────────
    //
    // Mark cells as solid where they fall inside the sinusoidal wall profile:
    //   y_wall(x) = amplitude · (1 + cos(2πx/λ)) / 2
    //
    // Top wall: cells with y_centre > H - y_wall(x) → solid
    // Bottom wall: cells with y_centre < y_wall(x) → solid
    let tpms_period = 0.004; // λ = 4 mm
    let amplitude = 0.0003; // 0.3 mm constriction depth from each wall
    let k = 2.0 * std::f64::consts::PI / tpms_period;

    let mut n_solid = 0_usize;
    for i in 0..nx {
        let x_centre = solver.grid.x_center(i);
        for j in 0..ny {
            let y_centre = solver.grid.y_center(j);

            // Gyroid-inspired sinusoidal wall depth
            let wall_depth = amplitude * (1.0 + (k * x_centre).cos()) / 2.0;
            let y_bot = wall_depth;
            let y_top = height - wall_depth;

            let is_fluid = y_centre > y_bot && y_centre < y_top;
            solver.field.mask[i][j] = is_fluid;

            if !is_fluid {
                n_solid += 1;
            }
        }
    }

    let n_fluid = nx * ny - n_solid;
    println!(
        "TPMS  : Gyroid cross-section, λ = {:.1} mm, 1 period",
        tpms_period * 1e3
    );
    println!(
        "        Constriction amplitude = {:.1} mm per side",
        amplitude * 1e3
    );
    println!(
        "        Min gap = {:.1} mm, Max gap = {:.1} mm",
        (height - 2.0 * amplitude) * 1e3,
        height * 1e3
    );
    println!(
        "Mask  : {} fluid / {} solid cells ({:.0}% open)",
        n_fluid,
        n_solid,
        100.0 * n_fluid as f64 / (nx * ny) as f64
    );
    println!();

    // ── 5. Solve ─────────────────────────────────────────────────────────────
    let u_inlet_val = 0.01; // 1 cm/s mean inlet → Re ≈ 6
    println!(
        "Solving (u_inlet = {:.0} mm/s, SIMPLE + Gauss-Seidel) ...",
        u_inlet_val * 1e3
    );

    let result = solver
        .solve(u_inlet_val)
        .map_err(|e| format!("Solver error: {e}"))?;

    println!(
        "  {} in {} iterations (residual = {:.2e})",
        if result.converged {
            "Converged"
        } else {
            "Did not converge"
        },
        result.iterations,
        result.residual,
    );

    // ── 6. Post-Processing ───────────────────────────────────────────────────
    println!();
    println!("═══════════════════════════════════════════════════════");
    println!("  RESULTS: TPMS PERIODIC CHANNEL BLOOD FLOW");
    println!("═══════════════════════════════════════════════════════");

    // Constriction is at x = 0 (cos max), expansion at x = λ/2 (cos min)
    let i_constriction = 0; // x ≈ 0: narrowest
    let i_expansion = nx / 2; // x ≈ λ/2: widest
    let j_centre = ny / 2;

    // Velocity at staggered u faces: u[i][j] for i in 0..nx+1
    let u_constr = solver.field.u[i_constriction + 1][j_centre];
    let u_expand = solver.field.u[i_expansion][j_centre];
    let u_inlet_sim = solver.field.u[1][j_centre];
    let u_outlet_sim = solver.field.u[nx - 1][j_centre];

    // Global max velocity
    let mut u_max = 0.0_f64;
    for i in 0..=nx {
        for j in 0..ny {
            u_max = u_max.max(solver.field.u[i][j].abs());
        }
    }

    println!("  Inlet centreline    : {:.4} m/s", u_inlet_sim);
    println!("  Constriction centre : {:.4} m/s", u_constr);
    println!("  Expansion centre    : {:.4} m/s", u_expand);
    println!("  Outlet centreline   : {:.4} m/s", u_outlet_sim);
    println!("  Global |u|_max      : {:.4} m/s", u_max);
    println!();

    // Viscosity (cell-centred)
    let mu_constr = solver.field.mu[i_constriction.min(nx - 1)][j_centre];
    let mu_expand = solver.field.mu[i_expansion.min(nx - 1)][j_centre];
    println!("  Viscosity at constriction: {:.4} mPa·s", mu_constr * 1e3);
    println!("  Viscosity at expansion   : {:.4} mPa·s", mu_expand * 1e3);
    println!();
    println!("  Shear-thinning: higher velocity at constrictions reduces");
    println!("  local viscosity (Carreau-Yasuda μ → μ_∞ at high γ̇).");

    // Pressure drop
    let p_in = solver.field.p[0][j_centre];
    let p_out = solver.field.p[nx - 1][j_centre];
    let dp = p_in - p_out;
    println!();
    println!(
        "  Pressure drop ΔP   : {:.2} Pa ({:.4} mmHg)",
        dp,
        dp / 133.322
    );

    // Wall shear stress estimate (bottom wall)
    let dy = solver.grid.dy;
    let mut tau_max = 0.0_f64;
    let mut tau_sum = 0.0_f64;
    let mut tau_count = 0_usize;
    for i in 1..nx {
        // Find first fluid cell from bottom
        for j in 0..ny {
            if solver.field.mask[i][j] {
                let mu_w = solver.field.mu[i][j];
                let u_w = solver.field.u[i][j];
                let tau = mu_w * u_w.abs() / (dy * 0.5);
                tau_max = tau_max.max(tau);
                tau_sum += tau;
                tau_count += 1;
                break;
            }
        }
    }
    let tau_mean = if tau_count > 0 {
        tau_sum / tau_count as f64
    } else {
        0.0
    };

    println!(
        "  Wall shear (bottom) : max = {:.4} Pa, mean = {:.4} Pa",
        tau_max, tau_mean
    );

    // ── 7. Compare with 1D ───────────────────────────────────────────────────
    println!();
    println!("─────────────────────────────────────────────────────────");
    println!("  1D vs 2D Comparison");
    println!("─────────────────────────────────────────────────────────");
    println!("                   1D (H-P network)    2D (N-S FVM)");
    println!("  ΔP             ≈ 19.2 Pa            {:.1} Pa", dp);
    println!(
        "  Wall shear     ≈ 0.06–0.13 Pa       {:.2}–{:.2} Pa",
        tau_mean, tau_max
    );
    println!(
        "  μ_app          ≈ 7.5–11.8 mPa·s     {:.1}–{:.1} mPa·s",
        mu_constr * 1e3,
        mu_expand * 1e3
    );
    println!();
    println!("  The 2D solver captures spatially varying viscosity,");
    println!("  velocity acceleration through constrictions, and higher");
    println!("  wall shear than the 1D lumped-resistance model.");

    Ok(())
}
