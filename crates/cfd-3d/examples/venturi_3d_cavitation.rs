//! 3D Venturi Flow with Cavitation Analysis
//!
//! Solves 3D flow through a venturi throat and computes cavitation number.
//! Validates against Bernoulli analytical solution for inviscid pressure drop.
//!
//! # Physics
//!
//! Cavitation number sigma = (p - p_v) / (0.5 * rho * u^2)
//! where p_v is the vapour pressure (~2340 Pa for water at 20C).
//!
//! # Reference
//!
//! ISO 5167-4:2022 — Measurement of fluid flow by means of pressure
//! differential devices (Venturi tubes).

use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::application::channel::venturi::VenturiMeshBuilder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D Venturi Flow — Cavitation Analysis");
    println!("======================================\n");

    // Venturi geometry (ISO 5167 like)
    let d_inlet: f64 = 2.0e-3;     // 2 mm inlet diameter
    let d_throat: f64 = 1.0e-3;    // 1 mm throat diameter
    let l_inlet = 3.0e-3;     // 3 mm inlet section
    let l_convergent = 2.0e-3;
    let l_throat = 1.0e-3;
    let l_divergent = 4.0e-3;
    let l_outlet = 3.0e-3;

    let builder = VenturiMeshBuilder::new(
        d_inlet, d_throat, l_inlet, l_convergent, l_throat, l_divergent, l_outlet,
    );

    println!("Geometry:");
    println!("  Inlet diameter:  {:.3e} m", d_inlet);
    println!("  Throat diameter: {:.3e} m", d_throat);
    let area_ratio_display: f64 = (d_throat / d_inlet).powi(2);
    println!("  Area ratio:      {:.3}", area_ratio_display);

    // Solver configuration
    let config = VenturiConfig3D {
        inlet_flow_rate: 1e-7, // 0.1 mL/s
        inlet_pressure: 101325.0,
        outlet_pressure: 101225.0,
        ..VenturiConfig3D::default()
    };

    let solver = VenturiSolver3D::new(builder, config);

    // Solve with water
    let water = ConstantPropertyFluid::<f64>::water_20c()?;
    let rho = water.density;
    println!("\nFluid: water at 20C (rho={:.1} kg/m^3, mu={:.4e} Pa.s)", rho, water.viscosity);

    println!("\nSolving...");
    let solution = solver.solve(water)?;

    println!("\nResults:");
    println!("  Inlet velocity:      {:.4} m/s", solution.u_inlet);
    println!("  Throat velocity:     {:.4} m/s", solution.u_throat);
    println!("  Inlet pressure:      {:.1} Pa", solution.p_inlet);
    println!("  Throat pressure:     {:.1} Pa", solution.p_throat);
    println!("  Outlet pressure:     {:.1} Pa", solution.p_outlet);
    println!("  Pressure drop (throat): {:.2} Pa", solution.dp_throat);
    println!("  Pressure recovery:      {:.2} Pa", solution.dp_recovery);
    println!("  Cp (throat):         {:.4}", solution.cp_throat);
    println!("  Mass error:          {:.2e}", solution.mass_error);

    // Bernoulli analytical comparison
    let area_ratio = (d_throat / d_inlet).powi(2);
    let u_throat_bernoulli = solution.u_inlet / area_ratio;
    let dp_bernoulli = 0.5 * rho * (u_throat_bernoulli.powi(2) - solution.u_inlet.powi(2));

    println!("\nBernoulli Comparison:");
    println!("  Analytical throat velocity: {:.4} m/s", u_throat_bernoulli);
    println!("  Analytical pressure drop:   {:.2} Pa", dp_bernoulli);
    println!("  Numerical pressure drop:    {:.2} Pa", solution.dp_throat);

    // Cavitation analysis
    let p_vapour = 2340.0; // Pa, water at 20C
    let sigma_inlet = (solution.p_inlet - p_vapour) / (0.5 * rho * solution.u_inlet.powi(2));
    let sigma_throat = (solution.p_throat - p_vapour) / (0.5 * rho * solution.u_throat.powi(2));

    println!("\nCavitation Number Table:");
    println!("  {:>12} {:>12} {:>12} {:>12}", "Location", "P [Pa]", "u [m/s]", "sigma");
    println!("  {:>12} {:>12.1} {:>12.4} {:>12.2}", "Inlet", solution.p_inlet, solution.u_inlet, sigma_inlet);
    println!("  {:>12} {:>12.1} {:>12.4} {:>12.2}", "Throat", solution.p_throat, solution.u_throat, sigma_throat);

    if sigma_throat < 1.0 {
        println!("\n  CAVITATION PREDICTED at throat (sigma < 1)");
    } else {
        println!("\n  No cavitation at throat (sigma > 1)");
    }

    println!("\nVenturi 3D example completed.");
    Ok(())
}
