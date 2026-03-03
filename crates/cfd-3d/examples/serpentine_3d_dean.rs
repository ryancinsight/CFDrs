//! 3D Serpentine Channel — Dean Vortex Validation
//!
//! Solves 3D flow through a serpentine (wavy) channel and validates
//! Dean vortex formation via the Dean number and pressure drop.
//!
//! # Physics
//!
//! In curved channels, centrifugal forces create secondary flow (Dean vortices).
//! The Dean number De = Re * sqrt(d / (2R)) characterises this effect,
//! where R is the radius of curvature.
//!
//! # Reference
//!
//! Dean, W.R. (1928). "Fluid motion in a curved channel".
//! *Proc. R. Soc. Lond. A* 121(787):402-420.

use cfd_3d::serpentine::{SerpentineConfig3D, SerpentineSolver3D};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::application::channel::serpentine::SerpentineMeshBuilder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D Serpentine Channel — Dean Vortex Validation");
    println!("===============================================\n");

    // Serpentine geometry
    let diameter = 1.0e-3; // 1 mm channel diameter
    let amplitude = 2.0e-3; // 2 mm wave amplitude
    let wavelength = 8.0e-3; // 8 mm wavelength

    let builder = SerpentineMeshBuilder::new(diameter, amplitude, wavelength)
        .with_periods(3)
        .with_resolution(40, 6);

    println!("Geometry:");
    println!("  Diameter:    {:.3e} m", diameter);
    println!("  Amplitude:   {:.3e} m", amplitude);
    println!("  Wavelength:  {:.3e} m", wavelength);
    println!("  Periods:     3");

    // Solver configuration
    let config = SerpentineConfig3D {
        inlet_flow_rate: 5e-7, // 0.5 mL/s
        inlet_pressure: 200.0,
        outlet_pressure: 0.0,
        ..SerpentineConfig3D::default()
    };

    let inlet_flow_rate = config.inlet_flow_rate;
    let solver = SerpentineSolver3D::new(builder, config);

    // Solve with water
    let water = ConstantPropertyFluid::<f64>::water_20c()?;
    let rho = water.density;
    let mu = water.viscosity;
    println!("\nFluid: water at 20C (rho={:.1}, mu={:.4e})", rho, mu);

    println!("\nSolving...");
    let solution = solver.solve(water)?;

    println!("\nResults:");
    println!("  Inlet velocity:      {:.4} m/s", solution.u_inlet);
    println!("  Inlet pressure:      {:.1} Pa", solution.p_inlet);
    println!("  Outlet pressure:     {:.1} Pa", solution.p_outlet);
    println!("  Total pressure drop: {:.2} Pa", solution.dp_total);
    println!("  Dean number:         {:.2}", solution.dean_number);

    // Analytical estimates
    let area = std::f64::consts::PI * (diameter / 2.0).powi(2);
    let u_mean = inlet_flow_rate / area;
    let re = rho * u_mean * diameter / mu;
    let r_curve = (wavelength / (2.0 * std::f64::consts::PI)).max(amplitude);
    let de_analytical = re * (diameter / (2.0 * r_curve)).sqrt();

    println!("\nAnalytical Estimates:");
    println!("  Mean velocity (from Q/A): {:.4} m/s", u_mean);
    println!("  Reynolds number:          {:.1}", re);
    println!("  Dean number (analytical): {:.2}", de_analytical);
    println!("  Dean number (solver):     {:.2}", solution.dean_number);

    // Dean vortex assessment
    if solution.dean_number > 40.0 {
        println!("\n  Dean vortices expected (De > 40)");
    } else {
        println!("\n  Weak secondary flow (De < 40)");
    }

    // Pressure drop per unit length
    let total_length = wavelength * 3.0; // approximate
    let dp_per_m = solution.dp_total / total_length;
    println!("  Pressure gradient: {:.1} Pa/m", dp_per_m);

    println!("\nSerpentine 3D example completed.");
    Ok(())
}
