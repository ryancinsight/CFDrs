//! Serpentine-channel mixing analysis.
//!
//! The analytical calculation uses the canonical transverse-diffusion model.
//! The second case drives the public discretized flow and scalar-transport
//! solver.  Each result is calculated from the selected geometry and material
//! parameters; this example does not use preset convergence or mixing values.

use std::error::Error;
use std::io;

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::serpentine_flow::{
    AdvectionDiffusionMixing, SerpentineGeometry, SerpentineSolver2D,
};

/// Mean inlet velocity, in m/s.
const INLET_VELOCITY: f64 = 0.01;
/// Density of the aqueous carrier fluid, in kg/m³.
const CARRIER_DENSITY: f64 = 1000.0;
/// Dynamic viscosity of water, in Pa·s.
const WATER_VISCOSITY: f64 = 1.0e-3;
/// Molecular diffusivity of a small aqueous solute, in m²/s.
const SMALL_SOLUTE_DIFFUSIVITY: f64 = 1.0e-9;
/// Diffusivity selected to resolve a measurable scalar field on the coarse grid.
const DISCRETIZED_DIFFUSIVITY: f64 = 1.0e-6;
/// Bisection error bound after the provider's 80 iterations.
const MIXING_FRACTION_TOLERANCE: f64 = 5.0e-12;

fn require(condition: bool, message: impl Into<String>) -> Result<(), Box<dyn Error>> {
    if condition {
        Ok(())
    } else {
        Err(io::Error::other(message.into()).into())
    }
}

fn report_analytical_mixing() -> Result<(), Box<dyn Error>> {
    let geometry = SerpentineGeometry::<f64>::microfluidic_standard();
    let model =
        AdvectionDiffusionMixing::new(geometry.width, INLET_VELOCITY, SMALL_SOLUTE_DIFFUSIVITY);
    let total_length = geometry.total_length();
    let mixing_length_90 = model.mixing_length_90_percent();
    let mixing_fraction_at_length = model.mixing_fraction(total_length);
    let mixing_fraction_at_target = model.mixing_fraction(mixing_length_90);

    require(
        (mixing_fraction_at_target - 0.9).abs() <= MIXING_FRACTION_TOLERANCE,
        format!("L90 must produce 90% mixing, got {mixing_fraction_at_target:.15}"),
    )?;

    println!("\nAnalytical transverse-diffusion model");
    println!("  channel width:       {:.0} μm", geometry.width * 1.0e6);
    println!("  total length:        {:.3} mm", total_length * 1.0e3);
    println!("  Peclet number:       {:.3e}", model.peclet_number());
    println!("  L90:                 {:.3} mm", mixing_length_90 * 1.0e3);
    println!(
        "  outlet mixing:       {:.3}%",
        mixing_fraction_at_length * 100.0
    );
    println!(
        "  mixing at L90:       {:.12}%",
        mixing_fraction_at_target * 100.0
    );

    Ok(())
}

fn report_discretized_mixing() -> Result<(), Box<dyn Error>> {
    let geometry = SerpentineGeometry::new(200e-6, 100e-6, 500e-6, 200e-6, 1);
    let mut solver = SerpentineSolver2D::new(
        geometry.clone(),
        BloodModel::Newtonian(WATER_VISCOSITY),
        CARRIER_DENSITY,
        40,
        20,
    );
    let solution = solver.solve(INLET_VELOCITY, DISCRETIZED_DIFFUSIVITY, 0.0, 1.0)?;

    require(
        solution.peclet.is_finite() && solution.peclet > 0.0,
        format!("solver produced invalid Peclet number {}", solution.peclet),
    )?;
    require(
        solution.l_mix_90.is_finite() && solution.l_mix_90 > 0.0,
        format!("solver produced invalid L90 {}", solution.l_mix_90),
    )?;
    require(
        solution.mixing_fraction_outlet.is_finite()
            && (0.0..=1.0).contains(&solution.mixing_fraction_outlet),
        format!(
            "solver produced invalid outlet mixing fraction {}",
            solution.mixing_fraction_outlet
        ),
    )?;
    require(
        solution.pressure_drop.is_finite(),
        format!(
            "solver produced invalid pressure drop {}",
            solution.pressure_drop
        ),
    )?;

    println!("\nDiscretized flow and scalar transport");
    println!(
        "  channel length:      {:.3} mm",
        geometry.total_length() * 1.0e3
    );
    println!("  grid:                40 × 20");
    println!("  Peclet number:       {:.3e}", solution.peclet);
    println!("  analytical L90:      {:.3} mm", solution.l_mix_90 * 1.0e3);
    println!(
        "  numerical mixing:    {:.3}%",
        solution.mixing_fraction_outlet * 100.0
    );
    println!("  pressure drop:       {:.6e} Pa", solution.pressure_drop);

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    println!("Serpentine mixing validation");
    report_analytical_mixing()?;
    report_discretized_mixing()?;
    println!("\nAll serpentine calculations completed.");
    Ok(())
}
