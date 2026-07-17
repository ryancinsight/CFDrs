//! Venturi analytical-model validation with blood rheology.
//!
//! The Venturi calculations use the public cfd-2d Bernoulli and viscous-loss
//! models.  The rheology calculation evaluates the public Casson and
//! Carreau-Yasuda constitutive models at the inlet and throat shear rates.

use std::error::Error;
use std::io;

use cfd_2d::solvers::venturi_flow::{BernoulliVenturi, VenturiGeometry, ViscousVenturi};
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};

/// Eight rounding operations bound the analytical continuity calculation.
const CONTINUITY_TOLERANCE: f64 = 8.0 * f64::EPSILON;

fn require(condition: bool, message: impl Into<String>) -> Result<(), Box<dyn Error>> {
    if condition {
        Ok(())
    } else {
        Err(io::Error::other(message.into()).into())
    }
}

fn validate_bernoulli_continuity() -> Result<(), Box<dyn Error>> {
    let geometry = VenturiGeometry::<f64>::iso_5167_standard();
    let inlet_velocity = 0.1;
    let inlet_pressure = 101_325.0;
    let density = 1000.0;
    let venturi = BernoulliVenturi::new(geometry.clone(), inlet_velocity, inlet_pressure, density);
    let throat_velocity = venturi.velocity_throat();
    let inlet_flow = geometry.area_inlet() * inlet_velocity;
    let throat_flow = geometry.area_throat() * throat_velocity;
    let relative_continuity_error = (throat_flow - inlet_flow).abs() / inlet_flow;

    require(
        relative_continuity_error <= CONTINUITY_TOLERANCE,
        format!(
            "Bernoulli continuity error {relative_continuity_error} exceeds {CONTINUITY_TOLERANCE}"
        ),
    )?;

    println!("\nBernoulli Venturi model");
    println!("  area ratio:          {:.6}", geometry.area_ratio());
    println!("  throat velocity:     {:.6} m/s", throat_velocity);
    println!("  throat pressure:     {:.3} Pa", venturi.pressure_throat());
    println!("  continuity error:    {:.3e}", relative_continuity_error);

    Ok(())
}

fn validate_viscous_loss() -> Result<(), Box<dyn Error>> {
    let geometry = VenturiGeometry::<f64>::iso_5167_standard();
    let inlet_pressure = 101_325.0;
    let loss_coefficient = 0.15;
    let venturi = ViscousVenturi::new(geometry, 0.5, inlet_pressure, 1000.0, loss_coefficient);
    let outlet_pressure = venturi.pressure_outlet_with_loss();
    let recovery_coefficient = venturi.pressure_recovery_coefficient();

    require(
        outlet_pressure.is_finite() && outlet_pressure < inlet_pressure,
        format!("viscous model produced invalid outlet pressure {outlet_pressure}"),
    )?;
    require(
        recovery_coefficient.is_finite() && recovery_coefficient < 0.0,
        format!("viscous model produced invalid recovery coefficient {recovery_coefficient}"),
    )?;

    println!("\nViscous-loss Venturi model");
    println!("  loss coefficient:    {:.3}", loss_coefficient);
    println!("  outlet pressure:     {:.3} Pa", outlet_pressure);
    println!("  recovery coefficient:{:.6}", recovery_coefficient);

    Ok(())
}

fn validate_blood_rheology() -> Result<(), Box<dyn Error>> {
    let geometry = VenturiGeometry::<f64>::iso_5167_standard();
    let inlet_velocity = 0.2;
    let throat_velocity = inlet_velocity / geometry.area_ratio();
    let inlet_shear_rate = 4.0 * inlet_velocity / geometry.w_inlet;
    let throat_shear_rate = 4.0 * throat_velocity / geometry.w_throat;
    let casson = CassonBlood::<f64>::normal_blood();
    let carreau_yasuda = CarreauYasudaBlood::<f64>::normal_blood();
    let casson_inlet_viscosity = casson.apparent_viscosity(inlet_shear_rate);
    let casson_throat_viscosity = casson.apparent_viscosity(throat_shear_rate);
    let carreau_inlet_viscosity = carreau_yasuda.apparent_viscosity(inlet_shear_rate);
    let carreau_throat_viscosity = carreau_yasuda.apparent_viscosity(throat_shear_rate);

    require(
        throat_shear_rate > inlet_shear_rate,
        format!(
            "Venturi throat shear rate {throat_shear_rate} does not exceed inlet {inlet_shear_rate}"
        ),
    )?;
    require(
        casson_throat_viscosity < casson_inlet_viscosity,
        format!(
            "Casson viscosity must decrease from inlet {casson_inlet_viscosity} to throat {casson_throat_viscosity}"
        ),
    )?;
    require(
        carreau_throat_viscosity < carreau_inlet_viscosity,
        format!(
            "Carreau-Yasuda viscosity must decrease from inlet {carreau_inlet_viscosity} to throat {carreau_throat_viscosity}"
        ),
    )?;

    println!("\nBlood rheology at the Venturi throat");
    println!("  inlet shear rate:    {:.3e} 1/s", inlet_shear_rate);
    println!("  throat shear rate:   {:.3e} 1/s", throat_shear_rate);
    println!(
        "  Casson viscosity:    {:.6e} → {:.6e} Pa·s",
        casson_inlet_viscosity, casson_throat_viscosity
    );
    println!(
        "  Carreau-Yasuda:      {:.6e} → {:.6e} Pa·s",
        carreau_inlet_viscosity, carreau_throat_viscosity
    );

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    println!("Venturi flow validation");
    validate_bernoulli_continuity()?;
    validate_viscous_loss()?;
    validate_blood_rheology()?;
    println!("\nAll Venturi calculations passed.");
    Ok(())
}
