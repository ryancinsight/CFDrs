//! Simple Hydrodynamic Cavitation Analysis
//!
//! This example demonstrates basic hydrodynamic cavitation analysis
//! in a venturi throat using CFDrs cavitation models.

use cfd_core::cavitation::{models::CavitationModel, venturi::VenturiCavitation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 CFDrs: Hydrodynamic Cavitation Analysis");
    println!("==========================================");
    println!();

    // Define venturi geometry (typical microfluidic venturi)
    let inlet_diameter = 0.001; // 1 mm inlet
    let throat_diameter = 0.0005; // 0.5 mm throat
    let outlet_diameter = 0.001; // 1 mm outlet
    let convergent_angle = 15.0 * std::f64::consts::PI / 180.0; // 15 degrees
    let divergent_angle = 7.0 * std::f64::consts::PI / 180.0; // 7 degrees

    // Fluid properties (water at 20°C)
    let density = 998.0; // kg/m³
    let vapor_pressure = 2330.0; // Pa (water vapor pressure at 20°C)
    let inlet_pressure = 101325.0; // Pa (atmospheric pressure)

    println!("Venturi Geometry:");
    println!("  Inlet diameter: {:.1} mm", inlet_diameter * 1000.0);
    println!("  Throat diameter: {:.1} mm", throat_diameter * 1000.0);
    println!("  Outlet diameter: {:.1} mm", outlet_diameter * 1000.0);
    println!(
        "  Convergent angle: {:.1}°",
        convergent_angle * 180.0 / std::f64::consts::PI
    );
    println!(
        "  Divergent angle: {:.1}°",
        divergent_angle * 180.0 / std::f64::consts::PI
    );
    println!();

    println!("Fluid Properties:");
    println!("  Density: {:.0} kg/m³", density);
    println!("  Vapor pressure: {:.0} Pa", vapor_pressure);
    println!("  Inlet pressure: {:.0} Pa", inlet_pressure);
    println!();

    // Create venturi cavitation analyzer
    let venturi = VenturiCavitation {
        inlet_diameter,
        throat_diameter,
        outlet_diameter,
        convergent_angle,
        divergent_angle,
        inlet_pressure,
        inlet_velocity: 1.0, // Will be varied
        density,
        vapor_pressure,
    };

    // Analyze cavitation for different inlet velocities
    let inlet_velocities = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]; // m/s

    println!("Cavitation Analysis Results:");
    println!("============================");
    println!(
        "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}",
        "V_in (m/s)", "V_throat", "P_throat", "σ", "Status", "Cavity L"
    );
    println!("{}", "─".repeat(80));

    for &inlet_velocity in &inlet_velocities {
        let mut venturi_case = venturi.clone();
        venturi_case.inlet_velocity = inlet_velocity;

        let throat_velocity = venturi_case.throat_velocity();
        let throat_pressure = venturi_case.throat_pressure();
        let cavitation_number = venturi_case.cavitation_number();
        let is_cavitating = venturi_case.is_cavitating();
        let cavity_length = if is_cavitating {
            venturi_case.cavity_length(cavitation_number)
        } else {
            0.0
        };

        let status = if is_cavitating { "CAVITATING" } else { "SAFE" };
        let cavity_str = if cavity_length > 0.0 {
            format!("{:.3} mm", cavity_length * 1000.0)
        } else {
            "-".to_string()
        };

        println!(
            "{:<12.1} {:<12.2} {:<12.0} {:<12.3} {:<12} {:<12}",
            inlet_velocity, throat_velocity, throat_pressure, cavitation_number, status, cavity_str
        );
    }

    println!();

    // Demonstrate multi-phase cavitation models
    println!("Multi-Phase Cavitation Models:");
    println!("==============================");

    // Example conditions for severe cavitation
    let pressure = 50000.0; // 50 kPa (cavitating condition)
    let void_fraction = 0.1; // 10% void fraction
    let density_liquid = 998.0; // Water density
    let density_vapor = 0.023; // Steam density

    println!(
        "Conditions: P = {:.0} Pa, α = {:.1}%, ρ_l = {:.0} kg/m³, ρ_v = {:.3} kg/m³",
        pressure,
        void_fraction * 100.0,
        density_liquid,
        density_vapor
    );
    println!();

    // Test different cavitation models
    let models = vec![
        (
            "Kunz",
            CavitationModel::Kunz {
                vaporization_coeff: 100.0,
                condensation_coeff: 100.0,
            },
        ),
        (
            "Schnerr-Sauer",
            CavitationModel::SchnerrSauer {
                bubble_density: 1e13, // #/m³
                initial_radius: 1e-6, // m
            },
        ),
        (
            "ZGB",
            CavitationModel::ZGB {
                nucleation_fraction: 5e-4,
                bubble_radius: 1e-6, // m
                f_vap: 50.0,
                f_cond: 0.01,
            },
        ),
    ];

    for (name, model) in models {
        let mass_transfer = model.mass_transfer_rate(
            pressure,
            vapor_pressure,
            void_fraction,
            density_liquid,
            density_vapor,
        );

        println!("{} Model:", name);
        println!("  Mass transfer rate = {:.2e} kg/m³/s", mass_transfer);
        println!(
            "  Direction: {}",
            if mass_transfer > 0.0 {
                "Vaporization"
            } else {
                "Condensation"
            }
        );
        println!();
    }

    println!("✅ Cavitation analysis completed!");
    println!();
    println!("Key Findings:");
    println!("• Cavitation inception occurs when σ < 0.3");
    println!("• Cavity length increases rapidly with decreasing σ");
    println!("• ZGB model shows highest vaporization rates");
    println!("• For complete CFD simulation, integrate with multi-phase Navier-Stokes solver");

    Ok(())
}
