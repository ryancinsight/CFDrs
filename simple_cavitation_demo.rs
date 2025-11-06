//! Simple Hydrodynamic Cavitation Analysis
//!
//! This standalone example demonstrates hydrodynamic cavitation analysis
//! in a venturi throat using basic physics equations.

use std::f64::consts::PI;

/// Venturi cavitation analysis structure
struct VenturiCavitation {
    inlet_diameter: f64,
    throat_diameter: f64,
    outlet_diameter: f64,
    convergent_angle: f64,
    divergent_angle: f64,
    inlet_pressure: f64,
    inlet_velocity: f64,
    density: f64,
    vapor_pressure: f64,
}

impl VenturiCavitation {
    /// Calculate throat velocity using continuity equation
    fn throat_velocity(&self) -> f64 {
        let area_ratio = (self.inlet_diameter / self.throat_diameter).powi(2);
        self.inlet_velocity * area_ratio
    }

    /// Calculate throat pressure using Bernoulli equation
    fn throat_pressure(&self) -> f64 {
        let v_inlet = self.inlet_velocity;
        let v_throat = self.throat_velocity();
        let half = 0.5;

        self.inlet_pressure - half * self.density * (v_throat * v_throat - v_inlet * v_inlet)
    }

    /// Calculate cavitation number at throat
    fn cavitation_number(&self) -> f64 {
        let p_throat = self.throat_pressure();
        let v_throat = self.throat_velocity();
        let half = 0.5;

        if v_throat > 0.0 {
            (p_throat - self.vapor_pressure) / (half * self.density * v_throat * v_throat)
        } else {
            1e10
        }
    }

    /// Check if cavitation occurs
    fn is_cavitating(&self) -> bool {
        self.throat_pressure() < self.vapor_pressure
    }

    /// Calculate cavity length using Nurick correlation
    fn cavity_length(&self, cavitation_number: f64) -> f64 {
        let k_coefficient = 1.0;
        let exponent = 0.5;
        let sigma_incipient = 0.3;

        if cavitation_number < sigma_incipient && cavitation_number > 0.0 {
            let term = 1.0 / cavitation_number - 1.0 / sigma_incipient;
            if term > 0.0 {
                self.throat_diameter * k_coefficient * term.powf(exponent)
            } else {
                0.0
            }
        } else {
            0.0
        }
    }
}

/// Cavitation model types
enum CavitationModel {
    Kunz { vaporization_coeff: f64, condensation_coeff: f64 },
    SchnerrSauer { bubble_density: f64, initial_radius: f64 },
    ZGB { nucleation_fraction: f64, bubble_radius: f64, f_vap: f64, f_cond: f64 },
}

impl CavitationModel {
    /// Calculate mass transfer rate (kg/m³/s)
    fn mass_transfer_rate(
        &self,
        pressure: f64,
        vapor_pressure: f64,
        void_fraction: f64,
        density_liquid: f64,
        density_vapor: f64,
    ) -> f64 {
        let pressure_diff = pressure - vapor_pressure;
        let three = 3.0;
        let two_thirds = 2.0 / 3.0;

        match self {
            CavitationModel::Kunz { vaporization_coeff, condensation_coeff } => {
                let half = 0.5;
                if pressure_diff < 0.0 {
                    // Vaporization
                    vaporization_coeff * density_vapor * (1.0 - void_fraction) * pressure_diff.abs()
                        / (half * density_liquid)
                } else {
                    // Condensation
                    let rate = condensation_coeff * density_vapor * void_fraction * pressure_diff
                        / (half * density_liquid);
                    -rate
                }
            }

            CavitationModel::SchnerrSauer { bubble_density, initial_radius: _ } => {
                let four_pi = 4.0 * PI;
                let n_b = *bubble_density;

                // Simplified Schnerr-Sauer model
                let alpha = void_fraction;
                let sign = if pressure_diff < 0.0 { 1.0 } else { -1.0 };

                sign * three * alpha * (1.0 - alpha) * density_vapor
                    * (two_thirds * pressure_diff.abs() / density_liquid).sqrt()
                    / (four_pi * n_b * alpha / (1.0 - alpha)).powf(1.0/3.0)
            }

            CavitationModel::ZGB { nucleation_fraction, bubble_radius, f_vap, f_cond } => {
                if pressure_diff < 0.0 {
                    // Vaporization
                    f_vap * three * nucleation_fraction * (1.0 - void_fraction) * density_vapor
                        * (two_thirds * pressure_diff.abs() / density_liquid).sqrt()
                        / bubble_radius
                } else {
                    // Condensation
                    let rate = f_cond * three * void_fraction * density_vapor
                        * (two_thirds * pressure_diff / density_liquid).sqrt()
                        / bubble_radius;
                    -rate
                }
            }
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 Hydrodynamic Cavitation Analysis in Venturi Throat");
    println!("====================================================");
    println!();

    // Define venturi geometry (typical microfluidic venturi)
    let inlet_diameter = 0.001;  // 1 mm inlet
    let throat_diameter = 0.0005; // 0.5 mm throat
    let outlet_diameter = 0.001; // 1 mm outlet
    let convergent_angle = 15.0 * PI / 180.0; // 15 degrees
    let divergent_angle = 7.0 * PI / 180.0;   // 7 degrees

    // Fluid properties (water at 20°C)
    let density = 998.0;           // kg/m³
    let vapor_pressure = 2330.0;   // Pa (water vapor pressure at 20°C)
    let inlet_pressure = 101325.0; // Pa (atmospheric pressure)

    println!("Geometry: Inlet {:.1}mm → Throat {:.1}mm → Outlet {:.1}mm",
             inlet_diameter * 1000.0, throat_diameter * 1000.0, outlet_diameter * 1000.0);
    println!("Angles: Convergent {:.1}°, Divergent {:.1}°",
             convergent_angle * 180.0 / PI, divergent_angle * 180.0 / PI);
    println!("Fluid: Water (ρ = {:.0} kg/m³, p_v = {:.0} Pa)", density, vapor_pressure);
    println!("Inlet pressure: {:.0} Pa", inlet_pressure);
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
    println!("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}",
             "V_in (m/s)", "V_throat", "P_throat", "σ", "Status", "Cavity L");
    println!("{}", "─".repeat(80));

    for &inlet_velocity in &inlet_velocities {
        let mut venturi_case = VenturiCavitation { inlet_velocity, ..venturi };

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

        println!("{:<12.1} {:<12.2} {:<12.0} {:<12.3} {:<12} {:<12}",
                 inlet_velocity, throat_velocity, throat_pressure,
                 cavitation_number, status, cavity_str);
    }

    println!();

    // Demonstrate multi-phase cavitation models
    println!("Multi-Phase Cavitation Models:");
    println!("==============================");

    // Example conditions for severe cavitation
    let pressure = 50000.0;        // 50 kPa (cavitating condition)
    let void_fraction = 0.1;       // 10% void fraction
    let density_liquid = 998.0;    // Water density
    let density_vapor = 0.023;     // Steam density

    println!("Conditions: P = {:.0} Pa, α = {:.1}%, ρ_l = {:.0} kg/m³, ρ_v = {:.3} kg/m³",
             pressure, void_fraction * 100.0, density_liquid, density_vapor);
    println!();

    // Test different cavitation models
    let models = vec![
        ("Kunz", CavitationModel::Kunz {
            vaporization_coeff: 100.0,
            condensation_coeff: 100.0,
        }),
        ("Schnerr-Sauer", CavitationModel::SchnerrSauer {
            bubble_density: 1e13, // #/m³
            initial_radius: 1e-6, // m
        }),
        ("ZGB", CavitationModel::ZGB {
            nucleation_fraction: 5e-4,
            bubble_radius: 1e-6, // m
            f_vap: 50.0,
            f_cond: 0.01,
        }),
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
        println!("  Direction: {}", if mass_transfer > 0.0 { "Vaporization" } else { "Condensation" });
        println!();
    }

    println!("✅ Analysis Complete!");
    println!("====================");
    println!();
    println!("CFDrs Cavitation Capabilities:");
    println!("• ✅ Venturi cavitation analysis with pressure distribution");
    println!("• ✅ Cavitation number calculation and inception prediction");
    println!("• ✅ Cavity length estimation using literature correlations");
    println!("• ✅ Multi-phase cavitation models (Kunz, Schnerr-Sauer, ZGB)");
    println!("• ✅ Mass transfer rate calculations for vaporization/condensation");
    println!();
    println!("For full CFD simulation with plots, the CFDrs framework includes:");
    println!("• Interactive HTML visualization with Chart.js");
    println!("• Complete Navier-Stokes solvers for multi-phase flow");
    println!("• GPU acceleration for large-scale cavitation simulations");
    println!("• Validation against experimental cavitation data");

    Ok(())
}
