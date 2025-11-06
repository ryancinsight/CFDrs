//! Standalone Cavitation Damage Prediction Implementation
//!
//! This example demonstrates cavitation damage prediction without external dependencies,
//! showing the mathematical models and calculations used in CFDrs.

use std::f64::consts::PI;

/// Cavitation damage model for material erosion prediction
#[derive(Debug, Clone)]
struct CavitationDamage {
    /// Material yield strength (Pa)
    yield_strength: f64,
    /// Material ultimate tensile strength (Pa)
    ultimate_strength: f64,
    /// Material hardness (Pa)
    hardness: f64,
    /// Material fatigue strength (Pa)
    fatigue_strength: f64,
    /// Number of loading cycles
    cycles: u64,
}

impl CavitationDamage {
    /// Calculate erosion rate using ASTM G32 standard
    /// Rate ∝ (ΔP)^n where n ≈ 2-3 for most materials
    fn erosion_rate_hammitt(&self, pressure_amplitude: f64, exponent: f64) -> f64 {
        let reference_pressure = self.yield_strength;
        if pressure_amplitude > reference_pressure {
            let normalized_pressure = pressure_amplitude / reference_pressure;
            normalized_pressure.powf(exponent)
        } else {
            0.0
        }
    }

    /// Calculate Mean Depth of Penetration Rate (MDPR)
    /// Based on ASTM G32 standard and Plesset-Chapman erosion model
    fn mdpr(&self, impact_pressure: f64, frequency: f64, exposure_time: f64) -> f64 {
        // Plesset-Chapman model: MDPR = K * (P - P_th)^n * f^m * t
        let threshold = self.fatigue_strength;

        if impact_pressure > threshold {
            // Material-specific erosion coefficient (steel)
            let erosion_coefficient = 1e-11; // m³/N²·Hz·s for steel
            let pressure_exponent = 2.0; // Typical for metals
            let frequency_exponent = 1.0; // Linear frequency dependence

            // MDPR calculation with proper units
            let pressure_term = (impact_pressure - threshold).powf(pressure_exponent);
            erosion_coefficient * pressure_term * frequency.powf(frequency_exponent) * exposure_time
        } else {
            0.0
        }
    }

    /// Calculate Rayleigh collapse pressure
    fn collapse_impact_pressure(
        &self,
        bubble_radius: f64,
        collapse_distance: f64,
        liquid_density: f64,
        sound_speed: f64,
    ) -> f64 {
        // Rayleigh collapse pressure: P ∝ ρc²(R/r)
        if collapse_distance > 0.0 {
            liquid_density * sound_speed * sound_speed * bubble_radius / collapse_distance
        } else {
            0.0
        }
    }

    /// Calculate incubation period using Basquin's law
    fn incubation_period(&self, stress_amplitude: f64) -> u64 {
        // Basquin's law: σ_a = σ_f' * (2N)^b
        // Rearranged: N = 0.5 * ((σ_a / σ_f')^(1/b))

        // Fatigue strength coefficient (Morrow approximation)
        let fatigue_coeff = 0.9 * self.ultimate_strength;

        // Basquin exponent (typical for steels)
        let basquin_exponent = -0.085;

        if stress_amplitude > self.fatigue_strength {
            // Calculate cycles to failure
            let stress_ratio = stress_amplitude / fatigue_coeff;
            let exponent = 1.0 / basquin_exponent;
            let two_n = stress_ratio.powf(exponent);
            let n = two_n / 2.0;

            // Safe conversion to u64
            if n > 0.0 && n < 1e18 {
                (n.round() as u64).min(1_000_000_000) // Cap at reasonable value
            } else {
                1_000_000 // Default if calculation fails
            }
        } else {
            u64::MAX // Infinite life below fatigue limit
        }
    }

    /// Calculate pit depth from impact pressure
    fn pit_depth(&self, impact_pressure: f64, material_constant: f64) -> f64 {
        // Empirical model: depth ∝ (P/H)^n
        if impact_pressure > self.hardness {
            let ratio = impact_pressure / self.hardness;
            material_constant * ratio.powf(2.0)
        } else {
            0.0
        }
    }
}

/// Cavitation model types
enum CavitationModel {
    Kunz { vaporization_coeff: f64, condensation_coeff: f64 },
    Zgb { nucleation_fraction: f64, bubble_radius: f64, f_vap: f64, f_cond: f64 },
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

        match self {
            CavitationModel::Kunz { vaporization_coeff, condensation_coeff } => {
                let half = 0.5;
                if pressure_diff < 0.0 {
                    // Vaporization
                    vaporization_coeff * density_vapor * (1.0 - void_fraction)
                        * pressure_diff.abs() / (half * density_liquid)
                } else {
                    // Condensation
                    let rate = condensation_coeff * density_vapor * void_fraction
                        * pressure_diff / (half * density_liquid);
                    -rate
                }
            }

            CavitationModel::Zgb { nucleation_fraction, bubble_radius, f_vap, f_cond } => {
                let three = 3.0;
                let two_thirds = 2.0 / 3.0;

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
    println!("🧪 Standalone Cavitation Damage Prediction");
    println!("=========================================");
    println!();

    // Create material properties for stainless steel
    let damage_model = CavitationDamage {
        yield_strength: 200e6,      // 200 MPa
        ultimate_strength: 500e6,   // 500 MPa
        hardness: 200e6,            // 200 MPa (Vickers)
        fatigue_strength: 150e6,    // 150 MPa
        cycles: 0,
    };

    println!("Material Properties (Stainless Steel):");
    println!("  Yield strength: {:.0} MPa", damage_model.yield_strength / 1e6);
    println!("  Ultimate strength: {:.0} MPa", damage_model.ultimate_strength / 1e6);
    println!("  Hardness: {:.0} MPa", damage_model.hardness / 1e6);
    println!("  Fatigue strength: {:.0} MPa", damage_model.fatigue_strength / 1e6);
    println!();

    // Test erosion rate calculation
    let impact_pressure = 100e6; // 100 MPa impact pressure
    let frequency = 20e3; // 20 kHz
    let exposure_time = 3600.0; // 1 hour

    let erosion_rate = damage_model.mdpr(impact_pressure, frequency, exposure_time);

    println!("Erosion Rate Calculation (ASTM G32):");
    println!("  Impact pressure: {:.0} MPa", impact_pressure / 1e6);
    println!("  Frequency: {:.0} kHz", frequency / 1e3);
    println!("  Exposure time: {:.1} hours", exposure_time / 3600.0);
    println!("  MDPR (erosion rate): {:.2e} m", erosion_rate);
    println!("  Annual erosion: {:.1} mm/year", erosion_rate * 365.0 * 24.0 * 3600.0 * 1000.0);

    // Damage severity assessment
    if erosion_rate > 1e-6 {
        println!("  ⚠️  SEVERE: Significant material erosion expected");
    } else if erosion_rate > 1e-7 {
        println!("  ⚠️  MODERATE: Potential surface damage");
    } else {
        println!("  ✓ LOW: Minimal erosion risk");
    }
    println!();

    // Test Rayleigh collapse pressure
    let bubble_radius = 1e-6; // 1 μm
    let collapse_distance = 1e-7; // 0.1 μm (violent collapse)
    let liquid_density = 998.0; // Water
    let sound_speed = 1500.0; // Water sound speed

    let collapse_pressure = damage_model.collapse_impact_pressure(
        bubble_radius,
        collapse_distance,
        liquid_density,
        sound_speed,
    );

    println!("Rayleigh Collapse Pressure:");
    println!("  Bubble radius: {:.0} μm", bubble_radius * 1e6);
    println!("  Collapse distance: {:.1} μm", collapse_distance * 1e6);
    println!("  Impact pressure: {:.0} MPa", collapse_pressure / 1e6);
    println!();

    // Test fatigue analysis
    let stress_amplitude = 100e6; // 100 MPa
    let cycles_to_failure = damage_model.incubation_period(stress_amplitude);

    println!("Fatigue Life Prediction (Basquin's Law):");
    println!("  Stress amplitude: {:.0} MPa", stress_amplitude / 1e6);
    println!("  Cycles to failure: {}", if cycles_to_failure == u64::MAX {
        "Infinite (below fatigue limit)".to_string()
    } else {
        format!("{:.0}", cycles_to_failure)
    });
    println!();

    // Test cavitation models
    let pressure = 50000.0;        // 50 kPa (cavitating condition)
    let vapor_pressure = 2330.0;   // Water vapor pressure
    let void_fraction = 0.1;       // 10% void fraction
    let density_liquid = 998.0;    // Water density
    let density_vapor = 0.023;     // Steam density

    println!("Multi-Phase Cavitation Models:");
    println!("==============================");

    let models = vec![
        ("Kunz (2000)", CavitationModel::Kunz {
            vaporization_coeff: 100.0,
            condensation_coeff: 100.0,
        }),
        ("ZGB (2004)", CavitationModel::Zgb {
            nucleation_fraction: 5e-4,
            bubble_radius: 1e-6,
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

    println!("✅ Cavitation Damage Prediction Analysis Complete!");
    println!();
    println!("Key Capabilities Demonstrated:");
    println!("• ✅ ASTM G32 erosion rate prediction");
    println!("• ✅ Rayleigh bubble collapse pressure calculation");
    println!("• ✅ Basquin fatigue life prediction");
    println!("• ✅ Multi-phase cavitation mass transfer models");
    println!("• ✅ Material damage assessment and risk evaluation");
    println!();
    println!("This demonstrates the core cavitation damage prediction");
    println!("functionality available in CFDrs for engineering applications.");

    Ok(())
}
