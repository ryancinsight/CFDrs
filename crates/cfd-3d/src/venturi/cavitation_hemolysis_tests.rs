//! Comprehensive validation tests for cavitation and hemolysis in venturi systems.
//!
//! This test suite validates the integration of:
//! - Cavitation regime classification (stable vs. inertial)
//! - Hemolysis prediction in cavitating blood flow
//! - Sonoluminescence estimation
//! - Material damage prediction

#[cfg(test)]
mod tests {
    use cfd_core::physics::cavitation::{
        CavitationRegime, CavitationRegimeClassifier, RayleighPlesset,
    };
    use cfd_core::physics::hemolysis::{HemolysisCalculator, HemolysisModel};

    /// Test cavitation regime identification in venturi throat
    #[test]
    fn test_venturi_cavitation_regime_classification() {
        // Venturi throat conditions: high velocity, low pressure
        let ambient_pressure = 101325.0; // 1 atm
        let throat_pressure = 2000.0; // Severe cavitation
        let vapor_pressure = 2339.0; // Water at 20°C

        let bubble_model = RayleighPlesset {
            initial_radius: 10e-6, // 10 μm nuclei
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure,
            polytropic_index: 1.4,
        };

        // Simulate hydrodynamic cavitation (no acoustic field)
        let classifier =
            CavitationRegimeClassifier::new(bubble_model, throat_pressure, None, None);

        let analysis = classifier.analyze().unwrap();

        // In venturi throat with P << P_v, expect inertial cavitation
        assert_eq!(
            analysis.regime,
            CavitationRegime::Inertial,
            "Venturi throat with P={} Pa should produce inertial cavitation",
            throat_pressure
        );

        // Blake threshold should be exceeded
        assert!(
            throat_pressure < analysis.blake_threshold,
            "Pressure below Blake threshold indicates cavitation inception"
        );

        // High damage potential expected
        assert!(
            analysis.damage_potential > 0.8,
            "Inertial cavitation should have high damage potential, got {}",
            analysis.damage_potential
        );

        // Sonoluminescence likely
        assert!(
            analysis.sonoluminescence_probability > 0.8,
            "Inertial cavitation should produce sonoluminescence"
        );

        println!("{}", analysis);
    }

    /// Test stable vs. inertial transition
    #[test]
    fn test_stable_to_inertial_transition() {
        let vapor_pressure = 2339.0;
        let bubble_model = RayleighPlesset {
            initial_radius: 5e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure,
            polytropic_index: 1.4,
        };

        // Test pressure sweep
        let pressures = vec![80e3, 50e3, 20e3, 5e3, 2e3]; // Pa

        let mut regimes = Vec::new();
        for &pressure in &pressures {
            let classifier = CavitationRegimeClassifier::new(bubble_model, pressure, None, None);
            regimes.push(classifier.classify_regime());
        }

        // Should transition from None -> Stable -> Inertial as pressure decreases
        println!("Pressure sweep results:");
        for (i, pressure) in pressures.iter().enumerate() {
            println!("  P = {} Pa -> {:?}", pressure, regimes[i]);
        }

        // At very low pressure, must be inertial
        assert_eq!(regimes[4], CavitationRegime::Inertial);
    }

    /// Test hemolysis prediction in cavitating venturi flow
    #[test]
    fn test_hemolysis_in_cavitating_venturi() {
        // Blood properties
        let hematocrit = 0.45;
        let hemoglobin_initial = 15.0; // g/dL
        let blood_volume = 1e-3; // 1 mL
        let flow_rate = 1e-6; // 1 mL/s

        // Use Giersiepen power law model
        let model = HemolysisModel::giersiepen_standard();

        let calc = HemolysisCalculator::new(
            model,
            hematocrit,
            hemoglobin_initial,
            blood_volume,
            flow_rate,
        )
        .unwrap();

        // Venturi throat conditions
        // Typical microfluidic venturi: D_throat ~ 50 μm, velocity ~ 10 m/s
        // Shear stress τ = μ du/dy ≈ μ U / (D/2) for wall-bounded flow
        let mu = 0.0035; // Pa·s, blood viscosity at high shear
        let velocity = 10.0; // m/s
        let d_throat = 50e-6; // 50 μm

        // Wall shear stress estimate
        let shear_stress = 4.0 * mu * velocity / d_throat; // Factor of 4 for parabolic profile
        // τ ≈ 2800 Pa - extreme shear!

        // Exposure time: residence time in throat
        let throat_length = 200e-6; // 200 μm
        let exposure_time = throat_length / velocity; // 20 μs

        let damage = calc
            .model
            .damage_index(shear_stress, exposure_time)
            .unwrap();

        println!("\nVenturi Hemolysis Analysis:");
        println!("  Throat diameter: {} μm", d_throat * 1e6);
        println!("  Velocity: {} m/s", velocity);
        println!("  Shear stress: {:.1} Pa", shear_stress);
        println!("  Exposure time: {:.1} μs", exposure_time * 1e6);
        println!("  Damage index: {:.2e}", damage);

        // High shear stress should cause significant damage
        assert!(
            damage > 1e-6,
            "Extreme shear stress ({}Pa) should cause measurable hemolysis",
            shear_stress
        );

        // Calculate hemoglobin release
        let delta_hb = calc.hemoglobin_release(damage);
        let nih = calc.normalized_index(delta_hb);

        println!("  ΔHb: {:.3} mg/dL", delta_hb * 100.0); // Convert g/dL to mg/dL
        println!("  NIH: {:.3}%", nih);

        // FDA guidance: NIH < 0.01 g/100L acceptable
        if delta_hb * 100.0 > 10.0 {
            println!("  ⚠️  WARNING: Hemolysis exceeds FDA guidance (10 mg/dL)");
        }
    }

    /// Test hemolysis with cavitation enhancement
    #[test]
    fn test_cavitation_enhanced_hemolysis() {
        // When cavitation occurs, hemolysis is dramatically increased
        // due to microjet impact and shock waves

        let model = HemolysisModel::giersiepen_standard();
        let calc = HemolysisCalculator::new(model, 0.45, 15.0, 1e-3, 1e-6).unwrap();

        // Baseline: high shear without cavitation
        let baseline_stress = 200.0; // Pa, high but subcavitation
        let baseline_time = 0.1; // s
        let baseline_damage = calc
            .model
            .damage_index(baseline_stress, baseline_time)
            .unwrap();

        // With cavitation: add collapse pressure contribution
        let bubble_model = RayleighPlesset {
            initial_radius: 10e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        // Estimate collapse impact pressure
        let collapse_radius = 1e-6; // Collapse to 1 μm
        let liquid_density = 998.0;
        let sound_speed = 1481.0; // m/s in water
        let impact_distance = 5e-6; // 5 μm from wall

        let impact_pressure =
            liquid_density * sound_speed * sound_speed * collapse_radius / impact_distance;
        // P ≈ 440 MPa for strong collapse!

        // Effective stress with cavitation
        let effective_stress = baseline_stress + impact_pressure * 0.1; // 10% coupling
        let cavitation_damage = calc
            .model
            .damage_index(effective_stress, baseline_time)
            .unwrap();

        println!("\nCavitation-Enhanced Hemolysis:");
        println!("  Baseline damage: {:.2e}", baseline_damage);
        println!("  Collapse impact pressure: {:.1} MPa", impact_pressure / 1e6);
        println!("  Cavitation-enhanced damage: {:.2e}", cavitation_damage);
        println!(
            "  Enhancement factor: {:.1}×",
            cavitation_damage / baseline_damage
        );

        // Cavitation should dramatically increase hemolysis
        assert!(
            cavitation_damage > baseline_damage * 10.0,
            "Cavitation should amplify hemolysis by at least 10×"
        );
    }

    /// Test sonoluminescence energy estimation in venturi
    #[test]
    fn test_sonoluminescence_in_venturi() {
        let bubble_model = RayleighPlesset {
            initial_radius: 50e-6, // Larger bubble for SBSL
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        // Venturi conditions causing strong collapse
        let ambient_pressure = 101325.0; // Pa
        let ambient_temperature = 293.15; // K
        let collapse_radius = 1e-6; // Collapse to 1 μm
        let emissivity = 1.0; // Blackbody approximation
        let flash_duration = 50e-12; // 50 ps, typical for SBSL

        let estimate = bubble_model
            .estimate_sonoluminescence(
                ambient_pressure,
                ambient_temperature,
                collapse_radius,
                emissivity,
                flash_duration,
            )
            .unwrap();

        println!("\nSonoluminescence Estimation:");
        println!("  Initial radius: {} μm", bubble_model.initial_radius * 1e6);
        println!("  Collapse radius: {} μm", collapse_radius * 1e6);
        println!(
            "  Compression ratio: {:.0}×",
            bubble_model.initial_radius / collapse_radius
        );
        println!("  Peak temperature: {:.0} K", estimate.peak_temperature);
        println!(
            "  Peak pressure: {:.1} GPa",
            estimate.peak_pressure / 1e9
        );
        println!("  Radiated energy: {:.2e} J", estimate.radiated_energy);
        println!(
            "  Radiated power: {:.2e} W",
            estimate.radiated_energy / flash_duration
        );

        // Sanity checks
        assert!(
            estimate.peak_temperature > ambient_temperature,
            "Peak temperature should exceed ambient"
        );
        assert!(
            estimate.peak_pressure > ambient_pressure,
            "Peak pressure should exceed ambient"
        );
        assert!(
            estimate.radiated_energy > 0.0,
            "Must radiate some energy"
        );

        // For strong collapse (50× compression), expect temperatures > 1000 K
        assert!(
            estimate.peak_temperature > 1000.0,
            "Strong collapse should produce high temperatures"
        );
    }

    /// Test scaling laws for cavitation damage
    #[test]
    fn test_cavitation_damage_scaling() {
        // Test that damage scales appropriately with flow conditions

        let bubble_model = RayleighPlesset {
            initial_radius: 10e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        // Test at different velocities (pressure scales as U²)
        let velocities = vec![5.0, 10.0, 20.0]; // m/s
        let d_throat = 100e-6; // 100 μm

        println!("\nCavitation Damage Scaling:");
        for &velocity in &velocities {
            // Dynamic pressure drop scales as 0.5 * ρ * U²
            let dynamic_pressure = 0.5 * 998.0 * velocity * velocity;
            let throat_pressure = (101325.0 - dynamic_pressure).max(2339.0);

            let classifier =
                CavitationRegimeClassifier::new(bubble_model, throat_pressure, None, None);
            let analysis = classifier.analyze().unwrap();

            println!("  U = {} m/s:", velocity);
            println!("    P_throat = {:.1} kPa", throat_pressure / 1e3);
            println!("    Regime: {:?}", analysis.regime);
            println!("    Damage potential: {:.2}", analysis.damage_potential);
            println!("    Hemolysis risk: {:.2}", analysis.hemolysis_risk);
        }

        // Higher velocity should produce more cavitation damage
        let low_v_classifier = CavitationRegimeClassifier::new(
            bubble_model,
            101325.0 - 0.5 * 998.0 * 5.0 * 5.0,
            None,
            None,
        );
        let high_v_classifier = CavitationRegimeClassifier::new(
            bubble_model,
            (101325.0 - 0.5 * 998.0 * 20.0 * 20.0).max(2339.0),
            None,
            None,
        );

        assert!(
            high_v_classifier.damage_potential() > low_v_classifier.damage_potential(),
            "Higher velocity should produce more damage"
        );
    }

    /// Integration test: full venturi device assessment
    #[test]
    fn test_venturi_device_blood_trauma_assessment() {
        use cfd_core::physics::hemolysis::BloodTrauma;

        // Simulate blood flow through millifluidic venturi
        // Inlet: D=2mm, U=1m/s
        // Throat: D=0.5mm, U=16m/s (continuity)
        // Outlet: D=2mm, U=1m/s

        let inlet_diameter = 2e-3;
        let throat_diameter = 0.5e-3;
        let inlet_velocity = 1.0;

        // Continuity: A₁U₁ = A₂U₂
        let area_ratio = (inlet_diameter / throat_diameter).powi(2);
        let throat_velocity = inlet_velocity * area_ratio; // 16 m/s

        // Pressure drop (Bernoulli): ΔP = 0.5ρ(U₂² - U₁²)
        let rho = 1060.0; // blood density
        let pressure_drop = 0.5 * rho * (throat_velocity.powi(2) - inlet_velocity.powi(2));
        let throat_pressure = 101325.0 - pressure_drop; // ~ -27 kPa gauge

        println!("\nVenturi Device Assessment:");
        println!("  Inlet: D={} mm, U={} m/s", inlet_diameter * 1e3, inlet_velocity);
        println!("  Throat: D={} mm, U={:.1} m/s", throat_diameter * 1e3, throat_velocity);
        println!("  Pressure drop: {:.1} kPa", pressure_drop / 1e3);
        println!("  Throat pressure: {:.1} kPa", throat_pressure / 1e3);

        // Cavitation analysis
        let bubble_model = RayleighPlesset {
            initial_radius: 10e-6,
            liquid_density: 1060.0,
            liquid_viscosity: 0.0035,
            surface_tension: 0.056,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        let classifier = CavitationRegimeClassifier::new(bubble_model, throat_pressure, None, None);
        let cav_analysis = classifier.analyze().unwrap();

        println!("\n{}", cav_analysis);

        // Hemolysis analysis
        let model = HemolysisModel::giersiepen_standard();
        let calc = HemolysisCalculator::new(model, 0.45, 15.0, 1e-3, 1e-6).unwrap();

        // Shear stress in throat
        let mu = 0.0035;
        let shear_stress = 8.0 * mu * throat_velocity / throat_diameter;
        let throat_length = 1e-3; // 1 mm
        let exposure_time = throat_length / throat_velocity;

        let damage = calc.model.damage_index(shear_stress, exposure_time).unwrap();
        let delta_hb = calc.hemolysis_release(damage);

        let trauma = BloodTrauma {
            hemolysis_level: delta_hb * 100.0, // Convert to mg/dL
            platelet_activation: cav_analysis.hemolysis_risk * 100.0,
            thrombosis_risk: cav_analysis.damage_potential,
            max_shear_stress: shear_stress,
            avg_exposure_time: exposure_time,
        };

        println!("\n{}", trauma);

        // Device acceptance criteria
        if trauma.meets_fda_guidance() {
            println!("\n✓ Device meets FDA guidance for blood-contacting devices");
        } else {
            println!("\n✗ Device FAILS FDA guidance - redesign required");
            println!("   Recommendations:");
            println!("   - Increase throat diameter to reduce velocity");
            println!("   - Optimize diverging angle to prevent cavitation");
            println!("   - Add flow conditioning upstream");
        }
    }
}
