#![allow(missing_docs)]
#![allow(clippy::print_stdout)]
//! Validation tests for cavitation-VOF solver integration
//!
//! Tests cover:
//! - Cavitation inception and void fraction evolution
//! - Mass transfer rate calculations
//! - Damage accumulation
//! - Bubble dynamics integration
//! - Conservation properties

use aequitas::systems::si::quantities::{
    Dimensionless, Length, MassDensity, NumberDensity, Pressure, SurfaceTension, Time, Velocity,
};
use cfd_3d::vof::{
    AdvectionMethod, BubbleDynamicsConfig, CavitationVofConfig, CavitationVofSolver,
    InterfaceReconstruction, VofConfig,
};
use cfd_core::physics::cavitation::{damage::CavitationDamage, models::CavitationModel};
use cfd_core::physics::fluid::BloodModel;
use leto::{geometry::Vector3, Array2};

fn length(value: f64) -> Length<f64> {
    Length::from_base(value)
}

fn mass_density(value: f64) -> MassDensity<f64> {
    MassDensity::from_base(value)
}

fn number_density(value: f64) -> NumberDensity<f64> {
    NumberDensity::from_base(value)
}

fn pressure(value: f64) -> Pressure<f64> {
    Pressure::from_base(value)
}

fn surface_tension(value: f64) -> SurfaceTension<f64> {
    SurfaceTension::from_base(value)
}

fn time(value: f64) -> Time<f64> {
    Time::from_base(value)
}

fn velocity(value: f64) -> Velocity<f64> {
    Velocity::from_base(value)
}

#[test]
fn test_cavitation_inception() {
    println!("Testing cavitation inception...");

    // Create solver with low inception threshold
    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: surface_tension(0.072),
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::Kunz {
            vaporization_coeff: 100.0,
            condensation_coeff: 100.0,
        },
        damage_model: None,
        bubble_dynamics: None,
        inception_threshold: 0.5, // Low threshold
        max_void_fraction: 0.8,
        relaxation_time: time(1e-6),
        vapor_pressure: pressure(2330.0),
        liquid_density: mass_density(998.0),
        liquid_blood_model: BloodModel::Newtonian(1.002e-3),
        vapor_density: mass_density(0.023),
        sound_speed: velocity(1500.0),
        nuclei_transport: None,
    };

    let mut solver = CavitationVofSolver::new(10, 10, 10, config).expect("expected value");

    // Create low pressure field (should trigger cavitation)
    let velocity_field = vec![Vector3::zeros(); 1000];
    let pressure_field = Array2::from_elem([10, 100], 2000.0); // Below vapor pressure
    let density_field = Array2::from_elem([10, 100], 998.0);

    // Step simulation
    solver
        .step(time(1e-5), &velocity_field, &pressure_field, &density_field)
        .expect("expected value");

    // Check that cavitation was detected
    let stats = solver.cavitation_statistics();
    assert!(
        stats.cavitating_cells > 0,
        "Cavitation should be detected in low pressure regions"
    );

    println!("✓ Cavitation inception test passed");
}

#[test]
fn test_damage_accumulation() {
    println!("Testing cavitation damage accumulation...");

    let damage_model = CavitationDamage {
        yield_strength: pressure(200e6),
        ultimate_strength: pressure(500e6),
        hardness: pressure(200e6),
        fatigue_strength: pressure(150e6),
        cycles: 0,
    };

    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: surface_tension(0.072),
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::ZGB {
            nucleation_fraction: Dimensionless::from_base(5e-4),
            bubble_radius: length(1e-6),
            f_vap: 50.0,
            f_cond: 0.0,
        },
        damage_model: Some(damage_model),
        bubble_dynamics: Some(BubbleDynamicsConfig {
            initial_radius: length(1e-4),
            number_density: number_density(1e13),
            polytropic_exponent: 1.4,
            surface_tension: surface_tension(0.072),
        }),
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: time(0.1),
        vapor_pressure: pressure(2330.0),
        liquid_density: mass_density(998.0),
        liquid_blood_model: BloodModel::Newtonian(1.002e-3),
        vapor_density: mass_density(0.023),
        sound_speed: velocity(1500.0),
        nuclei_transport: None,
    };

    let mut solver = CavitationVofSolver::new(10, 10, 10, config).expect("expected value");

    // Create cavitating conditions
    let velocity_field = vec![Vector3::new(20.0, 0.0, 0.0); 1000]; // High velocity
    let density_field = Array2::from_elem([10, 100], 998.0);

    // Initialize void fraction to ensure damage calculation is triggered
    {
        let mut volume_fraction = solver.volume_fraction();
        volume_fraction.fill(0.1);
        solver
            .set_volume_fraction(&volume_fraction)
            .expect("expected value");
    }

    // Run multiple steps to accumulate damage
    // We oscillate pressure to induce both growth (low pressure) and collapse (high pressure)
    // Damage physically occurs during collapse events
    for i in 0..20 {
        let p_val = if i % 2 == 0 { 1000.0 } else { 50000.0 };
        let pressure_field = Array2::from_elem([10, 100], p_val);
        solver
            .step(time(1e-5), &velocity_field, &pressure_field, &density_field)
            .expect("expected value");
    }

    // Check damage accumulation
    let stats = solver.cavitation_statistics();
    assert!(
        stats.max_damage > 0.0,
        "Damage should accumulate in cavitating flow"
    );

    if let Some(damage_field) = solver.damage_field() {
        let total_damage: f64 = damage_field.iter().sum();
        assert!(total_damage > 0.0, "Total damage should be positive");
    }

    println!("✓ Damage accumulation test passed");
}

#[test]
fn test_sonoluminescence_energy_field_requires_collapse_and_is_finite() {
    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: surface_tension(0.072),
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::Kunz {
            vaporization_coeff: 100.0,
            condensation_coeff: 100.0,
        },
        damage_model: None,
        bubble_dynamics: Some(BubbleDynamicsConfig {
            initial_radius: length(1e-6),
            number_density: number_density(1e13),
            polytropic_exponent: 1.4,
            surface_tension: surface_tension(0.072),
        }),
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: time(1e-6),
        vapor_pressure: pressure(2330.0),
        liquid_density: mass_density(998.0),
        liquid_blood_model: BloodModel::Newtonian(1.002e-3),
        vapor_density: mass_density(0.023),
        sound_speed: velocity(1500.0),
        nuclei_transport: None,
    };

    let mut solver = CavitationVofSolver::new(6, 4, 3, config).expect("expected value");

    let velocity_field = vec![Vector3::zeros(); 6 * 4 * 3];
    let pressure_field = Array2::from_elem([6, 4 * 3], 1.0e6);
    let density_field = Array2::from_elem([6, 4 * 3], 998.0);

    solver
        .step(
            time(1e-12),
            &velocity_field,
            &pressure_field,
            &density_field,
        )
        .expect("expected value");

    let energy = solver
        .sonoluminescence_energy_field(&pressure_field, 293.15, 50e-12, 1.0)
        .expect("expected value");

    assert_eq!(energy.shape(), [6, 4 * 3]);

    let mut any_positive = false;
    for e in energy.iter().copied() {
        assert!(e.is_finite());
        assert!(e >= 0.0);
        any_positive |= e > 0.0;
    }
    assert!(any_positive);
}

#[test]
fn test_mass_conservation() {
    println!("Testing mass conservation with cavitation...");

    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: surface_tension(0.072),
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::SchnerrSauer {
            bubble_density: number_density(1e13),
            initial_radius: length(1e-6),
        },
        damage_model: None,
        bubble_dynamics: None,
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: time(0.1),
        vapor_pressure: pressure(2330.0),
        liquid_density: mass_density(998.0),
        liquid_blood_model: BloodModel::Newtonian(1.002e-3),
        vapor_density: mass_density(0.023),
        sound_speed: velocity(1500.0),
        nuclei_transport: None,
    };

    let mut solver = CavitationVofSolver::new(20, 10, 10, config).expect("expected value");

    // Initialize with some void fraction
    {
        let mut volume_fraction = solver.volume_fraction();
        for i in 0..volume_fraction.shape()[0] {
            for j in 0..volume_fraction.shape()[1] {
                if i > 10 && i < 15 {
                    // Central region
                    volume_fraction[[i, j]] = 0.1; // 10% void fraction
                }
            }
        }
        solver
            .set_volume_fraction(&volume_fraction)
            .expect("expected value");
    }

    let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); 2000];
    let pressure_field = Array2::from_elem([20, 100], 2330.0); // Exactly vapor pressure
    let density_field = Array2::from_elem([20, 100], 998.0);

    // Measure initial total volume
    let initial_volume: f64 = solver.volume_fraction().iter().sum();

    // Run simulation steps
    for _ in 0..50 {
        solver
            .step(time(1e-5), &velocity_field, &pressure_field, &density_field)
            .expect("expected value");
    }

    // Check volume conservation (should be close, cavitation adds/removes mass)
    let final_volume: f64 = solver.volume_fraction().iter().sum();
    let volume_change = (final_volume - initial_volume).abs();

    // Volume change should be reasonable (not explosive growth)
    assert!(
        volume_change < initial_volume * 0.1,
        "Volume change too large: {}",
        volume_change
    );

    println!(
        "✓ Mass conservation test passed (volume change: {:.2e})",
        volume_change
    );
}

#[test]
fn test_bubble_dynamics_integration() {
    println!("Testing bubble dynamics integration...");

    let bubble_config = BubbleDynamicsConfig {
        initial_radius: length(2e-6), // Larger initial bubble
        number_density: number_density(1e12),
        polytropic_exponent: 1.4,
        surface_tension: surface_tension(0.072),
    };
    let bubble_initial_radius = bubble_config.initial_radius.into_base();

    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: surface_tension(0.072),
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::ZGB {
            nucleation_fraction: Dimensionless::from_base(5e-4),
            bubble_radius: length(2e-6),
            f_vap: 50.0,
            f_cond: 0.01,
        },
        damage_model: None,
        bubble_dynamics: Some(bubble_config),
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: time(0.1),
        vapor_pressure: pressure(2330.0),
        liquid_density: mass_density(998.0),
        liquid_blood_model: BloodModel::Newtonian(1.002e-3),
        vapor_density: mass_density(0.023),
        sound_speed: velocity(1500.0),
        nuclei_transport: None,
    };

    let mut solver = CavitationVofSolver::new(5, 5, 5, config).expect("expected value");

    // Test with varying pressure (should cause bubble oscillation)
    let velocity_field = vec![Vector3::zeros(); 125];
    let mut pressure_field = Array2::from_elem([5, 25], 50000.0);

    // Run simulation with pressure variations
    let mut last_radius = bubble_initial_radius;
    let mut radius_changed = false;

    for step in 0..20 {
        // Vary pressure sinusoidally
        let pressure_variation = 10000.0 * (step as f64 * 0.1).sin();
        for i in 0..pressure_field.shape()[0] {
            for j in 0..pressure_field.shape()[1] {
                pressure_field[[i, j]] = 50000.0 + pressure_variation;
            }
        }

        let density_field = Array2::from_elem([5, 25], 998.0);
        solver
            .step(time(1e-6), &velocity_field, &pressure_field, &density_field)
            .expect("expected value");

        // Check that bubble radii are updated
        if let Some(radius_field) = solver.bubble_radius_field() {
            let current_radius = radius_field[[0, 0]];
            if (current_radius - last_radius).abs() > 1e-12 {
                radius_changed = true;
            }
            last_radius = current_radius;
            assert!(
                current_radius.is_finite() && current_radius >= 0.0,
                "Bubble radii should remain finite and non-negative"
            );
        }
    }

    assert!(
        radius_changed,
        "Bubble radius should change with pressure variation"
    );

    println!("✓ Bubble dynamics integration test passed");
}

#[test]
fn test_cavitation_statistics() {
    println!("Testing cavitation statistics calculation...");

    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: surface_tension(0.072),
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::Kunz {
            vaporization_coeff: 100.0,
            condensation_coeff: 100.0,
        },
        damage_model: None,
        bubble_dynamics: None,
        inception_threshold: 0.5,
        max_void_fraction: 0.8,
        relaxation_time: time(1e-6),
        vapor_pressure: pressure(2330.0),
        liquid_density: mass_density(998.0),
        liquid_blood_model: BloodModel::Newtonian(1.002e-3),
        vapor_density: mass_density(0.023),
        sound_speed: velocity(1500.0),
        nuclei_transport: None,
    };

    let mut solver = CavitationVofSolver::new(10, 10, 10, config).expect("expected value");

    // Create mixed conditions
    let velocity_field = vec![Vector3::zeros(); 1000];
    let mut pressure_field = Array2::zeros([10, 100]);

    // Set some cells to cavitating pressure
    for i in 0..10 {
        for j in 0..100 {
            if i < 3 {
                // First 3 rows cavitating
                pressure_field[[i, j]] = 1000.0; // Low pressure
            } else {
                pressure_field[[i, j]] = 101325.0; // Atmospheric
            }
        }
    }

    let density_field = Array2::from_elem([10, 100], 998.0);

    // Add some void fraction to cavitating regions
    {
        let mut volume_fraction = solver.volume_fraction();
        for i in 0..10 {
            for j in 0..100 {
                if i < 3 {
                    volume_fraction[[i, j]] = 0.2; // 20% void fraction
                }
            }
        }
        solver
            .set_volume_fraction(&volume_fraction)
            .expect("expected value");
    }

    solver
        .step(time(1e-5), &velocity_field, &pressure_field, &density_field)
        .expect("expected value");

    let stats = solver.cavitation_statistics();

    // Check statistics are reasonable
    assert!(stats.cavitation_fraction >= 0.0 && stats.cavitation_fraction <= 1.0);
    assert!(stats.total_void_fraction >= 0.0);
    assert!(stats.max_void_fraction >= 0.0 && stats.max_void_fraction <= 1.0);
    assert!(stats.cavitating_cells <= stats.total_cells);

    // Should have cavitation in the low-pressure regions
    assert!(
        stats.cavitating_cells > 0,
        "Should detect cavitation in low-pressure regions"
    );

    println!("✓ Cavitation statistics test passed");
    println!("   Cavitation fraction: {:.3}", stats.cavitation_fraction);
    println!("   Cavitating cells: {}", stats.cavitating_cells);
}

#[test]
fn test_cavitation_model_comparison() {
    println!("Testing cavitation model comparison...");

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
                bubble_density: number_density(1e13),
                initial_radius: length(1e-6),
            },
        ),
        (
            "ZGB",
            CavitationModel::ZGB {
                nucleation_fraction: Dimensionless::from_base(5e-4),
                bubble_radius: length(1e-6),
                f_vap: 50.0,
                f_cond: 0.01,
            },
        ),
    ];

    for (name, model) in models {
        let config = CavitationVofConfig {
            vof_config: VofConfig {
                surface_tension_coefficient: surface_tension(0.072),
                interface_compression: 0.1,
                reconstruction_method: InterfaceReconstruction::PLIC,
                advection_method: AdvectionMethod::Geometric,
                max_iterations: 10,
                tolerance: 1e-6,
                cfl_number: 0.3,
                enable_compression: false,
            },
            cavitation_model: model,
            damage_model: None,
            bubble_dynamics: None,
            inception_threshold: 0.3,
            max_void_fraction: 0.8,
            relaxation_time: time(1e-6),
            vapor_pressure: pressure(2330.0),
            liquid_density: mass_density(998.0),
            liquid_blood_model: BloodModel::Newtonian(1.002e-3),
            vapor_density: mass_density(0.023),
            sound_speed: velocity(1500.0),
            nuclei_transport: None,
        };

        let mut solver = CavitationVofSolver::new(5, 5, 5, config).expect("expected value");

        // Test with cavitating conditions
        let velocity_field = vec![Vector3::new(10.0, 0.0, 0.0); 125];
        let pressure_field = Array2::from_elem([5, 25], 1000.0); // Very low pressure
        let density_field = Array2::from_elem([5, 25], 998.0);

        solver
            .step(time(1e-5), &velocity_field, &pressure_field, &density_field)
            .expect("expected value");

        let stats = solver.cavitation_statistics();

        // All models should predict some cavitation under these conditions
        assert!(
            stats.cavitating_cells > 0,
            "{} model should predict cavitation under low pressure",
            name
        );

        println!(
            "✓ {} model: {} cavitating cells",
            name, stats.cavitating_cells
        );
    }

    println!("✓ Cavitation model comparison test passed");
}
