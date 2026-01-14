//! Validation tests for cavitation-VOF solver integration
//!
//! Tests cover:
//! - Cavitation inception and void fraction evolution
//! - Mass transfer rate calculations
//! - Damage accumulation
//! - Bubble dynamics integration
//! - Conservation properties

use cfd_3d::vof::{
    AdvectionMethod, BubbleDynamicsConfig, CavitationVofConfig, CavitationVofSolver,
    InterfaceReconstruction, VofConfig,
};
use cfd_core::physics::cavitation::{damage::CavitationDamage, models::CavitationModel};
use nalgebra::{DMatrix, Vector3};

#[test]
fn test_cavitation_inception() {
    println!("Testing cavitation inception...");

    // Create solver with low inception threshold
    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: 0.072,
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
        relaxation_time: 1e-6,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        vapor_density: 0.023,
        sound_speed: 1500.0,
    };

    let mut solver = CavitationVofSolver::new(10, 10, 10, config).unwrap();

    // Create low pressure field (should trigger cavitation)
    let velocity_field = vec![Vector3::zeros(); 1000];
    let pressure_field = DMatrix::from_element(10, 100, 2000.0); // Below vapor pressure
    let density_field = DMatrix::from_element(10, 100, 998.0);

    // Step simulation
    solver
        .step(1e-5, &velocity_field, &pressure_field, &density_field)
        .unwrap();

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
        yield_strength: 200e6,
        ultimate_strength: 500e6,
        hardness: 200e6,
        fatigue_strength: 150e6,
        cycles: 0,
    };

    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: 0.072,
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::ZGB {
            nucleation_fraction: 5e-4,
            bubble_radius: 1e-6,
            f_vap: 50.0,
            f_cond: 0.01,
        },
        damage_model: Some(damage_model),
        bubble_dynamics: Some(BubbleDynamicsConfig {
            initial_radius: 1e-6,
            number_density: 1e13,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        }),
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: 0.1,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        vapor_density: 0.023,
        sound_speed: 1500.0,
    };

    let mut solver = CavitationVofSolver::new(10, 10, 10, config).unwrap();

    // Create cavitating conditions
    let velocity_field = vec![Vector3::new(20.0, 0.0, 0.0); 1000]; // High velocity
    let pressure_field = DMatrix::from_element(10, 100, 1000.0); // Low pressure
    let density_field = DMatrix::from_element(10, 100, 998.0);

    // Run multiple steps to accumulate damage
    for _ in 0..10 {
        solver
            .step(1e-5, &velocity_field, &pressure_field, &density_field)
            .unwrap();
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
            surface_tension_coefficient: 0.072,
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
            initial_radius: 1e-6,
            number_density: 1e13,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        }),
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: 1e-6,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        vapor_density: 0.023,
        sound_speed: 1500.0,
    };

    let mut solver = CavitationVofSolver::new(6, 4, 3, config).unwrap();

    let velocity_field = vec![Vector3::zeros(); 6 * 4 * 3];
    let pressure_field = DMatrix::from_element(6, 4 * 3, 1.0e6);
    let density_field = DMatrix::from_element(6, 4 * 3, 998.0);

    solver
        .step(1e-7, &velocity_field, &pressure_field, &density_field)
        .unwrap();

    let energy = solver
        .sonoluminescence_energy_field(&pressure_field, 293.15, 50e-12, 1.0)
        .unwrap();

    assert_eq!(energy.nrows(), 6);
    assert_eq!(energy.ncols(), 4 * 3);

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
            surface_tension_coefficient: 0.072,
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::SchnerrSauer {
            bubble_density: 1e13,
            initial_radius: 1e-6,
        },
        damage_model: None,
        bubble_dynamics: None,
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: 0.1,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        vapor_density: 0.023,
        sound_speed: 1500.0,
    };

    let mut solver = CavitationVofSolver::new(20, 10, 10, config).unwrap();

    // Initialize with some void fraction
    {
        let mut volume_fraction = solver.volume_fraction();
        for i in 0..volume_fraction.nrows() {
            for j in 0..volume_fraction.ncols() {
                if i > 10 && i < 15 {
                    // Central region
                    volume_fraction[(i, j)] = 0.1; // 10% void fraction
                }
            }
        }
        solver.set_volume_fraction(&volume_fraction).unwrap();
    }

    let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); 2000];
    let pressure_field = DMatrix::from_element(20, 100, 2330.0); // Exactly vapor pressure
    let density_field = DMatrix::from_element(20, 100, 998.0);

    // Measure initial total volume
    let initial_volume: f64 = solver.volume_fraction().iter().sum();

    // Run simulation steps
    for _ in 0..50 {
        solver
            .step(1e-5, &velocity_field, &pressure_field, &density_field)
            .unwrap();
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
        initial_radius: 2e-6, // Larger initial bubble
        number_density: 1e12,
        polytropic_exponent: 1.4,
        surface_tension: 0.072,
    };

    let config = CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: 0.072,
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::ZGB {
            nucleation_fraction: 5e-4,
            bubble_radius: 2e-6,
            f_vap: 50.0,
            f_cond: 0.01,
        },
        damage_model: None,
        bubble_dynamics: Some(bubble_config),
        inception_threshold: 0.3,
        max_void_fraction: 0.8,
        relaxation_time: 0.1,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        vapor_density: 0.023,
        sound_speed: 1500.0,
    };

    let mut solver = CavitationVofSolver::new(5, 5, 5, config).unwrap();

    // Test with varying pressure (should cause bubble oscillation)
    let velocity_field = vec![Vector3::zeros(); 125];
    let mut pressure_field = DMatrix::from_element(5, 25, 50000.0);

    // Run simulation with pressure variations
    let mut last_radius = 1e-6;
    let mut radius_changed = false;

    for step in 0..20 {
        // Vary pressure sinusoidally
        let pressure_variation = 10000.0 * (step as f64 * 0.1).sin();
        for i in 0..pressure_field.nrows() {
            for j in 0..pressure_field.ncols() {
                pressure_field[(i, j)] = 50000.0 + pressure_variation;
            }
        }

        let density_field = DMatrix::from_element(5, 25, 998.0);
        solver
            .step(1e-6, &velocity_field, &pressure_field, &density_field)
            .unwrap();

        // Check that bubble radii are updated
        if let Some(radius_field) = solver.bubble_radius_field() {
            let current_radius = radius_field[(0, 0)];
            if (current_radius - last_radius).abs() > 1e-12 {
                radius_changed = true;
            }
            last_radius = current_radius;
            assert!(current_radius > 0.0, "Bubble radii should be positive");
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
            surface_tension_coefficient: 0.072,
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
        damage_model: Some(CavitationDamage {
            yield_strength: 200e6,
            ultimate_strength: 500e6,
            hardness: 200e6,
            fatigue_strength: 150e6,
            cycles: 1000,
        }),
        bubble_dynamics: None,
        inception_threshold: 0.5,
        max_void_fraction: 0.8,
        relaxation_time: 1e-6,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        vapor_density: 0.023,
        sound_speed: 1500.0,
    };

    let mut solver = CavitationVofSolver::new(10, 10, 10, config).unwrap();

    // Create mixed conditions
    let velocity_field = vec![Vector3::zeros(); 1000];
    let mut pressure_field = DMatrix::zeros(10, 100);

    // Set some cells to cavitating pressure
    for i in 0..10 {
        for j in 0..100 {
            if i < 3 {
                // First 3 rows cavitating
                pressure_field[(i, j)] = 1000.0; // Low pressure
            } else {
                pressure_field[(i, j)] = 101325.0; // Atmospheric
            }
        }
    }

    let density_field = DMatrix::from_element(10, 100, 998.0);

    // Add some void fraction to cavitating regions
    {
        let mut volume_fraction = solver.volume_fraction();
        for i in 0..10 {
            for j in 0..100 {
                if i < 3 {
                    volume_fraction[(i, j)] = 0.2; // 20% void fraction
                }
            }
        }
        solver.set_volume_fraction(&volume_fraction).unwrap();
    }

    solver
        .step(1e-5, &velocity_field, &pressure_field, &density_field)
        .unwrap();

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
                bubble_density: 1e13,
                initial_radius: 1e-6,
            },
        ),
        (
            "ZGB",
            CavitationModel::ZGB {
                nucleation_fraction: 5e-4,
                bubble_radius: 1e-6,
                f_vap: 50.0,
                f_cond: 0.01,
            },
        ),
    ];

    for (name, model) in models {
        let config = CavitationVofConfig {
            vof_config: VofConfig {
                surface_tension_coefficient: 0.072,
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
            relaxation_time: 1e-6,
            vapor_pressure: 2330.0,
            liquid_density: 998.0,
            vapor_density: 0.023,
            sound_speed: 1500.0,
        };

        let mut solver = CavitationVofSolver::new(5, 5, 5, config).unwrap();

        // Test with cavitating conditions
        let velocity_field = vec![Vector3::new(10.0, 0.0, 0.0); 125];
        let pressure_field = DMatrix::from_element(5, 25, 1000.0); // Very low pressure
        let density_field = DMatrix::from_element(5, 25, 998.0);

        solver
            .step(1e-5, &velocity_field, &pressure_field, &density_field)
            .unwrap();

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
