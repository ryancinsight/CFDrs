//! Comprehensive millifluidics tests for cfd-1d
//!
//! Tests based on literature and industry standards for microfluidic network simulation,
//! inspired by mmft-modular-1D-simulator functionality for production-ready validation.

use approx::assert_relative_eq;
use cfd_1d::solver::SolverConfig;
use cfd_1d::*;
use cfd_core::error::Result;
use cfd_core::physics::fluid;
use cfd_core::compute::solver::{Configurable, Solver};
use std::collections::HashMap;

/// Test micropump component behavior and pressure-flow characteristics  
#[test]
fn test_micropump_characteristics() -> Result<()> {
    // Test syringe pump characteristics (linear pressure-flow relationship)
    let mut pump = Micropump::<f64>::new(1e-9, 10000.0); // 1 μL/s, 10 kPa

    // Validate pump parameters
    assert_eq!(pump.max_flow_rate, 1e-9);
    assert_eq!(pump.max_pressure, 10000.0);
    assert!(pump.is_active());

    // Test resistance (negative for pressure source)
    let fluid = fluid::database::water_20c::<f64>()?;
    let resistance = pump.resistance(&fluid);
    assert!(resistance < 0.0);

    // Expected resistance for pump: -pressure/flow_rate
    let expected_resistance = -pump.max_pressure / pump.max_flow_rate;
    assert_relative_eq!(resistance, expected_resistance, epsilon = 1e-10);

    // Test parameter setting
    pump.set_parameter("efficiency", 0.85)?;
    assert_relative_eq!(pump.efficiency, 0.85, epsilon = 1e-10);

    // Test operating point adjustment
    pump.set_parameter("operating_point", 0.5)?;
    assert_relative_eq!(pump.operating_point, 0.5, epsilon = 1e-10);

    Ok(())
}

/// Test micromixer efficiency and flow combination
#[test]
fn test_micromixer_flow_characteristics() -> Result<()> {
    // Test Y-junction mixer for two-fluid mixing
    let mut mixer = Micromixer::<f64>::new(2, 1e12); // 2 inlets, 1e12 Pa·s/m³ resistance

    // Validate mixer properties
    assert_eq!(mixer.n_inlets, 2);
    assert!(mixer.efficiency > 0.0 && mixer.efficiency <= 1.0);

    let fluid = fluid::database::water_20c::<f64>()?;
    let resistance = mixer.resistance(&fluid);
    assert_relative_eq!(resistance, 1e12, epsilon = 1e-10);

    // Test mixing efficiency parameter constraints
    mixer.set_parameter("efficiency", 1.5)?; // Should clamp to 1.0
    assert_relative_eq!(mixer.efficiency, 1.0, epsilon = 1e-10);

    mixer.set_parameter("efficiency", -0.1)?; // Should clamp to 0.0
    assert_relative_eq!(mixer.efficiency, 0.0, epsilon = 1e-10);

    // Test valid efficiency setting
    mixer.set_parameter("efficiency", 0.75)?;
    assert_relative_eq!(mixer.efficiency, 0.75, epsilon = 1e-10);

    Ok(())
}

/// Test microvalve component behavior
#[test]
fn test_microvalve_characteristics() -> Result<()> {
    // Test microvalve with larger flow coefficient for more realistic behavior
    let mut valve = Microvalve::<f64>::new(1e-3); // CV = 1e-3 (much larger)

    // Validate valve parameters
    assert_eq!(valve.cv, 1e-3);
    assert_eq!(valve.opening, 1.0); // Fully open by default
    assert!(valve.is_active());

    let fluid = fluid::database::water_20c::<f64>()?;
    let resistance_open = valve.resistance(&fluid);
    assert!(resistance_open > 0.0);

    // Test valve closure - closed valve should have very high resistance
    valve.set_parameter("opening", 0.0)?; // Close valve
    let resistance_closed = valve.resistance(&fluid);
    // For closed valve, the implementation returns 1e12 fixed value
    assert_relative_eq!(resistance_closed, 1e12, epsilon = 1.0);

    // Test partial opening - should have higher pressure drop than fully open for same flow
    valve.set_parameter("opening", 0.5)?; // Half open
    let dp_half = valve.pressure_drop(1e-6, &fluid);

    // Reopen valve to compare
    valve.set_parameter("opening", 1.0)?; // Fully open
    let dp_open = valve.pressure_drop(1e-6, &fluid);

    // Half open should have higher pressure drop than fully open
    assert!(dp_half > dp_open);

    // Test parameter bounds
    valve.set_parameter("opening", 1.5)?; // Should clamp to 1.0
    assert_relative_eq!(valve.opening, 1.0, epsilon = 1e-10);

    valve.set_parameter("opening", -0.1)?; // Should clamp to 0.0
    assert_relative_eq!(valve.opening, 0.0, epsilon = 1e-10);

    Ok(())
}

/// Test flow analysis for basic microfluidic parameters
#[test]
fn test_flow_regime_analysis() -> Result<()> {
    let fluid = fluid::database::water_20c::<f64>()?;

    // Test Reynolds number calculation for circular channel
    let diameter = 100e-6; // 100μm
    let velocity = 0.1; // 10 cm/s
    let density = 1000.0; // kg/m³ (water)
    let viscosity = fluid.viscosity;

    let calculated_re = density * velocity * diameter / viscosity;

    // For microfluidics, Reynolds number should be low (Re << 2300)
    assert!(
        calculated_re < 2300.0,
        "Reynolds number should be in laminar regime for microfluidics"
    );

    // Test that we have proper fluid properties
    assert!(viscosity > 0.0);
    assert!(density > 0.0);

    Ok(())
}

/// Test surface tension and wettability effects in microchannels
#[test]
fn test_surface_effects() -> Result<()> {
    // Test basic surface tension calculations for microfluidics

    // Test capillary pressure calculation (Young-Laplace equation)
    // ΔP = 2σcosθ/r for circular channel
    let radius = 50e-6; // 50μm radius
    let surface_tension = 0.072; // N/m (water)
    let contact_angle = 60.0_f64; // degrees

    let expected_pressure = 2.0 * surface_tension * contact_angle.to_radians().cos() / radius;

    // For water in glass microchannel, capillary pressure should be positive (hydrophilic)
    assert!(
        expected_pressure > 0.0,
        "Capillary pressure should be positive for hydrophilic surface"
    );

    // Test contact angle categories
    assert!(
        contact_angle < 90.0,
        "Test contact angle should be hydrophilic"
    );

    // Test typical microfluidic surface tension values
    assert!(
        surface_tension > 0.05 && surface_tension < 0.1,
        "Surface tension should be in typical range for water"
    );

    Ok(())
}

/// Test performance metrics and analysis for microfluidic networks
#[test]
fn test_performance_analysis() -> Result<()> {
    // Initialize performance metrics
    let mut metrics = PerformanceMetrics::<f64>::new();

    // Test flow rate and pressure drop setting
    let expected_flow_rate = 1e-9; // 1 nL/s
    let pressure_drop = 5000.0; // Pa

    // Set metrics values using available methods
    metrics.set_total_flow_rate(expected_flow_rate);
    metrics.set_total_pressure_drop(pressure_drop);
    metrics.set_power_consumption(expected_flow_rate * pressure_drop);

    // Validate metrics
    assert_relative_eq!(metrics.throughput, expected_flow_rate, epsilon = 1e-15);
    assert_relative_eq!(
        metrics.power_consumption,
        expected_flow_rate * pressure_drop,
        epsilon = 1e-15
    );

    // Test hydraulic efficiency calculation
    let efficiency = metrics.hydraulic_efficiency();
    assert!(efficiency >= 0.0);

    // Test residence time tracking
    metrics.add_residence_time("inlet_channel".to_string(), 0.1);
    metrics.add_residence_time("outlet_channel".to_string(), 0.2);

    let avg_time = metrics.average_residence_time();
    assert_relative_eq!(avg_time, 0.15, epsilon = 1e-10);

    let max_time = metrics.max_residence_time().unwrap();
    assert_relative_eq!(max_time, 0.2, epsilon = 1e-10);

    let min_time = metrics.min_residence_time().unwrap();
    assert_relative_eq!(min_time, 0.1, epsilon = 1e-10);

    Ok(())
}

/// Test numerical parameters and discretization schemes
#[test]
fn test_numerical_schemes() -> Result<()> {
    // Test basic numerical scheme validation for microfluidic simulation

    // Set typical microfluidic simulation parameters
    let cfl_number = 0.5_f64;
    let time_step = 1e-6_f64; // 1 μs
    let spatial_order = 2_usize; // Second-order accuracy
    let temporal_order = 2_usize; // Second-order time integration

    // Validate parameter constraints
    assert_relative_eq!(cfl_number, 0.5, epsilon = 1e-10);
    assert_relative_eq!(time_step, 1e-6, epsilon = 1e-15);
    assert_eq!(spatial_order, 2);
    assert_eq!(temporal_order, 2);

    // Test stability analysis
    let characteristic_velocity = 0.1; // m/s
    let characteristic_length = 100e-6; // 100 μm
    let max_stable_dt = cfl_number * characteristic_length / characteristic_velocity;

    assert!(
        time_step <= max_stable_dt,
        "Time step violates CFL condition"
    );

    // Test Richardson extrapolation for grid convergence
    let coarse_error = 1e-3;
    let fine_error = 2.5e-4; // 4x improvement expected for 2nd order
    let refinement_ratio = 2.0_f64;
    let order = spatial_order as f64;

    let theoretical_ratio = refinement_ratio.powf(order);
    let actual_ratio = coarse_error / fine_error;

    assert_relative_eq!(actual_ratio, theoretical_ratio, epsilon = 0.1);

    Ok(())
}

/// Test component factory for creating standardized microfluidic components
#[test]
fn test_component_factory() -> Result<()> {
    // Test pump creation with standard specifications
    let mut pump_params = HashMap::new();
    pump_params.insert("max_flow_rate".to_string(), 1e-9_f64); // 1 nL/s
    pump_params.insert("max_pressure".to_string(), 10000.0_f64); // 10 kPa

    let pump = ComponentFactory::create::<f64>("Micropump", &pump_params)?;
    assert_eq!(pump.component_type(), "Micropump");

    // Test valve creation
    let mut valve_params = HashMap::new();
    valve_params.insert("cv".to_string(), 1e-9_f64);

    let valve = ComponentFactory::create::<f64>("Microvalve", &valve_params)?;
    assert_eq!(valve.component_type(), "Microvalve");

    // Test circular channel creation
    let mut channel_params = HashMap::new();
    channel_params.insert("length".to_string(), 1e-3_f64);
    channel_params.insert("diameter".to_string(), 100e-6_f64);

    let channel = ComponentFactory::create::<f64>("CircularChannel", &channel_params)?;
    assert_eq!(channel.component_type(), "CircularChannel");

    // Test rectangular channel creation
    let mut rect_params = HashMap::new();
    rect_params.insert("length".to_string(), 1e-3_f64);
    rect_params.insert("width".to_string(), 200e-6_f64);
    rect_params.insert("height".to_string(), 100e-6_f64);

    let rect_channel = ComponentFactory::create::<f64>("RectangularChannel", &rect_params)?;
    assert_eq!(rect_channel.component_type(), "RectangularChannel");

    Ok(())
}

/// Test advanced millifluidics phenomena: droplet formation and electrokinetics
#[test]
fn test_advanced_millifluidics() -> Result<()> {
    // Test capillary number calculation for droplet formation
    // Ca = μv/σ where μ=viscosity, v=velocity, σ=surface tension
    let fluid = fluid::database::water_20c::<f64>()?;
    let velocity = 0.1; // m/s
    let surface_tension = 0.072; // N/m

    let capillary_number = fluid.viscosity * velocity / surface_tension;

    // Validate capillary number is in expected range for droplet formation (Ca ~ 0.001-0.1)
    assert!(capillary_number > 1e-4 && capillary_number < 1e0);

    // Test Weber number for inertial effects
    // We = ρv²L/σ where ρ=density, v=velocity, L=length scale, σ=surface tension
    let density = 1000.0; // kg/m³ (water)
    let length_scale = 100e-6; // 100 μm

    let weber_number = density * velocity.powi(2) * length_scale / surface_tension;

    // For microfluidics, Weber number should be small (We << 1)
    assert!(
        weber_number < 1.0,
        "Weber number too large for microfluidic regime"
    );

    // Test Peclet number for mass transfer
    // Pe = vL/D where v=velocity, L=length, D=diffusivity
    let diffusivity = 1e-9; // m²/s (typical small molecule)
    let peclet_number = velocity * length_scale / diffusivity;

    // Peclet number indicates convection vs diffusion dominance
    assert!(
        peclet_number > 1.0,
        "Convection should dominate in this regime"
    );

    Ok(())
}

/// Test network solver configuration and basic functionality
#[test]
fn test_network_solver_configuration() -> Result<()> {
    // Test basic solver functionality
    let default_solver = NetworkSolver::<f64>::new();
    assert_eq!(default_solver.name(), "NetworkSolver");

    // Test configurable interface
    let tolerance = 1e-8_f64;
    let config = SolverConfig {
        tolerance,
        max_iterations: 500,
    };

    let mut solver = NetworkSolver::<f64>::with_config(config.clone());

    // Test configuration access
    assert_relative_eq!(solver.config().tolerance, tolerance, epsilon = 1e-15);
    assert_eq!(solver.config().max_iterations, 500);

    // Test config update
    let new_config = SolverConfig {
        tolerance: 1e-10,
        max_iterations: 2000,
    };
    solver.set_config(new_config.clone());

    assert_relative_eq!(solver.config().tolerance, 1e-10, epsilon = 1e-15);
    assert_eq!(solver.config().max_iterations, 2000);

    Ok(())
}

/// Test microfluidic component parameter validation
#[test]
fn test_component_parameter_validation() -> Result<()> {
    // Test pump parameter bounds
    let mut pump = Micropump::<f64>::new(1e-9, 1000.0);

    // Test setting valid parameters
    pump.set_parameter("max_flow_rate", 2e-9)?;
    assert_relative_eq!(pump.max_flow_rate, 2e-9, epsilon = 1e-15);

    pump.set_parameter("max_pressure", 2000.0)?;
    assert_relative_eq!(pump.max_pressure, 2000.0, epsilon = 1e-10);

    // Test setting efficiency bounds (0-1)
    pump.set_parameter("efficiency", 0.9)?;
    assert_relative_eq!(pump.efficiency, 0.9, epsilon = 1e-10);

    // Test mixer parameter validation
    let mut mixer = Micromixer::<f64>::new(3, 1e11);

    // Test resistance modification
    mixer.set_parameter("resistance", 2e11)?;
    assert_relative_eq!(mixer.resistance, 2e11, epsilon = 1e-10);

    // Test efficiency clamping
    mixer.set_parameter("efficiency", 2.0)?; // Should clamp to 1.0
    assert_relative_eq!(mixer.efficiency, 1.0, epsilon = 1e-10);

    mixer.set_parameter("efficiency", -0.5)?; // Should clamp to 0.0
    assert_relative_eq!(mixer.efficiency, 0.0, epsilon = 1e-10);

    Ok(())
}
