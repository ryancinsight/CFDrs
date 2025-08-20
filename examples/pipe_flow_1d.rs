//! Pipe Flow 1D Example - Ergonomic API Design
//!
//! This example demonstrates the 1D CFD API with:
//! - Error-free builder pattern (single point of failure in .build())
//! - Self-documenting parameter names using ChannelProperties
//! - Clear method naming (.build() instead of .build_network())
//! - Standard library usage (windows() instead of custom windowed_operation)
//! - Cohesive workflow connecting simulation results to analysis

use cfd_suite::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("Pipe Flow 1D Example - Ergonomic API Design");
    println!("===========================================");

    // Demonstrate unified prelude and composition-based configuration
    let water = Fluid::<f64>::water();
    println!("Fluid Properties:");
    println!("  Name: {}", water.name);
    println!("  Density: {} kg/m³", water.density);
    println!("  Viscosity: {} Pa·s", water.viscosity);

    // Create 1D network using improved builder pattern (Issue 1: Error-free configuration)
    let mut network = NetworkBuilder::<f64>::new()
        .add_inlet_pressure("inlet", 0.0, 0.0, 101325.0)    // No Result<> needed
        .add_outlet_pressure("outlet", 1.0, 0.0, 101225.0)   // Fluent API
        .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(
            100.0,  // resistance (Pa·s/m³) - self-documenting!
            1.0,    // length (m)
            1e-6,   // area (m²)
        ))
        .build()?;  // Single point of failure

    println!("\nNetwork created with improved builder pattern");
    println!("✓ Error-free configuration methods");
    println!("✓ Self-documenting ChannelProperties");

    // Demonstrate improved configuration system (Issue 3: Clear method naming)
    let config = SolverConfig::<f64>::builder()
        .tolerance(1e-8)
        .max_iterations(1000)
        .relaxation_factor(0.9)
        .parallel(true)
        .verbosity(1)
        .build(); // Clear, unambiguous method name

    println!("✓ Configuration built using clear .build() method");

    // Create solver and solve the network
    let solver = NetworkSolver::with_config(config);
    let solution = solver.solve_steady_state(&mut network)?;

    println!("\nSolution Results:");
    println!("  Converged: {}", solution.converged);
    println!("  Iterations: {}", solution.iterations);
    println!("  Residual: {:.2e}", solution.residual);

    // Issue 5: Connect simulation results to analysis (cohesive workflow)
    println!("\nExtracting Flow Results for Analysis:");
    
    // Get the calculated flow rate from the solved network
    let channel = network.get_edge("ch1").unwrap();
    let flow_rate = channel.flow_rate.unwrap_or(0.0);  // m³/s
    
    // Calculate velocity from flow rate and channel area
    let channel_area = 1e-6;  // m² (from ChannelProperties)
    let calculated_velocity = flow_rate / channel_area;  // m/s
    
    println!("  Flow rate: {:.6e} m³/s", flow_rate);
    println!("  Calculated velocity: {:.4} m/s", calculated_velocity);

    // Create a velocity range around the calculated value for comprehensive analysis
    let base_velocity = calculated_velocity.max(0.01); // Ensure minimum velocity
    let velocity_range: Vec<f64> = (0..6)
        .map(|i| base_velocity * (0.5 + i as f64 * 0.2))  // Create range: 0.5x to 1.5x
        .collect();

    let diameter = (4.0 * channel_area / std::f64::consts::PI).sqrt(); // Hydraulic diameter

    println!("\nReynolds Number Analysis for Flow Range:");
    println!("Velocity (m/s) | Reynolds | Flow Regime     | Status");
    println!("---------------|----------|-----------------|--------");

    // Use iterator combinators for functional analysis
    let reynolds_analysis: Vec<_> = velocity_range.iter()
        .map(|&v| {
            let re = water.reynolds_number(v, diameter);
            let regime = match re {
                r if r < 2300.0 => "Laminar",
                r if r < 4000.0 => "Transitional",
                _ => "Turbulent",
            };
            let status = if (v - calculated_velocity).abs() < 1e-6 { 
                "← Calculated" 
            } else { 
                "" 
            };
            (v, re, regime, status)
        })
        .collect();

    // Display results using iterator patterns
    reynolds_analysis.iter()
        .for_each(|(v, re, regime, status)| {
            println!("{:13.4} | {:8.0} | {:15} | {}", v, re, regime, status);
        });

    // Issue 4: Use standard library instead of custom implementation
    println!("\nVelocity Gradient Analysis (using std::slice::windows):");
    
    // OLD (custom implementation):
    // let velocity_gradients: Vec<_> = SliceOps::windowed_operation(
    //     &velocities, 2, |window| window[1] - window[0]
    // );
    
    // NEW (standard library - idiomatic Rust):
    let velocity_gradients: Vec<_> = velocity_range
        .windows(2)
        .map(|window| window[1] - window[0])
        .collect();

    println!("✓ Using standard library windows() method");
    velocity_gradients.iter().enumerate().for_each(|(i, &grad)| {
        println!("  Gradient {}: {:.4} m/s per step", i + 1, grad);
    });

    // Demonstrate zero-copy statistical operations using iterators
    let re_values: Vec<f64> = reynolds_analysis.iter()
        .map(|(_, re, _, _)| *re)
        .collect();

    use cfd_math::MathIteratorExt;
    let mean_re = re_values.iter().cloned().mean().unwrap_or(0.0);
    let variance_re = re_values.iter().cloned().variance().unwrap_or(0.0);
    let std_dev_re = variance_re.sqrt();

    println!("\nStatistical Analysis:");
    println!("  Mean Reynolds number: {:.1}", mean_re);
    println!("  Standard deviation: {:.1}", std_dev_re);
    println!("  Flow regime classification: {}", 
        if mean_re < 2300.0 { "Primarily Laminar" }
        else if mean_re < 4000.0 { "Transitional Range" }
        else { "Primarily Turbulent" }
    );

    // Summary of improvements demonstrated
    println!("\n🎉 API Improvements Demonstrated:");
    println!("✓ Issue 1: Error-free builder pattern with single .build() failure point");
    println!("✓ Issue 2: Self-documenting ChannelProperties instead of magic numbers");
    println!("✓ Issue 3: Clear .build() method name instead of .build_network()");
    println!("✓ Issue 4: Standard library windows() instead of custom windowed_operation");
    println!("✓ Issue 5: Cohesive workflow connecting simulation to analysis");

    Ok(())
}