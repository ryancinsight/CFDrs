//! Pipe flow 1D example demonstrating unified SSOT prelude and design principles.
//!
//! This example showcases:
//! - Unified prelude usage (SSOT principle)
//! - Composition-based configuration (SOLID principles)
//! - Iterator-based operations (CUPID principles)
//! - Zero-copy abstractions and iterator patterns

use cfd_suite::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    println!("Pipe Flow 1D Example - Unified SSOT Design");
    println!("================================================");

    // Demonstrate unified prelude and composition-based configuration
    let water = Fluid::<f64>::water();
    println!("Fluid: {}", water.name);
    println!("Density: {} kg/m³", water.density);
    println!("Viscosity: {} Pa·s", water.viscosity);

    // Create 1D network using unified builder pattern
    let mut network = NetworkBuilder::<f64>::new()
        .add_inlet_pressure("inlet", 0.0, 0.0, 101325.0)?    // ID, x, y, pressure (Pa)
        .add_outlet_pressure("outlet", 1.0, 0.0, 101225.0)?   // ID, x, y, pressure (Pa)
        .add_channel("ch1", "inlet", "outlet", 100.0, 1.0, 1e-6)? // ID, from, to, resistance, length, area
        .build()?;

    println!("\nNetwork created with unified builder pattern");

    // Demonstrate configuration system using composition
    let config = SolverConfig::<f64>::builder()
        .tolerance(1e-8)
        .max_iterations(1000)
        .relaxation_factor(0.9)
        .parallel(true)
        .verbosity(1)
        .build_network(); // Build network-specific config

    println!("Configuration built using composition pattern");

    // Create solver with configuration
    let solver = NetworkSolver::with_config(config);
    let solution = solver.solve_steady_state(&mut network)?;

    println!("\nSolution Results:");
    println!("Converged: {}", solution.converged);
    println!("Iterations: {}", solution.iterations);
    println!("Residual: {:.2e}", solution.residual);

    // Demonstrate iterator patterns for analysis
    let velocities = vec![0.05, 0.1, 0.15, 0.2, 0.25, 0.3];
    let diameter = 0.01; // 10 mm pipe

    println!("\nReynolds Number Analysis using Iterator Combinators:");
    println!("Velocity (m/s) | Reynolds | Flow Regime");
    println!("---------------|----------|------------");

    // Use iterator combinators for functional analysis
    let reynolds_analysis: Vec<_> = velocities.iter()
        .map(|&v| {
            let re = water.reynolds_number(v, diameter);
            let regime = match re {
                r if r < 2300.0 => "Laminar",
                r if r < 4000.0 => "Transitional",
                _ => "Turbulent",
            };
            (v, re, regime)
        })
        .collect();

    // Display results using iterator patterns
    reynolds_analysis.iter()
        .for_each(|(v, re, regime)| {
            println!("{:13.2} | {:8.0} | {}", v, re, regime);
        });

    // Demonstrate zero-copy slice operations
    let re_values: Vec<f64> = reynolds_analysis.iter()
        .map(|(_, re, _)| *re)
        .collect();

    // Use mathematical operations
    use cfd_math::SliceOps;
    let mean_re = re_values.iter().cloned().mean().unwrap_or(0.0);
    let max_re = re_values.iter().fold(0.0_f64, |acc, &x| acc.max(x));

    println!("\nStatistical Analysis:");
    println!("Mean Reynolds number: {:.1}", mean_re);
    println!("Maximum Reynolds number: {:.1}", max_re);

    // Demonstrate windowed operations for gradient analysis
    let velocity_gradients: Vec<_> = SliceOps::windowed_operation(
        &velocities,
        2,
        |window| window[1] - window[0]
    );

    println!("\nVelocity gradients: {:?}", velocity_gradients);

    Ok(())
}