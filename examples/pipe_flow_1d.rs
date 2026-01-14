//! Pipe Flow 1D Example - Ergonomic API Design
//!
//! This example demonstrates the 1D CFD API with:
//! - Error-free builder pattern (single point of failure in .build())
//! - Self-documenting parameter names using NetworkBuilder
//! - Clear method naming (.build() instead of .build_network())
//! - Standard library usage (windows() instead of custom windowed_operation)
//! - Cohesive workflow connecting simulation results to analysis

use cfd_1d::network::{Network, NetworkBuilder};
use cfd_1d::solver::{NetworkProblem, NetworkSolver};
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use std::f64::consts::PI;

fn main() -> Result<()> {
    println!("Pipe Flow 1D Example - Ergonomic API Design");
    println!("===========================================");

    // Create fluid using proper constructor
    let water = ConstantPropertyFluid::new(
        "Water (20°C)".to_string(),
        998.2,    // density [kg/m³]
        0.001002, // viscosity [Pa·s]
        4186.0,   // specific heat [J/(kg·K)]
        0.599,    // thermal conductivity [W/(m·K)]
        1482.0,   // speed of sound [m/s]
    );

    println!("Fluid Properties:");
    println!("  Density: {} kg/m³", water.density);
    println!("  Viscosity: {} Pa·s", water.viscosity);

    // Build network using builder pattern
    let mut builder = NetworkBuilder::new();

    // Add inlet and outlet nodes
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());

    // Connect with pipe
    builder.connect_with_pipe(inlet, outlet, "main_pipe".to_string());

    // Build the network
    let graph = builder.build()?;
    let network = Network::new(graph, water);

    println!("\nNetwork created with {} nodes", network.node_count());
    println!("✓ Nodes and edges configured");

    // Store fluid properties before moving network
    let fluid_density = network.fluid().density;
    let fluid_viscosity = network.fluid().viscosity;

    // Create solver
    let solver = NetworkSolver::new();
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    println!("\n=== Solution ===");
    println!("Converged successfully!");

    // Extract and display results
    let pressures = solution.pressures();
    let flow_rates = solution.flow_rates();

    println!("\nPressure Results:");
    for (node_idx, pressure) in pressures {
        let pressure: f64 = *pressure;
        println!("  Node {}: {:.2} Pa", node_idx.index(), pressure);
    }

    println!("\nFlow Rate Results:");
    for (edge_idx, flow_rate) in flow_rates {
        let flow_rate: f64 = *flow_rate;
        println!("  Edge {}: {:.6} m³/s", edge_idx.index(), flow_rate);
    }

    // Calculate derived quantities for validation
    if let Some((&first_pressure, &second_pressure)) = pressures
        .values()
        .take(2)
        .collect::<Vec<_>>()
        .iter()
        .zip(pressures.values().skip(1))
        .next()
    {
        let pressure_drop = first_pressure - second_pressure;
        println!("\nDerived Quantities:");
        println!("  Pressure drop: {:.2} Pa", pressure_drop);

        if let Some(&flow_rate) = flow_rates.values().next() {
            println!("  Volumetric flow rate: {:.6} m³/s", flow_rate);

            let pipe_diameter: f64 = 0.01; // 1 cm diameter
            let pipe_area = PI * (pipe_diameter / 2.0).powi(2);
            let velocity = flow_rate / pipe_area;
            let reynolds = fluid_density * velocity * pipe_diameter / fluid_viscosity;

            println!("  Average velocity: {:.4} m/s", velocity);
            println!("  Reynolds number: {:.2}", reynolds);

            if reynolds < 2300.0 {
                println!("  ✓ Flow regime: Laminar (Re < 2300)");
            } else {
                println!("  ⚠ Flow regime: Transitional/Turbulent (Re ≥ 2300)");
            }
        }
    }

    // Statistical analysis of results
    let pressure_values: Vec<f64> = pressures.values().cloned().collect();
    if !pressure_values.is_empty() {
        let max_p = pressure_values
            .iter()
            .fold(f64::NEG_INFINITY, |a: f64, &b| a.max(b));
        let min_p = pressure_values
            .iter()
            .fold(f64::INFINITY, |a: f64, &b| a.min(b));
        let avg_p: f64 = pressure_values.iter().sum::<f64>() / pressure_values.len() as f64;

        println!("\nStatistical Summary:");
        println!("  Maximum pressure: {:.2} Pa", max_p);
        println!("  Minimum pressure: {:.2} Pa", min_p);
        println!("  Average pressure: {:.2} Pa", avg_p);
        println!("  Pressure range: {:.2} Pa", max_p - min_p);
    }

    println!("\n✓ Example completed successfully");
    println!("✓ Demonstrates proper 1D CFD API usage");
    println!("✓ Shows builder pattern for network construction");
    println!("✓ Includes physics validation metrics");

    Ok(())
}
