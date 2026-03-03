//! Pipe Flow 1D Example - Ergonomic API Design
//!
//! This example demonstrates the 1D CFD API with:
//! - Blueprint-driven topology via cfd-schematics
//! - NetworkBuilderSink bridge to cfd-1d solver
//! - Boundary condition setup and Kirchhoff solve
//! - Derived quantity computation (Re, velocity)

use cfd_1d::network::Network;
use cfd_1d::solver::{NetworkProblem, NetworkSolver};
use cfd_1d::NetworkBuilderSink;
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_schematics::application::use_cases::NetworkGenerationService;
use cfd_schematics::serpentine_chain;
use std::f64::consts::PI;

fn main() -> Result<()> {
    println!("Pipe Flow 1D Example - Ergonomic API Design");
    println!("===========================================");

    // Create fluid using proper constructor
    let water = ConstantPropertyFluid::new(
        "Water (20C)".to_string(),
        998.2,    // density [kg/m^3]
        0.001002, // viscosity [Pa.s]
        4186.0,   // specific heat [J/(kg.K)]
        0.599,    // thermal conductivity [W/(m.K)]
        1482.0,   // speed of sound [m/s]
    );

    println!("Fluid Properties:");
    println!("  Density: {} kg/m^3", water.density);
    println!("  Viscosity: {} Pa.s", water.viscosity);

    // Build network via cfd-schematics blueprint + NetworkBuilderSink bridge
    let pipe_diameter = 0.01; // 1 cm diameter
    let blueprint = serpentine_chain("main_pipe_1d", 1, 1.0, pipe_diameter);

    let sink = NetworkBuilderSink::<f64, _>::new(water.clone());
    let service = NetworkGenerationService::new(sink);
    let mut network: Network<f64, ConstantPropertyFluid<f64>> = service.generate(&blueprint)?;

    println!("\nNetwork created with {} nodes", network.node_count());

    // Find inlet and outlet nodes
    let inlet = network
        .graph
        .node_indices()
        .find(|idx| network.graph[*idx].id == "inlet")
        .expect("inlet node");
    let outlet = network
        .graph
        .node_indices()
        .find(|idx| network.graph[*idx].id == "outlet")
        .expect("outlet node");

    // Set boundary pressures
    network.set_pressure(inlet, 1000.0);
    network.set_pressure(outlet, 0.0);

    // Store fluid properties before moving network
    let fluid_density = water.density;
    let fluid_viscosity = water.viscosity;

    // Solve
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
        println!("  Edge {}: {:.6} m^3/s", edge_idx.index(), flow_rate);
    }

    // Calculate derived quantities for validation
    if let Some((&_edge_idx, &flow_rate)) = flow_rates.iter().next() {
        let pipe_area = PI * (pipe_diameter / 2.0).powi(2);
        let velocity = flow_rate.abs() / pipe_area;
        let reynolds = fluid_density * velocity * pipe_diameter / fluid_viscosity;

        println!("\nDerived Quantities:");
        println!("  Volumetric flow rate: {:.6} m^3/s", flow_rate);
        println!("  Average velocity: {:.4} m/s", velocity);
        println!("  Reynolds number: {:.2}", reynolds);

        if reynolds < 2300.0 {
            println!("  Flow regime: Laminar (Re < 2300)");
        } else {
            println!("  Flow regime: Transitional/Turbulent (Re >= 2300)");
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

    println!("\nExample completed successfully");

    Ok(())
}
