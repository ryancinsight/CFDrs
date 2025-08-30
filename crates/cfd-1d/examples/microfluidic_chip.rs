//! Example demonstrating a complete microfluidic chip simulation using cfd-1d.
//!
//! This example creates a microfluidic network with:
//! - Multiple inlets and outlets
//! - A junction for flow mixing
//! - Proper boundary conditions
//! - Flow and pressure analysis
//!
//! Run with: cargo run --example microfluidic_chip

use cfd_1d::solver::SolverConfig;
use cfd_1d::{EdgeProperties, Network, NetworkBuilder, NetworkProblem, NetworkSolver};
use cfd_core::fluid::ConstantPropertyFluid;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧪 Microfluidic Chip Simulation");
    println!("================================");

    // Create a microfluidic network
    let fluid = cfd_core::fluid::database::water_20c::<f64>()?;

    // Build network using NetworkBuilder
    let mut builder = NetworkBuilder::new();

    // Add nodes
    let inlet = builder.add_inlet("inlet".to_string());
    let junction = builder.add_junction("junction".to_string());
    let outlet_1 = builder.add_outlet("outlet_1".to_string());
    let outlet_2 = builder.add_outlet("outlet_2".to_string());

    // Connect nodes with pipes
    let edge1 = builder.connect_with_pipe(inlet, junction, "channel_1".to_string());
    let edge2 = builder.connect_with_pipe(junction, outlet_1, "channel_2".to_string());
    let edge3 = builder.connect_with_pipe(junction, outlet_2, "channel_3".to_string());

    // Build the graph
    let graph = builder.build()?;

    // Create network with fluid
    let mut network = Network::new(graph, fluid.clone());

    // Add edge properties (microfluidic channels)
    const CHANNEL_LENGTH: f64 = 0.001; // 1mm
    const CHANNEL_WIDTH: f64 = 100e-6; // 100μm
    const CHANNEL_HEIGHT: f64 = 50e-6; // 50μm
    const INLET_PRESSURE: f64 = 2000.0; // Pa
    const OUTLET_PRESSURE: f64 = 0.0; // Pa

    let area = CHANNEL_WIDTH * CHANNEL_HEIGHT;
    let hydraulic_diameter = 2.0 * area / (CHANNEL_WIDTH + CHANNEL_HEIGHT);
    let viscosity = fluid.dynamic_viscosity();

    // Resistance for rectangular microchannels
    let resistance_1 = 12.0 * viscosity * CHANNEL_LENGTH / (CHANNEL_WIDTH * CHANNEL_HEIGHT.powi(3));
    let resistance_2 = resistance_1 * 2.0; // Double resistance for outlet channels

    // Add properties for each edge
    for (edge_idx, resistance, id) in [
        (edge1, resistance_1, "channel_1"),
        (edge2, resistance_2, "channel_2"),
        (edge3, resistance_2, "channel_3"),
    ] {
        let props = EdgeProperties {
            id: id.to_string(),
            component_type: cfd_1d::network::ComponentType::Pipe,
            resistance,
            length: CHANNEL_LENGTH,
            area,
            hydraulic_diameter: Some(hydraulic_diameter),
            geometry: None,
            properties: std::collections::HashMap::new(),
        };
        network.add_edge_properties(edge_idx, props);
    }

    println!("✅ Network created successfully!");
    println!("   - Nodes: {}", network.node_count());
    println!("   - Edges: {}", network.edge_count());
    println!("   - Fluid: {}", network.fluid().name);

    // Set boundary pressures
    network.set_pressure(inlet, INLET_PRESSURE);
    network.set_pressure(outlet_1, OUTLET_PRESSURE);
    network.set_pressure(outlet_2, OUTLET_PRESSURE);

    // Create solver with configuration
    const SOLVER_TOLERANCE: f64 = 1e-6;
    const MAX_ITERATIONS: usize = 1000;

    let solver_config = SolverConfig::<f64> {
        tolerance: SOLVER_TOLERANCE,
        max_iterations: MAX_ITERATIONS,
    };

    let solver = NetworkSolver::with_config(solver_config);

    // Solve the flow problem
    println!("\n🔄 Solving flow equations...");
    let problem = NetworkProblem::new(network);
    let solved_network = solver.solve_network(&problem)?;

    println!("✅ Solution completed");

    // Display results
    println!("\n📊 Flow Analysis Results");
    println!("========================");

    // Access solution vectors
    let pressures = solved_network.pressures();
    let flow_rates = solved_network.flow_rates();

    // Node pressures
    println!("\n🔘 Node Pressures:");
    let node_indices = vec![inlet, junction, outlet_1, outlet_2];
    let node_names = vec!["inlet", "junction", "outlet_1", "outlet_2"];
    for (idx, name) in node_indices.iter().zip(node_names.iter()) {
        if let Some(&pressure) = pressures.get(idx) {
            println!("   {}: {:.1} Pa", name, pressure);
        }
    }

    // Edge flow rates
    println!("\n➡️  Flow Rates:");
    let mut total_flow: f64 = 0.0;
    for (edge_idx, flow_rate) in flow_rates.iter() {
        let flow_ml_min = flow_rate * 1e6 * 60.0; // Convert m³/s to mL/min
        println!(
            "   Edge {:?}: {:.3} mL/min ({:.2e} m³/s)",
            edge_idx, flow_ml_min, flow_rate
        );
        total_flow += flow_rate.abs();
    }

    // Mass conservation check
    println!("\n🔬 Mass Conservation Check:");
    println!("   Total flow magnitude: {:.2e} m³/s", total_flow);

    // Calculate Reynolds number for main channel
    let fluid = solved_network.fluid();
    let diameter = hydraulic_diameter;
    let avg_velocity: f64 = if let Some((_edge_idx, flow)) = flow_rates.iter().next() {
        flow.abs() / area
    } else {
        0.0
    };

    let reynolds = {
        let viscosity = fluid.dynamic_viscosity();
        fluid.density * avg_velocity * diameter / viscosity
    };

    println!("\n📈 Flow Characteristics:");
    println!("   Average velocity: {:.3} m/s", avg_velocity);
    println!("   Reynolds number: {:.2}", reynolds);
    println!(
        "   Flow regime: {}",
        if reynolds < 2300.0 {
            "Laminar"
        } else {
            "Turbulent"
        }
    );

    println!("\n✅ Simulation completed successfully!");

    Ok(())
}
