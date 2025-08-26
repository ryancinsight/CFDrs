//! Example demonstrating a complete microfluidic chip simulation using cfd-1d.
//!
//! This example creates a microfluidic network with:
//! - Multiple inlets and outlets
//! - A junction for flow mixing
//! - Proper boundary conditions
//! - Flow and pressure analysis
//! Run with: cargo run --example microfluidic_chip

use cfd_1d::solver::SolverConfig;
use cfd_1d::{ChannelProperties, NetworkBuilder, NetworkProblem, NetworkSolver, Node, NodeType};
use cfd_core::fluid::Fluid;
use cfd_core::BoundaryCondition;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧪 Microfluidic Chip Simulation");
    println!("================================");
    // Create a microfluidic network
    let fluid = Fluid::<f64>::water()?;
    let network = NetworkBuilder::new(fluid)
        // Standard T-junction network
        .add_node(Node::new("inlet".to_string(), NodeType::Inlet))
        .add_node(Node::new("junction".to_string(), NodeType::Junction))
        .add_node(Node::new("outlet_1".to_string(), NodeType::Outlet))
        .add_node(Node::new("outlet_2".to_string(), NodeType::Outlet))
        // Channels connecting the network
        .add_edge(
            "inlet",
            "junction",
            ChannelProperties::new(100.0, 0.001, 100e-6),
        )?
            "outlet_1",
            ChannelProperties::new(200.0, 0.001, 100e-6),
            "outlet_2",
        .build();
    println!("✅ Network created successfully!");
    println!("   - Nodes: {}", network.node_count());
    println!("   - Edges: {}", network.edge_count());
    println!("   - Fluid: {}", network.fluid().name);
    // Set boundary conditions
    let mut network = network;
    network.set_boundary_condition(
        "inlet",
        BoundaryCondition::PressureInlet { pressure: 2000.0 },
    )?;
        "outlet_1",
        BoundaryCondition::PressureOutlet { pressure: 0.0 },
        "outlet_2",
    // Create solver with configuration
    let solver_config = SolverConfig::<f64> {
        tolerance: 1e-6,
        max_iterations: 1000,
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
    for (i, node) in solved_network.nodes().enumerate() {
        if i < pressures.len() {
            println!("   {}: {:.1} Pa", node.id, pressures[i]);
        }
    }
    // Edge flow rates
    println!("\n➡️  Flow Rates:");
    let mut total_flow: f64 = 0.0;
    for (i, _edge) in solved_network.edges().enumerate() {
        if i < flow_rates.len() {
            let flow_rate: f64 = flow_rates[i];
            let flow_ml_min = flow_rate * 1e6 * 60.0; // Convert m³/s to mL/min
            println!(
                "   Edge {}: {:.3} mL/min ({:.2e} m³/s)",
                i, flow_ml_min, flow_rate
            );
            total_flow += flow_rate.abs();
    // Mass conservation check
    println!("\n🔬 Mass Conservation Check:");
    println!("   Total flow magnitude: {:.2e} m³/s", total_flow);
    // Calculate Reynolds number for main channel
    let fluid = solved_network.fluid();
    let diameter: f64 = 100e-6; // 100 μm channel
    let avg_velocity: f64 = if !flow_rates.is_empty() {
        let flow: f64 = flow_rates[0];
        flow.abs() / (std::f64::consts::PI * (diameter / 2.0).powi(2))
    } else {
        0.0
    let reynolds = if let Ok(viscosity) = fluid.dynamic_viscosity(20.0) {
        fluid.density * avg_velocity * diameter / viscosity
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
}
}
}
}
