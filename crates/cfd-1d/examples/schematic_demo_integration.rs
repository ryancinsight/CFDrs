//! Example demonstrating cfd-1d simulation using cfd-schematics specifications.
//!
//! This example shows how to:
//! 1. Define a network using `cfd-schematics` types (`NodeSpec`, `ChannelSpec`).
//! 2. Convert these specifications into `cfd-1d` simulation objects (`Node`, `Edge`).
//! 3. Build and solve the network.
//!
//! Run with: cargo run -p cfd-1d --example schematic_demo_integration

use cfd_1d::network::{Edge, Network, NetworkBuilder, Node};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::fluid::database;
use cfd_schematics::domain::model::{ChannelSpec, EdgeId, NodeId, NodeKind, NodeSpec};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔌 Schematic Integration Demo");
    println!("============================");

    // 1. Define Schematic Specifications (Design Phase)
    // In a real app, these might come from a JSON file or UI
    let inlet_spec = NodeSpec::new("inlet", NodeKind::Inlet);
    let junction_spec = NodeSpec::new("junction", NodeKind::Junction);
    let outlet_spec = NodeSpec::new("outlet", NodeKind::Outlet);

    // Define a pipe connecting Inlet -> Junction
    // Length: 10mm, Diameter: 1mm
    // Resistance calculated roughly for laminar flow (Poiseuille) or set explicitly
    // For this demo, we'll set resistance to 1.0 (arbitrary units or pre-calculated)
    let pipe_spec = ChannelSpec::new_pipe(
        "pipe1",
        "inlet",
        "junction",
        0.01,   // 10mm length
        0.001,  // 1mm diameter
        1.0e8,  // Resistance (Pa*s/m^3)
        0.0,    // Quadratic coeff
    );

    // Define a valve connecting Junction -> Outlet
    // Cv = 0.5 (flow coefficient)
    let valve_spec = ChannelSpec::new_valve(
        "valve1",
        "junction",
        "outlet",
        0.5,
    );

    println!("✅ Defined specifications:");
    println!("   - Nodes: {:?}, {:?}, {:?}", inlet_spec.id, junction_spec.id, outlet_spec.id);
    println!("   - Edges: {:?}, {:?}", pipe_spec.id, valve_spec.id);

    // 2. Build Simulation Network (Simulation Phase)
    let mut builder = NetworkBuilder::<f64>::new();

    // Convert specs to nodes and add them
    // Note: We need to keep track of NodeIndices to connect edges
    let inlet_idx = builder.add_node(Node::from(&inlet_spec));
    let junction_idx = builder.add_node(Node::from(&junction_spec));
    let outlet_idx = builder.add_node(Node::from(&outlet_spec));

    // Convert specs to edges and connect
    // Edge::from(&ChannelSpec) creates the Edge object
    // We connect them using the indices we got earlier
    builder.add_edge(inlet_idx, junction_idx, Edge::from(&pipe_spec));
    builder.add_edge(junction_idx, outlet_idx, Edge::from(&valve_spec));

    let graph = builder.build()?;
    println!("\n✅ Built simulation graph");

    // 3. Setup Physics and Solve
    let fluid = database::water_20c::<f64>()?;
    let mut network = Network::new(graph, fluid);

    // Set boundary conditions
    // Inlet pressure: 1000 Pa
    // Outlet pressure: 0 Pa
    network.set_pressure(inlet_idx, 1000.0);
    network.set_pressure(outlet_idx, 0.0);

    // Solver configuration
    let config = SolverConfig {
        tolerance: 1e-6,
        max_iterations: 100,
    };
    let solver = NetworkSolver::with_config(config);

    println!("\n🔄 Solving...");
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    println!("✅ Solution converged!");

    // 4. Analyze Results
    println!("\n📊 Results:");
    
    // Check flow through the valve
    // We need to find the edge index for "valve1"
    if let Some(edge_ref) = solution.graph.edge_references().find(|e| e.weight().id == "valve1") {
       let flow = solution.flow_rates.get(&edge_ref.id()).unwrap();
       println!("   Valve Flow Rate: {:.4e} m^3/s", flow);
       println!("   Valve Flow Rate: {:.2} mL/min", flow * 1e6 * 60.0);
    }

    // Check pressure at junction
    if let Some(node_idx) = solution.graph.node_indices().find(|&i| solution.graph[i].id == "junction") {
        let pressure = solution.pressures.get(&node_idx).unwrap();
        println!("   Junction Pressure: {:.2} Pa", pressure);
    }

    Ok(())
}

// Helper to access edge references
use petgraph::visit::EdgeRef;
