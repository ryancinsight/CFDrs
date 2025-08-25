//! Pipe flow example demonstrating the 1D CFD solver
//! 
//! This example follows SOLID principles and demonstrates proper error handling

use cfd_suite::prelude::*;
use cfd_suite::core::{Result, BoundaryCondition};
use cfd_suite::d1::{Network, Node, NodeType, ChannelProperties, NetworkSolver, NetworkProblem};

fn main() -> Result<()> {
    println!("Pipe Flow Example");
    println!("========================\n");

    // Create fluid with proper error handling
    let fluid = Fluid::<f64>::water()?;
    println!("Fluid: {}", fluid.name);
    println!("Density: {} kg/m³", fluid.density);
    
    // Build network with fluid
    let mut network = Network::new(fluid.clone());
    
    // Add nodes with proper IDs
    network.add_node(Node::new("inlet".to_string(), NodeType::Inlet));
    network.add_node(Node::new("outlet".to_string(), NodeType::Outlet));
    
    // Add a channel between nodes (1m long, 1mm² cross-section)
    let length = 1.0;
    let area = 1e-6; // 1mm²
    let viscosity = fluid.dynamic_viscosity(1.0)?; // Get viscosity (shear rate doesn't matter for Newtonian)
    let resistance = 8.0 * viscosity * length / (std::f64::consts::PI * area * area);
    let channel_props = ChannelProperties::new(resistance, length, area);
    network.add_edge("inlet", "outlet", channel_props)?;
    
    // Set boundary conditions
    network.set_boundary_condition(
        "inlet", 
        BoundaryCondition::PressureInlet { pressure: 101325.0 }
    )?;
    network.set_boundary_condition(
        "outlet",
        BoundaryCondition::PressureOutlet { pressure: 101225.0 }
    )?;
    
    // Create and configure solver
    let mut solver = NetworkSolver::<f64>::new();
    
    // Create problem and solve
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;
    
    // Display results
    println!("\nResults:");
    println!("--------");
    println!("Network solved with {} nodes", solution.node_count());
    
    // Access pressures from the solution
    let pressures = solution.pressures();
    if pressures.len() > 0 {
        println!("Inlet pressure: {:.2} Pa", pressures[0]);
        if pressures.len() > 1 {
            println!("Outlet pressure: {:.2} Pa", pressures[1]);
            let flow_rate = (pressures[0] - pressures[1]) / resistance;
            println!("Flow rate: {:.6} m³/s", flow_rate);
        }
    }
    
    println!("\nSimulation completed successfully!");
    Ok(())
}