//! Simple pipe flow example demonstrating the 1D CFD solver
//! 
//! This example follows SOLID principles and demonstrates proper error handling

use cfd_suite::prelude::*;
use cfd_suite::core::Result;

fn main() -> Result<()> {
    println!("Simple Pipe Flow Example");
    println!("========================\n");

    // Create fluid with proper error handling
    let fluid = Fluid::<f64>::water()?;
    println!("Fluid: {}", fluid.name);
    println!("Density: {} kg/m³", fluid.density);
    
    // Build network with fluid
    let mut network = NetworkBuilder::new(fluid)
        .add_node(Node::new(0, 0.0, 0.0, 0.0))  // Inlet
        .add_node(Node::new(1, 1.0, 0.0, 0.0))  // Outlet
        .build()?;
    
    // Set boundary conditions
    network.set_pressure(0, 101325.0)?;  // 1 atm at inlet
    network.set_pressure(1, 101225.0)?;  // Slightly lower at outlet
    
    // Create and configure solver
    let config = SolverConfig {
        tolerance: 1e-6,
        max_iterations: 1000,
    };
    
    let mut solver = NetworkSolver::with_config(config);
    
    // Create problem and solve
    let problem = NetworkProblem { network };
    let solution = solver.solve(&problem)?;
    
    // Display results
    println!("\nResults:");
    println!("--------");
    for i in 0..solution.node_count() {
        if let Some(node) = solution.get_node(i) {
            println!("Node {}: Pressure = {:.2} Pa", i, node.pressure);
        }
    }
    
    println!("\nSimulation completed successfully!");
    Ok(())
}