//! Working pipe flow example - pragmatic and functional
//! Following SOLID principles with actual API

use cfd_suite::core::{Fluid, Result};
use cfd_suite::d1::{Network, NetworkBuilder, NetworkSolver, NetworkProblem};
use cfd_suite::d1::Node;

fn main() -> Result<()> {
    println!("Working Pipe Flow Example");
    println!("=========================\n");

    // Create fluid - using the actual API
    let fluid = Fluid::<f64>::water()?;
    println!("Fluid: {}", fluid.name);
    println!("Density: {} kg/m³\n", fluid.density);
    
    // Build network using actual NetworkBuilder API
    let mut builder = NetworkBuilder::new(fluid);
    
    // Add nodes
    let inlet = Node::new(0, 0.0, 0.0, 0.0);
    let outlet = Node::new(1, 1.0, 0.0, 0.0);
    
    builder = builder.add_node(inlet);
    builder = builder.add_node(outlet);
    
    let mut network = builder.build()?;
    
    // Set boundary conditions
    network.set_pressure(0, 101325.0)?;  // 1 atm inlet
    network.set_pressure(1, 101225.0)?;  // Lower outlet
    
    // Create solver with actual config structure
    let mut solver = NetworkSolver::<f64>::new();
    
    // Create problem
    let problem = NetworkProblem { network };
    
    // Solve
    let solution = solver.solve(&problem)?;
    
    // Display results
    println!("Solution:");
    println!("---------");
    println!("Nodes in solution: {}", solution.node_count());
    
    for i in 0..solution.node_count() {
        if let Some(node) = solution.get_node(i) {
            println!("Node {}: Pressure = {:.2} Pa", i, node.pressure);
        }
    }
    
    println!("\n✅ Simulation completed successfully!");
    Ok(())
}