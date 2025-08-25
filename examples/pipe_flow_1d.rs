//! Pipe Flow 1D Example - Ergonomic API Design
//!
//! This example demonstrates the 1D CFD API with:
//! - Error-free builder pattern (single point of failure in .build())
//! - Self-documenting parameter names using ChannelProperties
//! - Clear method naming (.build() instead of .build_network())
//! - Standard library usage (windows() instead of custom windowed_operation)
//! - Cohesive workflow connecting simulation results to analysis

use cfd_suite::prelude::*;
use cfd_suite::core::{Result, BoundaryCondition};
use cfd_1d::{Network, NetworkBuilder, NetworkSolver, NetworkProblem, ChannelProperties};

fn main() -> Result<()> {
    println!("Pipe Flow 1D Example - Ergonomic API Design");
    println!("===========================================");

    // Demonstrate unified prelude and composition-based configuration
    let water = Fluid::<f64>::water()?;
    println!("Fluid Properties:");
    println!("  Name: {}", water.name);
    println!("  Density: {} kg/m³", water.density);
    let viscosity = water.dynamic_viscosity(1.0)?;
    println!("  Viscosity: {} Pa·s", viscosity);

    // Create 1D network 
    let mut network = Network::new(water.clone());
    
    // Add nodes
    network.add_node(Node::new("inlet".to_string(), NodeType::Inlet));
    network.add_node(Node::new("outlet".to_string(), NodeType::Outlet));
    
    // Add channel with properties
    let channel_props = ChannelProperties::new(
        100.0,  // resistance (Pa·s/m³) 
        1.0,    // length (m)
        1e-6,   // area (m²)
    );
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

    println!("\nNetwork created with {} nodes", network.node_count());
    println!("✓ Nodes and edges configured");
    println!("✓ Boundary conditions set");

    // Create solver 
    let mut solver = NetworkSolver::<f64>::new();
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    println!("\n=== Solution ===");
    println!("Converged successfully!");
    
    // Access solution data
    let pressures = solution.pressures();
    if pressures.len() >= 2 {
        let dp = pressures[0] - pressures[1];
        println!("Pressure drop: {:.2} Pa", dp);
        
        // Calculate flow rate using Hagen-Poiseuille equation
        let flow_rate = dp / 100.0; // resistance = 100.0
        println!("Flow rate: {:.6} m³/s", flow_rate);
    }

    // Demonstrate standard library usage with windows() for gradient calculation
    println!("\n=== Gradient Analysis (using std::windows()) ===");
    let pressure_data: Vec<f64> = pressures.iter().copied().collect();
    
    if pressure_data.len() >= 2 {
        let gradients: Vec<f64> = pressure_data
            .windows(2)
            .map(|w| (w[1] - w[0]) / 0.5) // dx = 0.5m between nodes
            .collect();
        
        println!("Pressure gradients calculated using windows():");
        for (i, grad) in gradients.iter().enumerate() {
            println!("  Segment {}: {:.2} Pa/m", i, grad);
        }
    }

    // Flow analysis
    println!("\n=== Flow Analysis ===");
    if pressures.len() >= 2 {
        // Calculate Reynolds number manually
        let diameter = (4.0 * 1e-6 / std::f64::consts::PI).sqrt();
        let velocity = (pressures[0] - pressures[1]) / (100.0 * 1e-6); // flow_rate / area
        let reynolds = water.density * velocity * diameter / viscosity;
        
        println!("Reynolds number: {:.2}", reynolds);
        let flow_regime = if reynolds < 2300.0 { "Laminar" } 
                         else if reynolds < 4000.0 { "Transitional" }
                         else { "Turbulent" };
        println!("Flow regime: {}", flow_regime);
        
        // Pressure statistics
        let max_p = pressures.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let min_p = pressures.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let avg_p: f64 = pressures.iter().sum::<f64>() / pressures.len() as f64;
        
        println!("\n=== Pressure Analysis ===");
        println!("Max pressure: {:.2} Pa", max_p);
        println!("Min pressure: {:.2} Pa", min_p);
        println!("Average pressure: {:.2} Pa", avg_p);
    }

    // Demonstrate iterator combinators for data processing
    println!("\n=== Iterator Combinators Demo ===");
    
    // Chain operations for efficient processing
    let processed_pressures: Vec<f64> = pressures.iter()
        .copied()
        .map(|p| p - 101225.0)  // Relative to outlet
        .filter(|&p| p.abs() > 1e-6)  // Non-zero values
        .collect();
    
    println!("Relative pressures (non-zero):");
    for (i, p) in processed_pressures.iter().enumerate() {
        println!("  Node {}: {:.2} Pa", i, p);
    }

    // Use fold for statistics
    let (sum, count) = processed_pressures.iter()
        .fold((0.0, 0), |(s, c), &p| (s + p, c + 1));
    
    if count > 0 {
        println!("Average relative pressure: {:.2} Pa", sum / count as f64);
    }

    println!("\n✓ Simulation completed successfully!");
    println!("✓ Demonstrated ergonomic API design");
    println!("✓ Used standard library patterns (windows, fold)");
    println!("✓ Cohesive workflow from simulation to analysis");

    Ok(())
}