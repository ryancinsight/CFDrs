//! Simple pipe flow example demonstrating the CFD suite structure.

use cfd_core::prelude::*;

fn main() {
    println!("CFD Simulation Suite - Simple Pipe Flow Example");
    println!("==============================================");
    
    // Create a simple fluid
    let water = Fluid::<f64>::water();
    println!("Fluid: {}", water.name);
    println!("  Density: {} kg/m³", water.density);
    println!("  Viscosity: {} Pa·s", water.viscosity);
    
    // Create a 1D domain
    let domain = cfd_core::domain::Domain1D::new(0.0, 1.0);
    println!("\nDomain: 1D pipe");
    println!("  Length: {} m", domain.length());
    
    // Define boundary conditions
    let mut bc_set = BoundaryConditionSet::new();
    bc_set.add("inlet", BoundaryCondition::pressure_inlet(101325.0));
    bc_set.add("outlet", BoundaryCondition::pressure_outlet(100000.0));
    
    println!("\nBoundary Conditions:");
    for (name, _bc) in &bc_set.conditions {
        println!("  {}", name);
    }
    
    // Plugin registry demonstration
    let registry = PluginRegistry::new();
    println!("\nPlugin System: Ready");
    println!("  Available plugins: {:?}", registry.list().unwrap_or_default());
    
    println!("\nSimulation setup complete!");
}