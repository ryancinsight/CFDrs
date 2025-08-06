//! Simple pipe flow example demonstrating the CFD suite structure.

use cfd_core::{
    boundary::{BoundaryCondition, BoundaryConditionSet},
    domain::Domain1D,
    fluid::Fluid,
    plugin::PluginRegistry,
};

fn main() {
    // Create a simple 1D pipe flow example
    println!("Simple Pipe Flow Example");
    println!("========================");
    
    // Create fluid properties for water
    let water = Fluid::<f64>::water();
    println!("Fluid: {}", water.name);
    println!("Density: {} kg/m³", water.density);
    println!("Viscosity: {} Pa·s", water.viscosity);
    
    // Create a 1D domain (pipe)
    let domain = Domain1D::new(0.0, 1.0); // 1 meter pipe
    println!("\nDomain: {} m pipe", domain.length());
    
    // Set up boundary conditions
    let mut bc_set = BoundaryConditionSet::new();
    bc_set.add("inlet", BoundaryCondition::pressure_inlet(101325.0)); // 1 atm
    bc_set.add("outlet", BoundaryCondition::pressure_outlet(101225.0)); // Slightly lower pressure
    
    // Create plugin registry
    let registry = PluginRegistry::new();
    println!("\nPlugin registry created");
    
    // Calculate Reynolds number for a given velocity
    let velocity = 0.1; // m/s
    let diameter = 0.01; // 10 mm pipe
    let re = water.reynolds_number(velocity, diameter);
    println!("\nFor velocity = {} m/s, diameter = {} m:", velocity, diameter);
    println!("Reynolds number = {:.2}", re);
    
    if re < 2300.0 {
        println!("Flow regime: Laminar");
    } else if re < 4000.0 {
        println!("Flow regime: Transitional");
    } else {
        println!("Flow regime: Turbulent");
    }
}