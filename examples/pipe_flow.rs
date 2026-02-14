//! Pipe flow example demonstrating the 1D CFD solver
//!
//! This example follows SOLID principles and demonstrates proper error handling

use cfd_1d::network::{ComponentType, EdgeProperties};
use cfd_1d::solver::{NetworkProblem, NetworkSolver};
use cfd_1d::Network;
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};
use cfd_fluidics::{serpentine_chain, FluidicDesigner};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

fn main() -> Result<()> {
    println!("Pipe Flow Example");
    println!("========================\n");

    // Create fluid with proper error handling
    let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
    println!("Fluid: {}", fluid.name());
    println!("Density: {} kg/m³", fluid.density);
    println!("Viscosity: {} Pa·s", fluid.viscosity);

    // Build network topology via cfd-fluidics (phase-2 architecture split)
    const PIPE_LENGTH: f64 = 1.0; // meters
    const PIPE_AREA: f64 = 1e-6; // m² (1mm²)
    let hydraulic_diameter = 2.0 * (PIPE_AREA / std::f64::consts::PI).sqrt();

    let designer = FluidicDesigner::new();
    let blueprint = serpentine_chain("single_pipe", 1, PIPE_LENGTH, hydraulic_diameter);
    let graph = designer.generate(&blueprint)?;

    let inlet = graph
        .node_indices()
        .find(|idx| graph[*idx].id == "inlet")
        .expect("inlet node exists in generated blueprint");
    let outlet = graph
        .node_indices()
        .find(|idx| graph[*idx].id == "outlet")
        .expect("outlet node exists in generated blueprint");
    let edge_idx = graph
        .edge_references()
        .find(|e| e.weight().id == "segment_1")
        .map(|e| e.id())
        .expect("segment_1 edge exists in generated blueprint");

    // Create network with fluid
    let mut network = Network::new(graph, fluid);

    // Add edge properties (1m long, 1mm² cross-section)
    const HAGEN_POISEUILLE_FACTOR: f64 = 8.0;
    const INLET_PRESSURE: f64 = 101325.0; // Pa (1 atm)
    const OUTLET_PRESSURE: f64 = 101225.0; // Pa

    let viscosity = network.fluid().viscosity;
    let resistance = HAGEN_POISEUILLE_FACTOR * viscosity * PIPE_LENGTH
        / (std::f64::consts::PI * PIPE_AREA * PIPE_AREA);

    let edge_props = EdgeProperties {
        id: "pipe1".to_string(),
        component_type: ComponentType::Pipe,
        resistance,
        length: PIPE_LENGTH,
        area: PIPE_AREA,
        hydraulic_diameter: Some(hydraulic_diameter),
        geometry: None,
        properties: HashMap::new(),
    };
    network.add_edge_properties(edge_idx, edge_props);

    // Set boundary conditions (inlet and outlet pressures)
    network.set_pressure(inlet, INLET_PRESSURE);
    network.set_pressure(outlet, OUTLET_PRESSURE);

    // Create and configure solver
    let solver = NetworkSolver::<f64>::new();

    // Create problem and solve
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    // Display results
    println!("\nResults:");
    println!("--------");
    println!("Network solved with {} nodes", solution.node_count());

    // Calculate flow rate using pressure difference
    let pressure_drop = INLET_PRESSURE - OUTLET_PRESSURE;
    let flow_rate = pressure_drop / resistance;
    println!("Flow rate: {:.6} m³/s", flow_rate);

    // Reynolds number calculation
    let diameter = 2.0 * (PIPE_AREA / std::f64::consts::PI).sqrt();
    let velocity = flow_rate / PIPE_AREA;
    let reynolds = solution.fluid().density * velocity * diameter / viscosity;
    println!("Reynolds number: {:.2}", reynolds);

    println!("\nSimulation completed successfully!");
    Ok(())
}
