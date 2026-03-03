//! Pipe flow example demonstrating the 1D CFD solver
//!
//! This example builds a single-pipe network via cfd-schematics blueprint
//! generation and solves for pressure / flow rate using the Kirchhoff solver.

use cfd_1d::network::Network;
use cfd_1d::solver::{NetworkProblem, NetworkSolver};
use cfd_1d::NetworkBuilderSink;
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};
use cfd_schematics::application::use_cases::NetworkGenerationService;
use cfd_schematics::serpentine_chain;
use petgraph::visit::EdgeRef;

fn main() -> Result<()> {
    println!("Pipe Flow Example");
    println!("========================\n");

    // Create fluid with proper error handling
    let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
    println!("Fluid: {}", fluid.name());
    println!("Density: {} kg/m^3", fluid.density);
    println!("Viscosity: {} Pa.s", fluid.viscosity);

    // Build network via cfd-schematics blueprint -> NetworkBuilderSink bridge
    const PIPE_LENGTH: f64 = 1.0; // meters
    const PIPE_DIAMETER: f64 = 0.001; // 1 mm

    // serpentine_chain with 1 segment produces a simple pipe
    let blueprint = serpentine_chain("single_pipe", 1, PIPE_LENGTH, PIPE_DIAMETER);

    let sink = NetworkBuilderSink::<f64, _>::new(fluid.clone());
    let service = NetworkGenerationService::new(sink);
    let mut network: Network<f64, ConstantPropertyFluid<f64>> = service.generate(&blueprint)?;

    let inlet = network
        .graph
        .node_indices()
        .find(|idx| network.graph[*idx].id == "inlet")
        .expect("inlet node exists in generated blueprint");
    let outlet = network
        .graph
        .node_indices()
        .find(|idx| network.graph[*idx].id == "outlet")
        .expect("outlet node exists in generated blueprint");

    const INLET_PRESSURE: f64 = 101325.0; // Pa (1 atm)
    const OUTLET_PRESSURE: f64 = 101225.0; // Pa

    // Set boundary conditions (inlet and outlet pressures)
    network.set_pressure(inlet, INLET_PRESSURE);
    network.set_pressure(outlet, OUTLET_PRESSURE);

    // Create and configure solver
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::new();

    // Create problem and solve
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    // Display results
    println!("\nResults:");
    println!("--------");
    println!("Network solved with {} nodes", solution.node_count());

    // Report flow rates on each edge
    for edge_ref in solution.graph.edge_references() {
        let edge = edge_ref.weight();
        if let Some(&flow) = solution.flow_rates.get(&edge_ref.id()) {
            println!("Edge '{}': flow = {:.6e} m^3/s", edge.id, flow);
        }
    }

    // Calculate Reynolds number using first edge flow rate
    if let Some((&_edge_idx, &flow_rate)) = solution.flow_rates.iter().next() {
        let pipe_area = std::f64::consts::PI * (PIPE_DIAMETER / 2.0).powi(2);
        let velocity = flow_rate.abs() / pipe_area;
        let reynolds = fluid.density * velocity * PIPE_DIAMETER / fluid.viscosity;
        println!("\nDerived:");
        println!("  Flow rate: {:.6e} m^3/s", flow_rate);
        println!("  Velocity: {:.4} m/s", velocity);
        println!("  Reynolds number: {:.2}", reynolds);
    }

    println!("\nSimulation completed successfully!");
    Ok(())
}
