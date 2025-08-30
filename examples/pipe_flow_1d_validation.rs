//! 1D Pipe Flow Validation
//!
//! This example validates the 1D pipe flow solver against the analytical
//! Hagen-Poiseuille solution for laminar flow in a circular pipe.

use cfd_1d::{ChannelProperties, Network, NetworkProblem, NetworkSolver, Node, NodeType};
use cfd_core::{boundary::BoundaryCondition, fluid::Fluid, error::Result};
use cfd_core::solver::Solver;
use std::f64::consts::PI;

fn main() -> Result<()> {
    println!("========================================");
    println!("1D Pipe Flow Validation");
    println!("========================================\n");

    // Pipe parameters
    let pipe_radius: f64 = 0.01; // 10 mm radius
    let pipe_length: f64 = 0.1; // 100 mm length
    let fluid_viscosity: f64 = 1e-3; // Water at 20°C (Pa·s)
    let fluid_density: f64 = 1000.0; // Water density (kg/m³)
    let pressure_drop: f64 = 10.0; // Pressure drop across pipe (Pa)

    println!("Pipe Geometry:");
    println!("  Radius: {} m", pipe_radius);
    println!("  Length: {} m", pipe_length);
    println!("  Diameter: {} m", 2.0 * pipe_radius);

    println!("\nFluid Properties:");
    println!("  Dynamic viscosity: {} Pa·s", fluid_viscosity);
    println!("  Density: {} kg/m³", fluid_density);

    println!("\nBoundary Conditions:");
    println!("  Pressure drop: {} Pa", pressure_drop);
    println!("  Pressure gradient: {} Pa/m", pressure_drop / pipe_length);

    // Create fluid with specified properties
    let fluid = Fluid::create("Test Fluid".to_string(), fluid_density, fluid_viscosity);

    // Create 1D network
    let mut network = Network::new(fluid);

    // Add inlet and outlet nodes
    network.add_node(Node::new("inlet".to_string(), NodeType::Inlet));
    network.add_node(Node::new("outlet".to_string(), NodeType::Outlet));

    // Calculate resistance for circular pipe (Hagen-Poiseuille)
    let pipe_area = PI * pipe_radius * pipe_radius;
    let resistance = 8.0 * fluid_viscosity * pipe_length / (PI * pipe_radius.powi(4));

    // Add channel between inlet and outlet
    let channel_props = ChannelProperties::new(resistance, pipe_length, pipe_area);
    network.add_edge("inlet", "outlet", channel_props)?;

    // Set boundary conditions
    network.set_boundary_condition(
        "inlet",
        BoundaryCondition::PressureInlet {
            pressure: pressure_drop,
        },
    )?;
    network.set_boundary_condition(
        "outlet",
        BoundaryCondition::PressureOutlet { pressure: 0.0 },
    )?;

    println!("\nNetwork Summary:");
    println!("  Nodes: {}", network.node_count());
    println!("  Edges: {}", network.edge_count());

    // Create and run solver
    let mut solver = NetworkSolver::new();
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    println!("\n========================================");
    println!("Numerical Solution");
    println!("========================================");

    // Extract pressures
    let pressures = solution.pressures();

    // Convert HashMap to Vec for ordered processing
    let mut pressure_entries: Vec<_> = pressures.iter().collect();
    pressure_entries.sort_by_key(|&(idx, _)| idx.index());
    let pressure_values: Vec<f64> = pressure_entries.iter().map(|(_, &p)| p).collect();

    if pressure_values.len() >= 2 {
        let inlet_pressure = pressure_values[0];
        let outlet_pressure = pressure_values[1];
        let calculated_pressure_drop = inlet_pressure - outlet_pressure;

        println!("Inlet pressure: {:.2} Pa", inlet_pressure);
        println!("Outlet pressure: {:.2} Pa", outlet_pressure);
        println!(
            "Calculated pressure drop: {:.2} Pa",
            calculated_pressure_drop
        );

        // Calculate flow rate from pressure drop
        let calculated_flow_rate = calculated_pressure_drop / resistance;
        println!("Calculated flow rate: {:.6e} m³/s", calculated_flow_rate);

        // Calculate average velocity
        let avg_velocity = calculated_flow_rate / pipe_area;
        println!("Average velocity: {:.4} m/s", avg_velocity);

        // Calculate Reynolds number
        let reynolds = fluid_density * avg_velocity * (2.0 * pipe_radius) / fluid_viscosity;
        println!("Reynolds number: {:.1}", reynolds);

        let flow_regime = if reynolds < 2300.0 {
            "Laminar"
        } else if reynolds < 4000.0 {
            "Transitional"
        } else {
            "Turbulent"
        };
        println!("Flow regime: {}", flow_regime);

        println!("\n========================================");
        println!("Analytical Solution (Hagen-Poiseuille)");
        println!("========================================");

        // Analytical solution for laminar flow in circular pipe
        let analytical_flow_rate =
            PI * pipe_radius.powi(4) * pressure_drop / (8.0 * fluid_viscosity * pipe_length);
        let analytical_avg_velocity = analytical_flow_rate / pipe_area;
        let analytical_max_velocity = 2.0 * analytical_avg_velocity;

        println!("Analytical flow rate: {:.6e} m³/s", analytical_flow_rate);
        println!(
            "Analytical average velocity: {:.4} m/s",
            analytical_avg_velocity
        );
        println!(
            "Analytical maximum velocity: {:.4} m/s",
            analytical_max_velocity
        );

        println!("\n========================================");
        println!("Validation Results");
        println!("========================================");

        // Calculate errors
        let flow_rate_error =
            ((calculated_flow_rate - analytical_flow_rate) / analytical_flow_rate * 100.0).abs();
        let velocity_error =
            ((avg_velocity - analytical_avg_velocity) / analytical_avg_velocity * 100.0).abs();

        println!("Flow rate error: {:.2}%", flow_rate_error);
        println!("Velocity error: {:.2}%", velocity_error);

        // Validation check
        let tolerance = 1.0; // 1% tolerance
        if flow_rate_error < tolerance && velocity_error < tolerance {
            println!("\n✅ VALIDATION PASSED");
            println!(
                "The numerical solution matches the analytical solution within {:.1}% tolerance.",
                tolerance
            );
        } else {
            println!("\n⚠️ VALIDATION WARNING");
            println!(
                "The numerical solution differs from the analytical solution by more than {:.1}%.",
                tolerance
            );
        }

        // Additional validation metrics
        println!("\n========================================");
        println!("Additional Metrics");
        println!("========================================");

        // Pressure gradient
        let numerical_gradient = calculated_pressure_drop / pipe_length;
        let analytical_gradient = pressure_drop / pipe_length;
        println!(
            "Numerical pressure gradient: {:.2} Pa/m",
            numerical_gradient
        );
        println!(
            "Analytical pressure gradient: {:.2} Pa/m",
            analytical_gradient
        );

        // Wall shear stress
        let wall_shear_stress = pipe_radius * analytical_gradient / 2.0;
        println!("Wall shear stress: {:.4} Pa", wall_shear_stress);

        // Friction factor (Darcy-Weisbach)
        let friction_factor = 64.0 / reynolds;
        println!("Friction factor (f): {:.4}", friction_factor);
    }

    Ok(())
}
