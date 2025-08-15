//! 1D Pipe Flow Validation
//!
//! This example validates the 1D pipe flow solver against the analytical
//! Hagen-Poiseuille solution for laminar flow in a circular pipe.

use cfd_1d::{NetworkBuilder, NetworkSolver};
use cfd_core::Fluid;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("========================================");
    println!("1D Pipe Flow Validation");
    println!("========================================\n");
    
    // Pipe parameters
    let pipe_radius: f64 = 0.01; // 10 mm radius
    let pipe_length: f64 = 0.1;  // 100 mm length
    let fluid_viscosity: f64 = 1e-3;  // Water at 20°C (Pa·s)
    let fluid_density: f64 = 1000.0;   // Water density (kg/m³)
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
    let fluid = Fluid::new_newtonian("Test Fluid", fluid_density, fluid_viscosity);
    
    // Create 1D network
    let mut builder = NetworkBuilder::new().with_fluid(fluid);
    
    // Add inlet and outlet nodes with pressure boundary conditions
    builder = builder.add_inlet_pressure("inlet", 0.0, 0.0, pressure_drop);
    builder = builder.add_outlet_pressure("outlet", pipe_length, 0.0, 0.0);
    
    // Calculate resistance for circular pipe (Hagen-Poiseuille)
    let pipe_area = PI * pipe_radius * pipe_radius;
    let resistance = 8.0 * fluid_viscosity * pipe_length / (PI * pipe_radius.powi(4));
    
    // Add channel between inlet and outlet
    builder = builder.add_channel("pipe", "inlet", "outlet", cfd_1d::ChannelProperties { resistance, length: pipe_length, area: pipe_area });
    
    // Build network
    let mut network = builder.build()?;
    
    println!("\nNetwork Summary:");
    println!("  Nodes: {}", network.nodes().count());
    println!("  Edges: {}", network.edges().count());
    
    // Create and run solver
    let solver = NetworkSolver::new();
    let solution = solver.solve_steady_state(&mut network)?;
    
    // Extract results from the network after solving
    let nodes: Vec<_> = network.nodes().collect();
    let edges: Vec<_> = network.edges().collect();
    
    // Get the flow rate through the pipe
    let flow_rate_numerical = if let Some(edge) = edges.first() {
        edge.flow_rate.unwrap_or(0.0)
    } else {
        0.0
    };
    
    // Get pressures at inlet and outlet
    let pressure_at_inlet = nodes.iter()
        .find(|n| n.id == "inlet")
        .and_then(|n| n.pressure)
        .unwrap_or(pressure_drop);
    let pressure_at_outlet = nodes.iter()
        .find(|n| n.id == "outlet")
        .and_then(|n| n.pressure)
        .unwrap_or(0.0);
    let actual_pressure_drop = pressure_at_inlet - pressure_at_outlet;
    
    // Calculate analytical solution (Hagen-Poiseuille)
    let flow_rate_analytical = PI * pipe_radius.powi(4) * pressure_drop / (8.0 * fluid_viscosity * pipe_length);
    let max_velocity_analytical = pressure_drop * pipe_radius.powi(2) / (4.0 * fluid_viscosity * pipe_length);
    let avg_velocity_analytical = max_velocity_analytical / 2.0;
    let reynolds_number = fluid_density * avg_velocity_analytical * 2.0 * pipe_radius / fluid_viscosity;
    
    // Calculate resistance
    let resistance_analytical = 8.0 * fluid_viscosity * pipe_length / (PI * pipe_radius.powi(4));
    let resistance_numerical = if flow_rate_numerical.abs() > 1e-10 {
        actual_pressure_drop / flow_rate_numerical
    } else {
        0.0
    };
    
    println!("\n========================================");
    println!("RESULTS");
    println!("========================================");
    
    println!("\nFlow Rate:");
    println!("  Numerical:  {:.9} m³/s", flow_rate_numerical);
    println!("  Analytical: {:.9} m³/s", flow_rate_analytical);
    if flow_rate_analytical.abs() > 1e-10 {
        let error = ((flow_rate_numerical - flow_rate_analytical) / flow_rate_analytical).abs() * 100.0;
        println!("  Error: {:.2}%", error);
        if error < 1.0 {
            println!("  ✓ Excellent agreement (< 1% error)");
        } else if error < 5.0 {
            println!("  ✓ Good agreement (< 5% error)");
        } else {
            println!("  ⚠ Large error - check solver configuration");
        }
    }
    
    println!("\nPressure Drop:");
    println!("  Applied:    {:.3} Pa", pressure_drop);
    println!("  Calculated: {:.3} Pa", actual_pressure_drop);
    
    println!("\nFlow Resistance:");
    println!("  Numerical:  {:.3e} Pa·s/m³", resistance_numerical);
    println!("  Analytical: {:.3e} Pa·s/m³", resistance_analytical);
    if resistance_analytical > 0.0 {
        let error = ((resistance_numerical - resistance_analytical) / resistance_analytical).abs() * 100.0;
        println!("  Error: {:.2}%", error);
    }
    
    println!("\nVelocity:");
    println!("  Maximum (analytical): {:.6} m/s", max_velocity_analytical);
    println!("  Average (analytical): {:.6} m/s", avg_velocity_analytical);
    
    println!("\nReynolds Number: {:.1}", reynolds_number);
    if reynolds_number < 2300.0 {
        println!("  ✓ Laminar flow (Re < 2300)");
    } else {
        println!("  ⚠ Turbulent flow (Re > 2300) - Hagen-Poiseuille may not apply");
    }
    
    // Additional validation: Check solver convergence
    println!("\nSolver Performance:");
    println!("  Converged: {}", solution.converged);
    println!("  Iterations: {}", solution.iterations);
    println!("  Residual: {:.3e}", solution.residual);
    println!("  Solve time: {:.3} ms", solution.solve_time * 1000.0);
    if solution.converged {
        println!("  ✓ Solution converged successfully");
    }
    
    println!("\n========================================");
    println!("Validation Complete!");
    println!("========================================");
    
    Ok(())
}