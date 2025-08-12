//! 1D Pipe Flow Validation
//!
//! This example validates the 1D pipe flow solver against the analytical
//! Hagen-Poiseuille solution for laminar flow in a circular pipe.

use cfd_1d::{NetworkBuilder, NetworkSolver, SolverConfig, Channel, ChannelShape};
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("========================================");
    println!("1D Pipe Flow Validation");
    println!("========================================\n");
    
    // Pipe parameters
    let pipe_radius = 0.01; // 10 mm radius
    let pipe_length = 0.1;  // 100 mm length
    let fluid_viscosity = 1e-3;  // Water at 20°C (Pa·s)
    let fluid_density = 1000.0;   // Water density (kg/m³)
    let pressure_drop = 10.0; // Pressure drop across pipe (Pa)
    
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
    
    // Create 1D network
    let mut builder = NetworkBuilder::new();
    
    // Add inlet and outlet nodes
    let inlet = builder.add_node("inlet", 0.0, 0.0, 0.0);
    let outlet = builder.add_node("outlet", pipe_length, 0.0, 0.0);
    
    // Create circular channel
    let channel = Channel {
        id: 0,
        inlet_node: inlet,
        outlet_node: outlet,
        shape: ChannelShape::Circular { diameter: 2.0 * pipe_radius },
        length: pipe_length,
        roughness: Some(0.0), // Smooth pipe
    };
    
    // Add channel to network
    builder.add_channel(channel);
    
    // Set boundary conditions
    builder.set_pressure_bc(inlet, pressure_drop);
    builder.set_pressure_bc(outlet, 0.0);
    
    // Build network
    let network = builder.build()?;
    
    println!("\nNetwork Summary:");
    println!("  Nodes: {}", network.nodes.len());
    println!("  Channels: {}", network.channels.len());
    
    // Configure solver
    let config = SolverConfig {
        tolerance: 1e-8,
        max_iterations: 1000,
        fluid_viscosity,
        fluid_density,
        use_entrance_effects: false, // Disable for comparison with fully developed flow
    };
    
    // Create and run solver
    let mut solver = NetworkSolver::new(config);
    let solution = solver.solve(&network)?;
    
    // Extract results
    let flow_rate_numerical = solution.flow_rates.get(&0).copied().unwrap_or(0.0);
    let pressure_at_inlet = solution.pressures.get(&inlet).copied().unwrap_or(0.0);
    let pressure_at_outlet = solution.pressures.get(&outlet).copied().unwrap_or(0.0);
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
    
    // Additional validation: Check continuity
    println!("\nContinuity Check:");
    let inlet_flow = solution.flow_rates.values().filter(|&&q| q > 0.0).sum::<f64>();
    let outlet_flow = solution.flow_rates.values().filter(|&&q| q < 0.0).sum::<f64>().abs();
    let continuity_error = (inlet_flow - outlet_flow).abs();
    println!("  Inlet flow:  {:.9} m³/s", inlet_flow);
    println!("  Outlet flow: {:.9} m³/s", outlet_flow);
    println!("  Error: {:.3e} m³/s", continuity_error);
    if continuity_error < 1e-10 {
        println!("  ✓ Mass conservation satisfied");
    }
    
    println!("\n========================================");
    println!("Validation Complete!");
    println!("========================================");
    
    Ok(())
}