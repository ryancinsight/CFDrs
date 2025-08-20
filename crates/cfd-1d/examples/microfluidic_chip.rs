//! Example demonstrating a complete microfluidic chip simulation using cfd-1d.
//!
//! This example creates a microfluidic chip with:
//! - Multiple inlets with different pressures
//! - A mixing chamber
//! - Flow sensors
//! - Valves for flow control
//! - Multiple outlets
//!
//! Run with: cargo run --example microfluidic_chip

use cfd_1d::{NetworkBuilder, NetworkSolver, ChannelProperties};
use cfd_core::fluid::Fluid;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧪 Microfluidic Chip Simulation");
    println!("================================");

    // Create a simple but functional microfluidic network
    let mut network = NetworkBuilder::new()
        .with_fluid(Fluid::water())

        // Simple T-junction network
        .add_inlet_pressure("inlet", 0.0, 0.0, 2000.0)  // 2 kPa
        .add_junction("junction", 2.0, 0.0)
        .add_outlet_pressure("outlet_1", 4.0, 1.0, 0.0)
        .add_outlet_pressure("outlet_2", 4.0, -1.0, 0.0)

        // Channels connecting the network with self-documenting properties
        .add_channel("input_ch", "inlet", "junction", ChannelProperties::new(100.0, 0.001, 100e-6))
        .add_channel("output_ch1", "junction", "outlet_1", ChannelProperties::new(200.0, 0.001, 100e-6))
        .add_channel("output_ch2", "junction", "outlet_2", ChannelProperties::new(200.0, 0.001, 100e-6))

        .build()?;

    println!("✅ Network created successfully!");
    println!("   - Nodes: {}", network.node_count());
    println!("   - Edges: {}", network.edge_count());
    println!("   - Fluid: {}", network.fluid().name);

    // Create solver with custom configuration using builder pattern
    let solver_config = cfd_core::SolverConfig::<f64>::builder()
        .max_iterations(1000)
        .tolerance(1e-6)
        .verbosity(2) // verbose = true means verbosity level 2
        .build(); // Clean, unambiguous method name
    
    let solver = NetworkSolver::with_config(solver_config);

    // Solve the flow problem
    println!("\n🔄 Solving flow equations...");
    let result = solver.solve_steady_state(&mut network)?;

    if result.converged {
        println!("✅ Solution converged in {} iterations", result.iterations);
        println!("   - Final residual: {:.2e}", result.residual);
        println!("   - Solve time: {:.3} ms", result.solve_time * 1000.0);
    } else {
        println!("❌ Solution did not converge");
        return Ok(());
    }

    // Display results
    println!("\n📊 Flow Analysis Results");
    println!("========================");

    // Node pressures
    println!("\n🔘 Node Pressures:");
    for node in network.nodes() {
        if let Some(pressure) = node.pressure {
            println!("   {}: {:.1} Pa", node.id, pressure);
        }
    }

    // Edge flow rates
    println!("\n➡️  Flow Rates:");
    let mut total_flow_in = 0.0;
    let mut total_flow_out = 0.0;

    for edge in network.edges() {
        if let Some(flow_rate) = edge.flow_rate {
            let flow_ml_min = flow_rate * 1e6 * 60.0; // Convert m³/s to mL/min
            println!("   {}: {:.3} mL/min ({:.2e} m³/s)", edge.id, flow_ml_min, flow_rate);

            // Track total flows
            match edge.id.as_str() {
                "input_ch" => total_flow_in += flow_rate,
                "output_ch1" | "output_ch2" => total_flow_out += flow_rate,
                _ => {}
            }
        }
    }

    // Flow conservation check
    println!("\n⚖️  Flow Conservation:");
    println!("   Total inflow:  {:.3} mL/min", total_flow_in * 1e6 * 60.0);
    println!("   Total outflow: {:.3} mL/min", total_flow_out * 1e6 * 60.0);
    let conservation_error = {
        let error: f64 = (total_flow_in - total_flow_out) / total_flow_in * 100.0;
        error.abs()
    };
    println!("   Conservation error: {:.2}%", conservation_error);

    // Pressure drops
    println!("\n📉 Pressure Drops:");
    for edge in network.edges() {
        if let Some(pressure_drop) = edge.pressure_drop {
            println!("   {}: {:.1} Pa", edge.id, pressure_drop);
        }
    }

    // Component analysis
    println!("\n🔧 Component Analysis:");

    // Channel resistances
    if let Some(ch1) = network.get_edge("output_ch1") {
        println!("   Output channel 1 resistance: {:.1e} Pa·s/m³", ch1.resistance);
    }

    if let Some(ch2) = network.get_edge("output_ch2") {
        println!("   Output channel 2 resistance: {:.1e} Pa·s/m³", ch2.resistance);
    }

    // Flow splitting analysis
    println!("\n🌀 Flow Splitting Analysis:");
    let input_flow = network.get_edge("input_ch").and_then(|e| e.flow_rate).unwrap_or(0.0);
    let output1_flow = network.get_edge("output_ch1").and_then(|e| e.flow_rate).unwrap_or(0.0);
    let output2_flow = network.get_edge("output_ch2").and_then(|e| e.flow_rate).unwrap_or(0.0);

    if input_flow > 0.0 {
        let ratio1 = output1_flow / input_flow * 100.0;
        let ratio2 = output2_flow / input_flow * 100.0;
        println!("   Output 1 split: {:.1}%", ratio1);
        println!("   Output 2 split: {:.1}%", ratio2);
    }

    // Performance metrics
    println!("\n📈 Performance Metrics:");
    let main_flow = network.get_edge("input_ch").and_then(|e| e.flow_rate).unwrap_or(0.0);
    let throughput = main_flow * 1e6 * 3600.0; // mL/hour
    println!("   Chip throughput: {:.1} mL/hour", throughput);
    
    // Reynolds number estimation (assuming circular channel)
    let fluid = network.fluid();
    let diameter = 200e-6; // 200 μm diameter
    let velocity = main_flow / (std::f64::consts::PI * (diameter/2.0_f64).powi(2));
    let reynolds = fluid.density * velocity * diameter / fluid.dynamic_viscosity(20.0);
    println!("   Reynolds number: {:.2}", reynolds);
    
    if reynolds < 1.0 {
        println!("   Flow regime: Stokes flow (very low Re)");
    } else if reynolds < 100.0 {
        println!("   Flow regime: Laminar flow");
    } else {
        println!("   Flow regime: Transitional/Turbulent");
    }

    println!("\n🎉 Simulation completed successfully!");
    
    Ok(())
}
