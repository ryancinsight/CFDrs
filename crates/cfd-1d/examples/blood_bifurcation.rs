//! 1D Blood Flow in a Bifurcation (Carreau-Yasuda Model)
//!
//! This example demonstrates simulating blood flow in a bifurcating artery
//! using the Carreau-Yasuda non-Newtonian viscosity model.
//!
//! Geometry:
//! - Parent vessel: D=4mm, L=20mm
//! - Daughter vessel 1: D=3mm, L=15mm
//! - Daughter vessel 2: D=2.5mm, L=15mm
//!
//! Conditions:
//! - Inlet Flow: 5 mL/s (pulsatile average)
//! - Outlet Pressure: 100 mmHg

use cfd_1d::channel::{ChannelGeometry, ChannelType, CrossSection, SurfaceProperties, Wettability};
use cfd_1d::network::{Network, NetworkBuilder};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::fluid::non_newtonian::CarreauYasuda;
use cfd_core::physics::fluid::FluidTrait;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Define Blood Properties (Carreau-Yasuda)
    let blood = CarreauYasuda::<f64>::blood();
    println!("Fluid: {}", blood.name());
    println!("Density: {} kg/m^3", blood.density);
    println!("Viscosity (inf): {} Pa.s", blood.viscosity_inf);
    println!("Viscosity (0): {} Pa.s", blood.viscosity_zero);

    // 2. Build Network Topology
    let mut builder = NetworkBuilder::<f64>::new();

    let inlet_node = builder.add_inlet("Inlet".to_string());
    let junction_node = builder.add_junction("Bifurcation".to_string());
    let outlet1_node = builder.add_outlet("Outlet1".to_string());
    let outlet2_node = builder.add_outlet("Outlet2".to_string());

    // Parent Artery
    let parent_edge = builder.connect_with_pipe(inlet_node, junction_node, "Parent".to_string());

    // Daughter Arteries
    let daughter1_edge = builder.connect_with_pipe(junction_node, outlet1_node, "Daughter1".to_string());
    let daughter2_edge = builder.connect_with_pipe(junction_node, outlet2_node, "Daughter2".to_string());

    let graph = builder.build()?;

    // 3. Create Network with Fluid
    let mut network = Network::new(graph, blood.clone());

    // Helper to create channel geometry
    let create_geometry = |diameter: f64, length: f64| -> ChannelGeometry<f64> {
        ChannelGeometry {
            channel_type: ChannelType::Straight,
            length,
            cross_section: CrossSection::Circular { diameter },
            surface: SurfaceProperties {
                roughness: 0.0,
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }
    };

    // 4. Assign Geometry (Diameter, Length)
    // Parent: 4mm diameter, 5cm length (L/D > 10)
    network.add_edge_properties(parent_edge, cfd_1d::network::EdgeProperties {
        id: "Parent".to_string(),
        component_type: cfd_1d::network::ComponentType::Pipe,
        length: 0.05,
        area: std::f64::consts::PI * (0.002 * 0.002), // r=2mm
        hydraulic_diameter: Some(0.004),
        resistance: 0.0, // Calculated automatically
        geometry: Some(create_geometry(0.004, 0.05)),
        properties: std::collections::HashMap::new(),
    });

    // Daughter 1: 3mm diameter, 4cm length
    network.add_edge_properties(daughter1_edge, cfd_1d::network::EdgeProperties {
        id: "Daughter1".to_string(),
        component_type: cfd_1d::network::ComponentType::Pipe,
        length: 0.04,
        area: std::f64::consts::PI * (0.0015 * 0.0015),
        hydraulic_diameter: Some(0.003),
        resistance: 0.0,
        geometry: Some(create_geometry(0.003, 0.04)),
        properties: std::collections::HashMap::new(),
    });

    // Daughter 2: 2.5mm diameter, 3cm length
    network.add_edge_properties(daughter2_edge, cfd_1d::network::EdgeProperties {
        id: "Daughter2".to_string(),
        component_type: cfd_1d::network::ComponentType::Pipe,
        length: 0.03,
        area: std::f64::consts::PI * (0.00125 * 0.00125),
        hydraulic_diameter: Some(0.0025),
        resistance: 0.0,
        geometry: Some(create_geometry(0.0025, 0.03)),
        properties: std::collections::HashMap::new(),
    });

    // 5. Boundary Conditions
    // Inlet Flow: 5 mL/s = 5e-6 m^3/s
    network.set_neumann_flow(inlet_node, 5e-6);

    // Outlet Pressures: 100 mmHg = 13332.2 Pa (relative to ambient? Or absolute)
    // Let's use 0 Pa gauge pressure at outlets for simplicity (driving pressure)
    // Or realistic: P_out = 13332.2 Pa
    let p_out = 13332.2;
    network.set_pressure(outlet1_node, p_out);
    network.set_pressure(outlet2_node, p_out);

    // Initial guess for flow rates (needed for non-Newtonian viscosity start)
    network.set_flow_rate(parent_edge, 5e-6);
    network.set_flow_rate(daughter1_edge, 3e-6);
    network.set_flow_rate(daughter2_edge, 2e-6);

    // Initial resistance update
    network.update_resistances()?;

    // 6. Setup Problem and Solver
    let problem = NetworkProblem::new(network);

    let config = SolverConfig {
        tolerance: 1e-6,
        max_iterations: 100, // Non-linear iterations
    };
    let solver = NetworkSolver::<f64, CarreauYasuda<f64>>::with_config(config);

    // 7. Solve
    println!("Solving network...");
    let solution = solver.solve(&problem)?;

    // 8. Analyze Results
    println!("Solution Converged!");
    println!("Node Pressures:");
    for (idx, p) in solution.pressures() {
        // Find node name
        let node = solution.graph.node_weight(*idx).unwrap();
        println!("  {}: {:.2} Pa ({:.2} mmHg)", node.id, p, p / 133.322);
    }

    println!("Edge Flow Rates:");
    for (idx, q) in solution.flow_rates() {
        let edge = solution.graph.edge_weight(*idx).unwrap();
        // Convert to mL/s
        let q_ml = q * 1e6;

        // Calculate shear rate and viscosity
        let props = solution.properties.get(idx).unwrap();

        let v = (*q).abs() / props.area;

        let shear_rate = if let Some(geometry) = &props.geometry {
             match geometry.cross_section {
                 CrossSection::Circular { diameter } => {
                     8.0 * v / diameter
                 },
                 _ => 0.0,
             }
        } else {
            0.0
        };

        let viscosity = blood.viscosity_at_shear(shear_rate, 310.15, 101325.0)?; // 37°C

        println!("  {}: {:.4} mL/s", edge.id, q_ml);
        println!("     Velocity: {:.2} cm/s", v * 100.0);
        println!("     Shear Rate: {:.1} 1/s", shear_rate);
        println!("     Apparent Viscosity: {:.2} mPa.s", viscosity * 1000.0);

        // Compare with Newtonian (high shear limit)
        let mu_inf = blood.viscosity_inf;
        let ratio = viscosity / mu_inf;
        println!("     Non-Newtonian Factor: {:.2}x (vs high-shear limit)", ratio);
    }

    Ok(())
}
