//! Example demonstrating scheme -> cfd-1d bridge.
//!
//! Shows two approaches:
//! 1. Blueprint-based: Use preset blueprints (serpentine_chain, venturi_rect)
//!    and feed them through NetworkGenerationService + NetworkBuilderSink.
//! 2. Geometry-based: Use create_geometry to inspect schematic topology.
//!
//! Run with: `cargo run --example scheme_integration_demo`

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    use cfd_1d::network::Network;
    use cfd_1d::solver::{NetworkProblem, NetworkSolver};
    use cfd_1d::NetworkBuilderSink;
    use cfd_core::compute::solver::Solver;
    use cfd_core::physics::fluid::ConstantPropertyFluid;
    use cfd_schematics::application::use_cases::NetworkGenerationService;
    use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig};
    use cfd_schematics::geometry::types::SplitType;
    use cfd_schematics::geometry::create_geometry;
    use cfd_schematics::{serpentine_chain, symmetric_bifurcation};
    use petgraph::visit::EdgeRef;

    println!("Scheme -> CFD-1D Bridge Demo");
    println!("===========================\n");

    // -- 1. Simple bifurcation via blueprint ----------------------------------
    println!("1. Creating a simple bifurcation blueprint...");
    let blueprint = symmetric_bifurcation("bif", 1e-3, 1e-3, 0.5e-3, 0.4e-3);

    let water = ConstantPropertyFluid::<f64>::water_20c()?;
    let sink = NetworkBuilderSink::<f64, _>::new(water.clone());
    let service = NetworkGenerationService::new(sink);
    let mut network: Network<f64, ConstantPropertyFluid<f64>> = service.generate(&blueprint)?;

    println!(
        "   Network: {} nodes, {} edges",
        network.node_count(),
        network.edge_count()
    );

    // Show node info
    println!("\n2. Node details:");
    for node in network.nodes() {
        println!(
            "   {} ({:?}) at ({:.4e}, {:.4e}) m",
            node.id, node.node_type, node.position.0, node.position.1
        );
    }

    // Set boundary conditions and solve
    let inlet = network
        .graph
        .node_indices()
        .find(|idx| network.graph[*idx].id == "inlet")
        .expect("inlet");
    let outlet = network
        .graph
        .node_indices()
        .find(|idx| {
            let nid = &network.graph[*idx].id;
            nid.contains("outlet") || nid.contains("out")
        })
        .expect("outlet");
    network.set_pressure(inlet, 1000.0);
    network.set_pressure(outlet, 0.0);

    let solver = NetworkSolver::new();
    let problem = NetworkProblem::new(network);
    let solution = solver.solve(&problem)?;

    println!("\n3. Solution:");
    for (node_idx, &pressure) in solution.pressures() {
        println!("   Node {}: {:.2} Pa", node_idx.index(), pressure);
    }

    // -- 2. Complex serpentine geometry inspection ----------------------------
    println!("\n4. Creating complex serpentine bifurcation + trifurcation geometry...");
    let serpentine_config = SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        ..SerpentineConfig::default()
    };

    let complex_system = create_geometry(
        (300.0, 150.0),
        &[SplitType::Bifurcation, SplitType::Trifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllSerpentine(serpentine_config),
    );
    println!(
        "   Scheme: {} nodes, {} channels",
        complex_system.nodes.len(),
        complex_system.channels.len()
    );

    // -- 3. JSON round-trip ---------------------------------------------------
    println!("\n5. JSON round-trip test...");
    let json = complex_system.to_json()?;
    let restored = cfd_schematics::geometry::types::ChannelSystem::from_json(&json)?;
    assert_eq!(complex_system.nodes.len(), restored.nodes.len());
    assert_eq!(complex_system.channels.len(), restored.channels.len());
    println!(
        "   Round-trip OK: {} nodes, {} channels",
        restored.nodes.len(),
        restored.channels.len()
    );

    // -- 4. Serpentine-chain blueprint via NetworkBuilderSink ------------------
    println!("\n6. Serpentine-chain blueprint -> solver...");
    let serp_bp = serpentine_chain("serp_demo", 3, 0.01, 0.001);
    let sink2 = NetworkBuilderSink::<f64, _>::new(water.clone());
    let service2 = NetworkGenerationService::new(sink2);
    let mut net2: Network<f64, _> = service2.generate(&serp_bp)?;
    println!(
        "   Network: {} nodes, {} edges",
        net2.node_count(),
        net2.edge_count()
    );

    // Edge resistance summary
    println!("\n7. Edge resistance summary:");
    for edge_ref in net2.graph.edge_references() {
        let e = edge_ref.weight();
        println!(
            "   {} ({:?}): R = {:.4e} Pa.s/m^3",
            e.id, e.edge_type, e.resistance
        );
    }

    println!("\nScheme -> CFD-1D bridge demo completed successfully!");
    Ok(())
}
