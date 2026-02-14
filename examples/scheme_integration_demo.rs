//! Example demonstrating scheme → cfd-1d bridge.
//!
//! Run with: `cargo run --example scheme_integration_demo --features scheme-integration`

#[cfg(not(feature = "scheme-integration"))]
fn main() {
    println!("Scheme Integration Demo");
    println!("======================");
    println!();
    println!("The 'scheme-integration' feature is not enabled.");
    println!();
    println!("Run with the feature enabled:");
    println!("  cargo run --example scheme_integration_demo --features scheme-integration");
    println!();
    println!("Alternatively, design microfluidic networks programmatically");
    println!("using the NetworkBuilder API without 2D schematic input.");
}

#[cfg(feature = "scheme-integration")]
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    use cfd_1d::scheme_bridge::SchemeNetworkConverter;
    use scheme::{
        config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig},
        geometry::{generator::create_geometry, SplitType},
    };

    println!("Scheme → CFD-1D Bridge Demo");
    println!("===========================\n");

    // ── 1. Simple bifurcation ────────────────────────────────────────────
    println!("1. Creating a simple bifurcation schematic...");
    let bifurcation_system = create_geometry(
        (200.0, 100.0),
        &[SplitType::Bifurcation],
        &GeometryConfig::default(),
        &ChannelTypeConfig::AllStraight,
    );
    println!(
        "   Scheme: {} nodes, {} channels",
        bifurcation_system.nodes.len(),
        bifurcation_system.channels.len()
    );

    // Convert to cfd-1d network
    println!("\n2. Converting to CFD-1D network...");
    let converter = SchemeNetworkConverter::new(&bifurcation_system);
    let summary = converter.summary();
    println!("{}", summary);

    let network = converter.build_network_with_water()?;
    println!(
        "   Network: {} nodes, {} edges",
        network.node_count(),
        network.edge_count()
    );

    // Show node info
    println!("\n3. Node details:");
    for node in network.nodes() {
        println!(
            "   {} ({:?}) at ({:.4e}, {:.4e}) m",
            node.id, node.node_type, node.position.0, node.position.1
        );
    }

    // ── 2. Complex serpentine system ─────────────────────────────────────
    println!("\n4. Creating complex serpentine bifurcation+trifurcation...");
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

    // Convert with millimetre scale (scheme coordinates in mm)
    let converter = SchemeNetworkConverter::with_scale(&complex_system, 1e-3);
    let summary = converter.summary();
    println!("{}", summary);

    let network = converter.build_network_with_water()?;
    println!(
        "   Network: {} nodes, {} edges",
        network.node_count(),
        network.edge_count()
    );

    // Show resistance info
    println!("\n5. Edge resistance summary:");
    for edge_ref in network.graph.edge_references() {
        let e = edge_ref.weight();
        println!(
            "   {} ({:?}): R = {:.4e} Pa·s/m³, k = {:.4e}",
            e.id, e.edge_type, e.resistance, e.quad_coeff
        );
    }

    // ── 3. JSON round-trip ──────────────────────────────────────────────
    println!("\n6. JSON round-trip test...");
    let json = bifurcation_system.to_json()?;
    let restored: scheme::geometry::ChannelSystem =
        scheme::geometry::ChannelSystem::from_json(&json)?;
    let conv2 = SchemeNetworkConverter::new(&restored);
    let net2 = conv2.build_network_with_water()?;
    assert_eq!(network.node_count(), net2.node_count());
    println!("   Round-trip OK: {} nodes, {} edges", net2.node_count(), net2.edge_count());

    println!("\nScheme → CFD-1D bridge demo completed successfully!");
    Ok(())
}
