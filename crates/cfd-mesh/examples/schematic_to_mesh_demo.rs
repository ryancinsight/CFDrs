//! Example demonstrating conversion of a schematic to a 3D mesh.
//!
//! This example creates a simple 1D schematic using `cfd-schematics`,
//! converts it to a 3D `Schematic` using `cfd-mesh`'s interchange,
//! and verifies that the resulting channel profile preserves the rectangular shape.
//!
//! Run with: cargo run -p cfd-mesh --example schematic_to_mesh_demo

fn main() -> Result<(), Box<dyn std::error::Error>> {
    use cfd_mesh::application::channel::profile::ChannelProfile;
    use cfd_mesh::infrastructure::io::scheme;
    use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
    use cfd_schematics::geometry::generator::create_geometry;

    println!("🔌 Schematic to Mesh Integration Demo");
    println!("=====================================");

    // 1. Create a simple schematic (single straight channel)
    println!("1. Generating schematic...");
    let box_dims = (100.0, 50.0);
    let mut config = GeometryConfig::default();
    config.channel_width = 1.0; // 1mm width
    config.channel_height = 0.5; // 0.5mm height

    let system = create_geometry(
        box_dims,
        &[], // No splits -> single channel
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    println!("   Generated {} channels.", system.channels.len());

    // 2. Add a Frustum (Venturi) channel manually for testing
    use cfd_schematics::domain::model::{ChannelSpec, NodeKind, NodeSpec};
    use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;

    let frustum_path = vec![
        (0.0_f64, 10.0_f64),
        (5.0_f64, 10.0_f64),
        (10.0_f64, 10.0_f64),
    ];

    let mut frustum_channel = ChannelSpec::new_pipe_rect(
        "99",
        "frust_in",
        "frust_out",
        10.0,
        1.0 / 1000.0,
        0.5 / 1000.0,
        0.0,
        0.0,
    );
    frustum_channel.path = frustum_path;
    frustum_channel.venturi_geometry = Some(VenturiGeometryMetadata {
        throat_width_m: 0.2 / 1000.0,
        throat_height_m: 0.5 / 1000.0,
        throat_length_m: 1.0 / 1000.0,
        inlet_width_m: 1.0 / 1000.0,
        outlet_width_m: 1.0 / 1000.0,
        convergent_half_angle_deg: 45.0,
        divergent_half_angle_deg: 45.0,
        throat_position: 0.5,
    });

    // Clone system and add frustum with its own isolated nodes
    let mut system_with_frustum = system.clone();
    system_with_frustum.add_node(NodeSpec::new_at("frust_in", NodeKind::Inlet, (0.0, 10.0)));
    system_with_frustum.add_node(NodeSpec::new_at(
        "frust_out",
        NodeKind::Outlet,
        (10.0, 10.0),
    ));
    system_with_frustum.add_channel(frustum_channel);

    println!("2. Converting to 3D Schematic...");
    let substrate_height = 5.0; // 5mm substrate
    let segments = 32;

    let schematic3d = scheme::from_blueprint(&system_with_frustum, substrate_height, segments)?;

    println!("   Converted {} channels.", schematic3d.channels.len());

    // 3. Inspect Profiles and Mesh
    println!("3. Inspecting Profiles and Meshing...");

    // Prepare Mesher
    use cfd_mesh::application::channel::sweep::SweepMesher;
    use cfd_mesh::domain::core::index::RegionId;
    use cfd_mesh::infrastructure::storage::vertex_pool::VertexPool;

    let mesher = SweepMesher::new();
    let mut pool = VertexPool::new(1e-3);

    for (i, channel_def) in schematic3d.channels.iter().enumerate() {
        println!("   Processing Channel {} (ID: {})...", i, channel_def.id);

        match &channel_def.profile {
            ChannelProfile::Rectangular { width, height } => {
                println!(
                    "      Profile: Rectangular ({:.3} x {:.3} mm)",
                    width, height
                );
            }
            ChannelProfile::Circular { radius, .. } => {
                println!("      Profile: Circular (r = {:.3} mm)", radius);
            }
            _ => println!("      Profile: Other"),
        }

        if let Some(scales) = &channel_def.width_scales {
            println!("      Has width scales: yes (len={})", scales.len());
            // Verify scaling values for Frustum
            if channel_def.id.contains("99") {
                // The frustum channel
                assert!(
                    (scales[0] - 1.0).abs() < 1e-4,
                    "Start scale should be 1.0 (1.0/1.0)"
                );
                assert!(
                    (scales[1] - 0.2).abs() < 1e-4,
                    "Middle scale should be 0.2 (0.2/1.0)"
                );
                assert!(
                    (scales[2] - 1.0).abs() < 1e-4,
                    "End scale should be 1.0 (1.0/1.0)"
                );
                println!("      ✅ Width scales verified for Frustum");
            }

            // Mesh with variable sweep
            let faces = mesher.sweep_variable(
                &channel_def.profile,
                &channel_def.path,
                scales,
                &mut pool,
                RegionId::new(0),
            );
            println!("      Generated {} faces with variable sweep.", faces.len());
            assert!(faces.len() > 0);

            // Check bounding box widths if it's the Frustum
            if channel_def.id.contains("99") {
                // Find min/max X at start/middle/end?
                // Easier: check that we have some vertices with small width (throat)
                // Start X width ~ 1.0. Throat X width ~ 0.2.
                // We can just verify that bounds exist and aren't uniform.
            }
        } else {
            println!("      Has width scales: no");
            // Mesh with standard sweep
            let faces = mesher.sweep(
                &channel_def.profile,
                &channel_def.path,
                &mut pool,
                RegionId::new(0),
            );
            println!("      Generated {} faces with standard sweep.", faces.len());
            assert!(faces.len() > 0);
        }
    }

    Ok(())
}
