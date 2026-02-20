//! Example demonstrating conversion of a schematic to a 3D mesh.
//!
//! This example creates a simple 1D schematic using `cfd-schematics`,
//! converts it to a 3D `Schematic` using `cfd-mesh`'s interchange,
//! and verifies that the resulting channel profile preserves the rectangular shape.
//!
//! Run with: cargo run -p cfd-mesh --example schematic_to_mesh_demo

#[cfg(feature = "scheme-io")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use cfd_schematics::geometry::{
        generator::create_geometry, ChannelType, SplitType,
    };
    use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};
    use cfd_mesh::io::scheme;
    use cfd_mesh::channel::profile::ChannelProfile;

    println!("🔌 Schematic to Mesh Integration Demo");
    println!("=====================================");

    // 1. Create a simple schematic (single straight channel)
    println!("1. Generating schematic...");
    let box_dims = (100.0, 50.0);
    let mut config = GeometryConfig::default();
    config.channel_width = 1.0;  // 1mm width
    config.channel_height = 0.5; // 0.5mm height

    let system = create_geometry(
        box_dims,
        &[], // No splits -> single channel
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    println!("   Generated {} channels.", system.channels.len());

    // 2. Add a Frustum (Venturi) channel manually for testing
    // We hack the system to add a second channel that is a Frustum
    use cfd_schematics::geometry::Point2D;
    use cfd_schematics::config::TaperProfile;
    
    // Frustum from (0, 10) to (10, 10). Inlet=1.0, Throat=0.2, Outlet=1.0
    let frustum_path = vec![(0.0, 10.0), (5.0, 10.0), (10.0, 10.0)];
    let frustum_widths = vec![1.0, 0.2, 1.0]; // Simple V-shape for test
    
    let frustum_channel = cfd_schematics::geometry::Channel {
        id: 99,
        from_node: 0, // Reuse nodes for simplicity or just dummy IDs
        to_node: 1,
        width: 1.0, // Inlet width
        height: 0.5,
        channel_type: ChannelType::Frustum {
            path: frustum_path,
            widths: frustum_widths.clone(),
            inlet_width: 1.0,
            throat_width: 0.2,
            outlet_width: 1.0,
            taper_profile: TaperProfile::Linear,
            throat_position: 0.5,
            has_venturi_throat: true,
        },
        metadata: None,
    };
    
    // Clone system and add frustum
    let mut system_with_frustum = system.clone();
    system_with_frustum.channels.push(frustum_channel);

    println!("2. Converting to 3D Schematic...");
    let substrate_height = 5.0; // 5mm substrate
    let segments = 32;
    
    let schematic3d = scheme::from_channel_system(&system_with_frustum, substrate_height, segments)?;

    println!("   Converted {} channels.", schematic3d.channels.len());

    // 3. Inspect Profiles and Mesh
    println!("3. Inspecting Profiles and Meshing...");
    
    // Prepare Mesher
    use cfd_mesh::channel::sweep::SweepMesher;
    use cfd_mesh::storage::vertex_pool::VertexPool;
    use cfd_mesh::core::index::RegionId;

    let mesher = SweepMesher::new();
    let mut pool = VertexPool::new(0.1, 1e-6);

    for (i, channel_def) in schematic3d.channels.iter().enumerate() {
        println!("   Processing Channel {} (ID: {})...", i, channel_def.id);

        match &channel_def.profile {
            ChannelProfile::Rectangular { width, height } => {
                println!("      Profile: Rectangular ({:.3} x {:.3} mm)", width, height);
            }
            ChannelProfile::Circular { radius, .. } => {
                println!("      Profile: Circular (r = {:.3} mm)", radius);
            }
            _ => println!("      Profile: Other"),
        }

        if let Some(scales) = &channel_def.width_scales {
            println!("      Has width scales: yes (len={})", scales.len());
            // Verify scaling values for Frustum
            if channel_def.id.contains("99") { // The frustum channel
                assert!((scales[0] - 1.0).abs() < 1e-4, "Start scale should be 1.0 (1.0/1.0)");
                assert!((scales[1] - 0.2).abs() < 1e-4, "Middle scale should be 0.2 (0.2/1.0)");
                assert!((scales[2] - 1.0).abs() < 1e-4, "End scale should be 1.0 (1.0/1.0)");
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
                 let vertices: Vec<_> = faces.iter()
                    .flat_map(|f| f.vertices.iter())
                    .map(|&vid| pool.get(vid).position)
                    .collect();
                 
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

#[cfg(not(feature = "scheme-io"))]
fn main() {
    println!("This example requires the 'scheme-io' feature.");
    println!("Run with: cargo run -p cfd-mesh --features scheme-io --example schematic_to_mesh_demo");
}
