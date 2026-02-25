#[test]
fn test_schematics2mesh_stl_output() {
    use std::fs::File;
    use std::path::PathBuf;
    use cfd_mesh::application::channel::profile::ChannelProfile;
    use cfd_mesh::application::channel::sweep::SweepMesher;
    use cfd_mesh::domain::core::index::RegionId;
    use cfd_mesh::domain::mesh::IndexedMesh;
    use cfd_mesh::infrastructure::io::scheme;
    use cfd_mesh::infrastructure::io::stl::write_stl_binary;
    use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
    use cfd_schematics::geometry::{generator::create_geometry, SplitType};

    // 1. Generate 2D geometry (a simple symmetric bifurcation)
    let box_dims = (20.0, 10.0);
    let mut config = GeometryConfig::default();
    config.channel_width = 4.0; // 4mm
    config.channel_height = 4.0;

    let system = create_geometry(
        box_dims,
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    // 2. Convert to 3D Schematic representation
    let substrate_height = 5.0; // 5mm substrate
    let segments = 32;
    let schematic3d = scheme::from_channel_system(&system, substrate_height, segments).unwrap();

    // 3. Sweep and collect into an IndexedMesh
    let mesher = SweepMesher::new();
    let mut mesh = IndexedMesh::new();

    for channel_def in &schematic3d.channels {
        if let Some(scales) = &channel_def.width_scales {
            let faces = mesher.sweep_variable(
                &channel_def.profile,
                &channel_def.path,
                scales,
                &mut mesh.vertices,
                RegionId::new(0),
            );
            for face in faces {
                mesh.faces.push(face);
            }
        } else {
            let faces = mesher.sweep(
                &channel_def.profile,
                &channel_def.path,
                &mut mesh.vertices,
                RegionId::new(0),
            );
            for face in faces {
                mesh.faces.push(face);
            }
        }
    }

    assert!(mesh.faces.len() > 0, "Mesh should have faces");

    // 4. Write to STL
    let mut stl_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    stl_path.push("outputs");
    std::fs::create_dir_all(&stl_path).unwrap();
    stl_path.push("schematics2mesh_bifurcation.stl");

    let mut file = File::create(&stl_path).unwrap();
    write_stl_binary(&mut file, &mesh).unwrap();
    
    // Check that STL exists and isn't empty
    let metadata = std::fs::metadata(&stl_path).unwrap();
    assert!(metadata.len() > 84); // Header + triangle count
}
