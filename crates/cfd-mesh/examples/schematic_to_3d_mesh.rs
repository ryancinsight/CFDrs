use cfd_mesh::application::channel::sweep::SweepMesher;
use cfd_mesh::domain::core::index::RegionId;
use cfd_mesh::domain::mesh::IndexedMesh;
use cfd_mesh::infrastructure::io::scheme;
use cfd_mesh::infrastructure::io::stl::write_stl_binary;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Schematic to 3D Mesh Converter");
    println!("==============================");

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    // Since we're in cfd-mesh, the sibling is at ../cfd-schematics
    let json_path = manifest_dir.join("../cfd-schematics/outputs/mirrored_bifurcation/schematic.json");
    
    if !json_path.exists() {
        println!("❌ Missing input file: {}", json_path.display());
        println!("Please run the cfd-schematics mirrored_bifurcation example first:");
        println!("  cargo run --example mirrored_bifurcation -p cfd-schematics");
        return Ok(());
    }

    println!("Loading JSON: {}", json_path.display());

    let mut file = File::open(&json_path)?;
    let mut json_str = String::new();
    file.read_to_string(&mut json_str)?;

    let interchange: cfd_schematics::geometry::InterchangeChannelSystem =
        serde_json::from_str(&json_str)?;

    let substrate_height = 5.0_f64; // 5mm depth
    let segments = 32;              // Cross-section segments
    let schematic3d = scheme::from_interchange(
        &interchange,
        substrate_height as cfd_mesh::domain::core::scalar::Real,
        segments,
    )?;

    use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
    use cfd_mesh::domain::geometry::primitives::{Cube, PrimitiveMesh};
    use cfd_mesh::domain::core::scalar::Point3r;

    // Use the actual schematic bounding box, centered in Z.
    let (bw, bd) = interchange.box_dims_mm;
    let bw = bw as f64;
    let bd = bd as f64;
    let half_h = substrate_height / 2.0;

    println!("Constructing substrate: {}x{}x{}mm...", bw, bd, substrate_height);

    // Build a bounding block matching exact schematic XY, centered at Z=0.
    // Extend slightly in all directions so that over-piercing tubes are trimmed.
    let substrate = Cube {
        origin: Point3r::new(-5.0, -5.0, -half_h),
        width: (bw + 10.0) as _,
        height: (bd + 10.0) as _,
        depth: substrate_height as _,
    }
    .build()?;

    let mesher = SweepMesher::new();
    let mut final_solid = substrate;
    let mut all_channels = IndexedMesh::new();

    for channel_def in &schematic3d.channels {
        println!("Meshing channel {}...", channel_def.id);
        let mut current_ch = IndexedMesh::new();
        if let Some(scales) = &channel_def.width_scales {
            let faces = mesher.sweep_variable(
                &channel_def.profile,
                &channel_def.path,
                scales,
                &mut current_ch.vertices,
                RegionId::new(0),
            );
            for face in faces { current_ch.faces.push(face); }
        } else {
            let faces = mesher.sweep(
                &channel_def.profile,
                &channel_def.path,
                &mut current_ch.vertices,
                RegionId::new(0),
            );
            for face in faces { current_ch.faces.push(face); }
        }

        // Raw-concat for the channels-only visualization file.
        // Map each vertex through the VertexPool insert API, then remap each face.
        {
            let positions: Vec<_> = current_ch.vertices.positions().collect::<Vec<_>>();
            let mut id_map: Vec<cfd_mesh::domain::core::index::VertexId> = Vec::with_capacity(positions.len());
            for pos in &positions {
                id_map.push(all_channels.add_vertex_pos(nalgebra::Point3::from(pos.coords)));
            }
            for (_, face) in current_ch.faces.iter_enumerated() {
                let v0 = id_map[face.vertices[0].as_usize()];
                let v1 = id_map[face.vertices[1].as_usize()];
                let v2 = id_map[face.vertices[2].as_usize()];
                all_channels.add_face(v0, v1, v2);
            }
        }

        println!("  Subtracting from substrate...");
        match csg_boolean_indexed(BooleanOp::Difference, &final_solid, &current_ch) {
            Ok(m) => final_solid = m,
            Err(e) => {
                println!("  Warning: Channel subtraction failed: {}", e);
            }
        }
    }

    println!(
        "Final subtracted solid has {} vertices and {} faces.",
        final_solid.vertices.len(),
        final_solid.faces.len()
    );

    // Export the 3D mesh as STL
    let out_dir = manifest_dir.join("outputs/schematic_to_3d");
    std::fs::create_dir_all(&out_dir)?;
    
    let solid_path = out_dir.join("mirrored_bifurcation_solid.stl");
    let mut solid_file = File::create(&solid_path)?;
    write_stl_binary(&mut solid_file, &final_solid)?;
    
    let channels_path = out_dir.join("mirrored_bifurcation_channels.stl");
    let mut channels_file = File::create(&channels_path)?;
    write_stl_binary(&mut channels_file, &all_channels)?;

    println!("✅ Successfully exported Solid 3D Mesh to: {}", solid_path.display());
    println!("✅ Successfully exported Channels 3D Mesh to: {}", channels_path.display());

    // Copy the schematics over to bundle them
    let input_dir = manifest_dir.join("../cfd-schematics/outputs/mirrored_bifurcation");
    let json_src = input_dir.join("schematic.json");
    let svg_src = input_dir.join("schematic.svg");
    let png_src = input_dir.join("schematic.png");
    
    if json_src.exists() {
        std::fs::copy(&json_src, out_dir.join("schematic.json"))?;
        println!("✅ Copied schematic.json to 3D output dir");
    }
    if svg_src.exists() {
        std::fs::copy(&svg_src, out_dir.join("schematic.svg"))?;
        println!("✅ Copied schematic.svg to 3D output dir");
    }
    if png_src.exists() {
        std::fs::copy(&png_src, out_dir.join("schematic.png"))?;
        println!("✅ Copied schematic.png to 3D output dir");
    }

    Ok(())
}
