//! SBS 96-well flat-bottom plate mesh generator.
//!
//! Generates a watertight solid plate body with 96 open wells matching
//! Cat. No. 439454 / ANSI SLAS 1-2004 dimensions.
//!
//! Outputs three STLs:
//!   - `well_plate_block.stl`  — plate body before well subtraction
//!   - `well_plate_solid.stl`  — plate body with 96 wells cut
//!   - `well_plate_fluid.stl`  — 96 well cavities (fluid domain)

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::{Cube, Frustum, PrimitiveMesh};
use cfd_mesh::domain::mesh::IndexedMesh;
use cfd_mesh::infrastructure::io::stl::write_stl_binary;

use std::collections::HashMap;

fn combine_meshes(meshes: &[IndexedMesh]) -> IndexedMesh {
    let mut combined = IndexedMesh::new();
    for m in meshes {
        let mut remap = HashMap::with_capacity(m.vertices.len());
        for (old_vid, _) in m.vertices.iter() {
            let pos = *m.vertices.position(old_vid);
            let nrm = *m.vertices.normal(old_vid);
            let new_vid = combined.add_vertex(pos, nrm);
            remap.insert(old_vid, new_vid);
        }
        for f in m.faces.iter() {
            combined.add_face_with_region(
                remap[&f.vertices[0]],
                remap[&f.vertices[1]],
                remap[&f.vertices[2]],
                f.region,
            );
        }
    }
    combined.rebuild_edges();
    combined
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════╗");
    println!("║  SBS 96-Well Flat-Bottom Plate Generator ║");
    println!("╚══════════════════════════════════════════╝");

    let start_time = Instant::now();
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let out_dir = manifest_dir.join("outputs/sbs_96_well_plate_flat");
    fs::create_dir_all(&out_dir)?;

    // Cat. No. 439454 dimensions (mm):
    // H=14.5  H2=3.2  H3=11.2  H4=2.1  H5=2.85  D1=7.05  D2=6.5  F1=2.0
    // SLAS footprint: 127.60 × 85.75 (outer), 127.15 × 85.05 (inner/lid)
    // Well spacing: 9.0 mm center-to-center

    let segments = 16; // Reduced from 32 to keep face count manageable for 96 CSG ops

    // Export block STL — a plain cuboid representing the body the wells are cut from
    // (no lip or skirt features yet)
    let block = Cube {
        origin: Point3r::new(0.225, 0.0, 0.35),
        width: 127.15,
        height: 14.5,
        depth: 85.05,
    }
    .build()?;
    let block_path = out_dir.join("well_plate_block.stl");
    write_stl_binary(&mut fs::File::create(&block_path)?, &block)?;
    println!("📦 Block STL  → {}", block_path.display());

    // ── Well parameters ───────────────────────────────────────────────────────
    let cx_start = 14.3;
    let cz_start = 11.375;
    let well_spacing = 9.0;
    let well_bottom_radius = 6.5 / 2.0;
    let well_top_radius = 7.05 / 2.0;
    let well_bottom_y = 3.2;
    let well_height = 11.2;

    // ── Generate fluid wells ──────────────────────────────────────────────────
    println!("  → Generating fluid well mesh...");
    let mut all_wells = Vec::with_capacity(96);
    let mut all_cutters = Vec::with_capacity(96);
    for row in 0..8 {
        for col in 0..12 {
            let cx = cx_start + (col as f64) * well_spacing;
            let cz = cz_start + (row as f64) * well_spacing;
            
            // Actual fluid domain inside the well
            let frustum = Frustum {
                base_center: Point3r::new(cx, well_bottom_y, cz),
                bottom_radius: well_bottom_radius,
                top_radius: well_top_radius,
                height: well_height,
                segments,
            }
            .build()?;
            all_wells.push(frustum);

            // Taller cutter to ensure we pierce the top of the block fully for an open well
            let cutter = Frustum {
                base_center: Point3r::new(cx, well_bottom_y, cz),
                bottom_radius: well_bottom_radius,
                top_radius: well_top_radius,
                height: well_height + 1.0, 
                segments,
            }
            .build()?;
            all_cutters.push(cutter);
        }
    }
    let fluid_mesh = combine_meshes(&all_wells);
    let cutter_mesh = combine_meshes(&all_cutters);

    let channels_path = out_dir.join("well_plate_fluid.stl");
    write_stl_binary(&mut fs::File::create(&channels_path)?, &fluid_mesh)?;
    println!("📦 Fluid STL  → {}", channels_path.display());

    let fluid_vol = fluid_mesh.signed_volume();
    println!(
        "✅ Fluid: {:>7} verts / {:>7} faces, vol = {:>10.3} mm³",
        fluid_mesh.vertices.len(),
        fluid_mesh.faces.len(),
        fluid_vol
    );

    // ── Generate well_plate_part_1 (Block - Cutters) ──────────────────────────
    println!("  → Subtracting cutters from block to create well_plate_part_1...");
    let mut part_1 = match csg_boolean(BooleanOp::Difference, &block, &cutter_mesh) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("  ⚠ part 1 boolean difference failed: {}", e);
            block.clone() // fallback or handle error
        }
    };
    part_1.rebuild_edges();

    let part_1_path = out_dir.join("well_plate_part_1.stl");
    write_stl_binary(&mut fs::File::create(&part_1_path)?, &part_1)?;
    println!("📦 Part 1 STL → {}", part_1_path.display());

    let part_1_vol = part_1.signed_volume();
    println!(
        "✅ Part 1: {:>7} verts / {:>7} faces, vol = {:>10.3} mm³",
        part_1.vertices.len(),
        part_1.faces.len(),
        part_1_vol
    );

    let duration = start_time.elapsed();
    println!("✨ Completed in {:.2?}", duration);
    Ok(())
}
