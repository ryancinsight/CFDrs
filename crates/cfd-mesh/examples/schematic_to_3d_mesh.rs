//! JSON schematic → watertight 3-D mesh (STL + OpenFOAM).
//!
//! Reads any `InterchangeChannelSystem` JSON produced by `cfd-schematics`,
//! sweeps each channel centerline into a closed tube mesh, subtracts every
//! tube from a bounding substrate via CSG, and writes:
//!
//! - `*_solid.stl`    — chip body (substrate − channel voids) for manufacturing
//! - `*_channels.stl` — channel void surfaces for inspection
//! - `constant/polyMesh/` — OpenFOAM surface mesh for snappyHexMesh / CFD
//!
//! ## Usage
//!
//! ```sh
//! # Default: uses mirrored_bifurcation schematic
//! cargo run -p cfd-mesh --example schematic_to_3d_mesh
//!
//! # Custom JSON path:
//! cargo run -p cfd-mesh --example schematic_to_3d_mesh -- path/to/schematic.json
//! ```
//!
//! If no JSON path is given, the example falls back to the mirrored_bifurcation
//! output produced by `cfd-schematics`:
//! ```sh
//! cargo run -p cfd-schematics --example mirrored_bifurcation
//! ```

use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

use cfd_mesh::application::channel::sweep::SweepMesher;
use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::domain::core::index::RegionId;
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::{Cube, PrimitiveMesh};
use cfd_mesh::domain::mesh::IndexedMesh;
use cfd_mesh::domain::topology::halfedge::PatchType;
use cfd_mesh::infrastructure::io::openfoam::write_openfoam_polymesh;
use cfd_mesh::infrastructure::io::scheme;
use cfd_mesh::infrastructure::io::stl::write_stl_binary;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════╗");
    println!("║  JSON Schematic → 3D Mesh (STL+OpenFOAM) ║");
    println!("╚══════════════════════════════════════════╝");

    // ── Resolve JSON path ─────────────────────────────────────────────────────
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            manifest_dir.join("../cfd-schematics/outputs/mirrored_bifurcation/schematic.json")
        });

    if !json_path.exists() {
        eprintln!("❌ Missing input file: {}", json_path.display());
        eprintln!("Please run the cfd-schematics mirrored_bifurcation example first:");
        eprintln!("  cargo run --example mirrored_bifurcation -p cfd-schematics");
        return Ok(());
    }

    let design_name = json_path
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|n| n.to_str())
        .unwrap_or("schematic");

    println!("📄 Loading JSON: {}", json_path.display());

    let mut file = File::open(&json_path)?;
    let mut json_str = String::new();
    file.read_to_string(&mut json_str)?;

    let interchange: cfd_schematics::geometry::InterchangeChannelSystem =
        serde_json::from_str(&json_str)?;

    // ── Convert to 3-D Schematic ──────────────────────────────────────────────
    let substrate_height = 5.0_f64; // mm
    let segments = 32;

    let schematic3d = scheme::from_interchange(
        &interchange,
        substrate_height as cfd_mesh::domain::core::scalar::Real,
        segments,
    )?;

    println!(
        "🔧 {} channels parsed (substrate: {}×{}×{} mm)",
        schematic3d.channels.len(),
        interchange.box_dims_mm.0,
        interchange.box_dims_mm.1,
        substrate_height,
    );

    // ── Build substrate block ─────────────────────────────────────────────────
    let (bw, bd) = interchange.box_dims_mm;
    let bw = bw as f64;
    let bd = bd as f64;
    let half_h = substrate_height / 2.0;

    let substrate = Cube {
        origin: Point3r::new(-5.0, -5.0, -half_h),
        width: (bw + 10.0) as _,
        height: (bd + 10.0) as _,
        depth: substrate_height as _,
    }
    .build()?;

    // ── Sweep channels + CSG subtract ────────────────────────────────────────
    let mesher = SweepMesher::new();
    let mut final_solid = substrate;
    let mut all_channels = IndexedMesh::new();

    for channel_def in &schematic3d.channels {
        println!("  → sweeping channel {} …", channel_def.id);

        let mut current_ch = IndexedMesh::new();
        if let Some(scales) = &channel_def.width_scales {
            let faces = mesher.sweep_variable(
                &channel_def.profile,
                &channel_def.path,
                scales,
                &mut current_ch.vertices,
                RegionId::new(0),
            );
            for face in faces {
                current_ch.faces.push(face);
            }
        } else {
            let faces = mesher.sweep(
                &channel_def.profile,
                &channel_def.path,
                &mut current_ch.vertices,
                RegionId::new(0),
            );
            for face in faces {
                current_ch.faces.push(face);
            }
        }

        // Accumulate channels-only mesh (for inspection STL)
        {
            let positions: Vec<_> = current_ch.vertices.positions().collect::<Vec<_>>();
            let mut id_map = Vec::with_capacity(positions.len());
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

        match csg_boolean_indexed(BooleanOp::Difference, &final_solid, &current_ch) {
            Ok(m) => final_solid = m,
            Err(e) => eprintln!("  ⚠  channel {} CSG failed: {}", channel_def.id, e),
        }
    }

    println!(
        "✅ Solid: {} vertices / {} faces",
        final_solid.vertices.len(),
        final_solid.faces.len()
    );

    // ── Output directories ────────────────────────────────────────────────────
    let out_dir = manifest_dir.join(format!("outputs/schematic_to_3d/{design_name}"));
    std::fs::create_dir_all(&out_dir)?;

    // ── STL export ────────────────────────────────────────────────────────────
    let solid_path = out_dir.join(format!("{design_name}_solid.stl"));
    write_stl_binary(&mut File::create(&solid_path)?, &final_solid)?;
    println!("📦 Solid STL  → {}", solid_path.display());

    let channels_path = out_dir.join(format!("{design_name}_channels.stl"));
    write_stl_binary(&mut File::create(&channels_path)?, &all_channels)?;
    println!("📦 Channels STL → {}", channels_path.display());

    // ── OpenFOAM export (solid mesh = chip body for snappyHexMesh) ───────────
    // All faces in the raw-CSG solid are unlabeled wall faces (region 0).
    // Export the chip body as a wall-only surface — snappyHexMesh will use it
    // as the geometry input and label its own patches from castellated-mesh BCs.
    let of_dir = out_dir.join("constant/polyMesh");
    write_openfoam_polymesh(
        &final_solid,
        &of_dir,
        &[(RegionId::new(0), "walls", PatchType::Wall)],
    )?;
    println!("🌊 OpenFOAM   → {}/", of_dir.display());

    // ── Copy originating schematics for reference ─────────────────────────────
    if let Some(src_dir) = json_path.parent() {
        for ext in &["json", "svg", "png"] {
            let src = src_dir.join(format!("schematic.{ext}"));
            if src.exists() {
                std::fs::copy(&src, out_dir.join(format!("schematic.{ext}")))?;
                println!("📎 Copied schematic.{ext}");
            }
        }
    }

    println!("\n✅ All outputs written to: {}", out_dir.display());
    Ok(())
}
