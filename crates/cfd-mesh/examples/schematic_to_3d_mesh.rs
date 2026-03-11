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
use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
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

    let interchange: cfd_schematics::NetworkBlueprint = serde_json::from_str(&json_str)?;

    // ── Convert to 3-D Schematic ──────────────────────────────────────────────
    let substrate_height = 5.0_f64; // mm
    let segments = 32;

    let schematic3d = scheme::from_blueprint(
        &interchange,
        substrate_height as cfd_mesh::domain::core::scalar::Real,
        segments,
    )?;

    println!(
        "🔧 {} channels parsed (substrate: {}×{}×{} mm)",
        schematic3d.channels.len(),
        interchange.box_dims.0,
        interchange.box_dims.1,
        substrate_height,
    );

    // ── Build substrate block ─────────────────────────────────────────────────
    let (bw, bd) = interchange.box_dims;
    let bw = bw as f64;
    let bd = bd as f64;
    let _half_h = substrate_height / 2.0;

    // Offset slightly so that the X=0 and X=100 boundaries of the substrate
    // don't perfectly hit the sweep segment boundaries causing CSG degeneracies
    // when subtracting to form open ports.
    let substrate = Cube {
        origin: Point3r::new(-0.1, 0.0, 0.0),
        width: (bw + 0.2) as _,
        height: bd as _,
        depth: substrate_height as _,
    }
    .build()?;

    // ── Topology-aware sweep: suppress caps at T/Y-junction nodes ────────────
    //
    // When two or more channels share an endpoint (e.g., the split node at X=25),
    // a capped sweep leaves a flat circular disk coplanar with the junction.
    // CSG Union at these coplanar seams is numerically ambiguous (GWN ≈ 0.5),
    // producing phantom "fan" triangles. Solution: detect shared endpoints and
    // suppress the cap on those ends — open-ended tubes T-intersect cleanly.

    // Collect all path endpoints (first and last point of each channel).
    let tol = 1e-4_f64;
    let endpoint_positions: Vec<(Point3r, Point3r)> = schematic3d
        .channels
        .iter()
        .map(|ch| {
            let pts = ch.path.points();
            let first = *pts.first().unwrap();
            let last = *pts.last().unwrap();
            (first, last)
        })
        .collect();

    // A position is a "junction" if it appears as an endpoint more than once.
    let is_junction = |pt: &Point3r| -> bool {
        let count = endpoint_positions
            .iter()
            .flat_map(|(a, b)| [a, b])
            .filter(|q| (q.coords - pt.coords).norm() < tol)
            .count();
        count > 1
    };

    let mut final_solid = substrate;
    let mut all_channels: Option<IndexedMesh> = None;

    for channel_def in &schematic3d.channels {
        println!("  → sweeping channel {} …", channel_def.id);

        let pts = channel_def.path.points();
        let first = *pts.first().unwrap();
        let last = *pts.last().unwrap();

        // Only cap ends that are open ports (not shared with another channel).
        let cap_start = !is_junction(&first);
        let cap_end = !is_junction(&last);

        let sweep = SweepMesher { cap_start, cap_end };

        let mut current_ch = IndexedMesh::new();
        if let Some(scales) = &channel_def.width_scales {
            let faces = sweep.sweep_variable(
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
            let faces = sweep.sweep(
                &channel_def.profile,
                &channel_def.path,
                &mut current_ch.vertices,
                RegionId::new(0),
            );
            for face in faces {
                current_ch.faces.push(face);
            }
        }

        current_ch.rebuild_edges();

        // Union into accumulated fluid mesh.
        if let Some(existing) = all_channels.take() {
            match csg_boolean(BooleanOp::Union, &existing, &current_ch) {
                Ok(m) => all_channels = Some(m),
                Err(_) => match csg_boolean(BooleanOp::Union, &existing, &current_ch) {
                    Ok(m) => all_channels = Some(m),
                    Err(e) => {
                        eprintln!("  ⚠ union for channel {} failed: {}", channel_def.id, e);
                        all_channels = Some(existing);
                    }
                },
            }
        } else {
            all_channels = Some(current_ch);
        }
    }

    let mut all_channels = all_channels.expect("Must have at least one channel");

    // Repair orientation and clean up floating artifact faces from overlapping Union
    println!("  → post-processing fluid mesh repairs …");
    println!(
        "      volume before repair: {}",
        all_channels.signed_volume()
    );
    all_channels.orient_outward();
    all_channels.retain_largest_component();
    all_channels.rebuild_edges();
    println!(
        "      volume after repair: {}",
        all_channels.signed_volume()
    );

    // Single difference: substrate minus the consolidated watertight channel mesh.
    println!("  → applying boolean difference (substrate − channels)");
    match csg_boolean(BooleanOp::Difference, &final_solid, &all_channels) {
        Ok(mut m) => {
            // Post-difference repair: the Difference can leave phantom seam triangles
            // at T-junctions (flat caps of channels that intersect the new solid face).
            // orient_outward() + retain_largest_component() removes disconnected artifacts.
            m.orient_outward();
            m.retain_largest_component();
            m.rebuild_edges();
            final_solid = m;
        }
        Err(e) => eprintln!("  ⚠ final CSG difference failed: {}", e),
    }

    let solid_vol = final_solid.signed_volume();
    let channels_vol = all_channels.signed_volume();

    println!(
        "  → [DEBUG] final_solid AABB: {:?}",
        final_solid.bounding_box()
    );
    println!(
        "  → [DEBUG] all_channels AABB: {:?}",
        all_channels.bounding_box()
    );

    println!(
        "✅ Solid: {:>7} vertices / {:>7} faces, vol = {:>10.3} mm³",
        final_solid.vertices.len(),
        final_solid.faces.len(),
        solid_vol
    );
    println!(
        "✅ Chans: {:>7} vertices / {:>7} faces, vol = {:>10.3} mm³",
        all_channels.vertices.len(),
        all_channels.faces.len(),
        channels_vol
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
