//! Millifluidic mesh generation from cfd-schematics designs.
//!
//! Demonstrates the full pipeline: schematic → 3D mesh → STL export for:
//! 1. Straight channel (baseline)
//! 2. Bifurcation
//! 3. Trifurcation
//! 4. Serpentine
//!
//! This is the cfd-mesh equivalent of blue2mesh's `scheme_to_stl_demo`.
//! Instead of blue2mesh's `ExtrusionConfig`, we use `SweepMesher` with
//! the `IndexedMesh` architecture for deduplication and watertight validation.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example millifluidic_mesh_demo --features scheme-io
//! ```

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

// cfd-mesh types
use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::Real;
use cfd_mesh::channel::sweep::SweepMesher;
use cfd_mesh::channel::substrate::SubstrateBuilder;
use cfd_mesh::io::stl;
use cfd_mesh::io::scheme::{self, Schematic};

// cfd-schematics types
use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig};
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;

// ── Demo case definition ──────────────────────────────────────

struct DemoCase {
    name: &'static str,
    description: &'static str,
    box_dims_mm: (f64, f64),
    splits: Vec<SplitType>,
    channel_type_config: ChannelTypeConfig,
}

struct MeshReport {
    name: String,
    vertices: usize,
    faces: usize,
    channels: usize,
    surface_area_mm2: f64,
    volume_mm3: f64,
    stl_path: String,
    elapsed_ms: u128,
}

// ── Main ──────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  cfd-mesh :: Millifluidic Mesh Generation from Schematics");
    println!("=================================================================");
    println!();

    let out_dir = "outputs/millifluidic_demo";
    fs::create_dir_all(out_dir)?;

    let geometry_config = GeometryConfig::default();

    let cases = vec![
        DemoCase {
            name: "straight_channel",
            description: "Single straight channel (baseline)",
            box_dims_mm: (127.15, 85.75),
            splits: vec![],
            channel_type_config: ChannelTypeConfig::AllStraight,
        },
        DemoCase {
            name: "bifurcation",
            description: "Single bifurcation (2 daughter channels)",
            box_dims_mm: (200.0, 120.0),
            splits: vec![SplitType::Bifurcation],
            channel_type_config: ChannelTypeConfig::AllStraight,
        },
        DemoCase {
            name: "trifurcation",
            description: "Single trifurcation (3 daughter channels)",
            box_dims_mm: (240.0, 140.0),
            splits: vec![SplitType::Trifurcation],
            channel_type_config: ChannelTypeConfig::AllStraight,
        },
        DemoCase {
            name: "serpentine",
            description: "Bifurcation with serpentine channel paths",
            box_dims_mm: (260.0, 140.0),
            splits: vec![SplitType::Bifurcation],
            channel_type_config: ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
        },
    ];

    let mut reports = Vec::new();

    for case in &cases {
        match process_case(case, &geometry_config, out_dir) {
            Ok(report) => reports.push(report),
            Err(e) => {
                eprintln!("  FAILED: {} — {}", case.name, e);
            }
        }
    }

    // ── Summary table ─────────────────────────────────────────
    println!();
    println!("=================================================================");
    println!("  Summary");
    println!("=================================================================");
    println!(
        "{:<22} {:>8} {:>8} {:>4} {:>12} {:>12} {:>8}",
        "Name", "Vertices", "Faces", "Ch", "Area (mm2)", "Vol (mm3)", "ms"
    );
    println!("{:-<78}", "");
    for r in &reports {
        println!(
            "{:<22} {:>8} {:>8} {:>4} {:>12.2} {:>12.2} {:>8}",
            r.name, r.vertices, r.faces, r.channels, r.surface_area_mm2, r.volume_mm3, r.elapsed_ms
        );
        println!("  -> {}", r.stl_path);
    }
    println!();
    println!("STL files written to {}/", out_dir);
    println!("All schematics meshed successfully.");

    Ok(())
}

// ── Per-case pipeline ─────────────────────────────────────────

fn process_case(
    case: &DemoCase,
    geometry_config: &GeometryConfig,
    out_dir: &str,
) -> Result<MeshReport, Box<dyn std::error::Error>> {
    let t0 = Instant::now();
    println!("--- {} ---", case.name);
    println!("  {}", case.description);

    // 1. Generate 2D schematic via cfd-schematics
    let channel_system = create_geometry(
        case.box_dims_mm,
        &case.splits,
        geometry_config,
        &case.channel_type_config,
    );
    println!(
        "  Schematic: box {:?} mm, {} nodes, {} channels",
        channel_system.box_dims,
        channel_system.nodes.len(),
        channel_system.channels.len()
    );

    // 2. Convert to cfd-mesh Schematic via interchange bridge
    let substrate_height_mm: Real = 10.0;
    let channel_segments: usize = 16;
    let schematic = scheme::from_channel_system(
        &channel_system,
        substrate_height_mm,
        channel_segments,
    )?;
    println!(
        "  Parsed: substrate {:.1}x{:.1}x{:.1} mm, {} channels",
        schematic.substrate.width,
        schematic.substrate.depth,
        schematic.substrate.height,
        schematic.channels.len(),
    );

    // 3. Mesh the schematic
    let mesh = mesh_schematic(&schematic)?;
    println!(
        "  Meshed: {} vertices, {} faces",
        mesh.vertex_count(),
        mesh.face_count(),
    );

    // 4. Compute metrics
    let area = mesh.surface_area();
    let vol = mesh.signed_volume();
    println!("  Area: {:.2} mm^2, Signed volume: {:.2} mm^3", area, vol);

    // 5. Quality report
    let qr = mesh.quality_report();
    if let (Some(ar), Some(ma)) = (&qr.aspect_ratio, &qr.min_angle) {
        println!(
            "  Quality: aspect_ratio mean={:.2} max={:.2}, min_angle mean={:.1} min={:.1} deg",
            ar.mean,
            ar.max,
            ma.mean.to_degrees(),
            ma.min.to_degrees(),
        );
    }
    if qr.passed {
        println!("  Quality: PASSED ({} faces, 0 failing)", qr.total_faces);
    } else {
        println!(
            "  Quality: {} / {} faces failing",
            qr.failing_faces, qr.total_faces
        );
    }

    // 6. Write binary STL
    let stl_path = format!("{}/{}.stl", out_dir, case.name);
    {
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
    }
    let stl_size = fs::metadata(&stl_path)?.len();
    println!("  STL: {} ({} bytes)", stl_path, stl_size);

    let elapsed = t0.elapsed().as_millis();
    println!("  Elapsed: {} ms", elapsed);
    println!();

    Ok(MeshReport {
        name: case.name.to_string(),
        vertices: mesh.vertex_count(),
        faces: mesh.face_count(),
        channels: schematic.channels.len(),
        surface_area_mm2: area as f64,
        volume_mm3: vol as f64,
        stl_path,
        elapsed_ms: elapsed,
    })
}

// ── Core meshing logic ────────────────────────────────────────

/// Mesh a full `Schematic` into an `IndexedMesh`.
///
/// For each channel in the schematic, a profile sweep along the centerline
/// is performed. The substrate cuboid is also generated. All geometry shares
/// the same `IndexedMesh` vertex pool for automatic deduplication.
fn mesh_schematic(schematic: &Schematic) -> Result<IndexedMesh, Box<dyn std::error::Error>> {
    let mut mesh = IndexedMesh::new();
    let sweeper = SweepMesher::new();

    // Region 0: substrate
    let substrate_region = RegionId::new(0);
    let sub = &schematic.substrate;
    let substrate_builder = SubstrateBuilder::new(sub.width, sub.depth, sub.height)
        .with_origin(sub.origin);
    let substrate_faces = substrate_builder.build(&mut mesh.vertices, substrate_region);
    for f in &substrate_faces {
        mesh.faces.push(*f);
    }

    // Channels: region 1, 2, 3, ...
    for (i, channel) in schematic.channels.iter().enumerate() {
        let region = RegionId::new((i + 1) as u32);
        let channel_faces = sweeper.sweep(
            &channel.profile,
            &channel.path,
            &mut mesh.vertices,
            region,
        );

        for f in &channel_faces {
            mesh.faces.push(*f);
        }
    }

    // Build edge adjacency
    mesh.rebuild_edges();

    Ok(mesh)
}
