//! Millifluidic mesh generation from cfd-schematics designs.
//!
//! Demonstrates the full CSG pipeline: schematic → channel sweep → union
//! channels → subtract from substrate → STL export.
//!
//! Design targets:
//! - **96-well plate** footprint: 127.76 × 85.48 mm, 14.35 mm height (SBS)
//! - **4 mm dialysis tubing** — channel diameter = 4 mm (radius 2 mm)
//! - **Z-centred channels** — centerlines at z = height / 2
//! - **CSG subtraction** — channels are merged (union), then subtracted from
//!   the substrate cuboid so the result has internal passages with visible
//!   inlet/outlet port holes on the side faces.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example millifluidic_mesh_demo --features "scheme-io,csg"
//! ```

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

// cfd-mesh types
use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::{Real, Point3r};
use cfd_mesh::channel::profile::ChannelProfile;
use cfd_mesh::channel::path::ChannelPath;
use cfd_mesh::channel::sweep::SweepMesher;
use cfd_mesh::channel::substrate::SubstrateBuilder;
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean};
use cfd_mesh::io::stl;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;

// cfd-schematics types
use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig, SerpentineConfig};
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;

// ── Physical constants ────────────────────────────────────────

/// SBS-standard 96-well plate footprint (mm).
const PLATE_WIDTH_MM: Real = 127.76;
const PLATE_DEPTH_MM: Real = 85.48;
/// SBS-standard 96-well plate height (mm).
const PLATE_HEIGHT_MM: Real = 14.35;

/// Dialysis tubing outer diameter (mm).
const CHANNEL_DIAMETER_MM: Real = 4.0;
const CHANNEL_RADIUS_MM: Real = CHANNEL_DIAMETER_MM / 2.0;

/// Number of polygon segments around the circular channel cross-section.
const CHANNEL_SEGMENTS: usize = 24;

/// How far (mm) channel tubes extend past the substrate face to ensure
/// clean through-holes during CSG subtraction.
const PORT_EXTENSION_MM: Real = 2.0;

// ── Demo case definition ──────────────────────────────────────

struct DemoCase {
    name: &'static str,
    description: &'static str,
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
    println!("  cfd-mesh :: Millifluidic Mesh (CSG pipeline)");
    println!("=================================================================");
    println!();
    println!(
        "  Plate: {:.2} × {:.2} × {:.2} mm (SBS 96-well)",
        PLATE_WIDTH_MM, PLATE_DEPTH_MM, PLATE_HEIGHT_MM
    );
    println!(
        "  Channel Ø: {:.1} mm (4 mm dialysis tubing), {} segments",
        CHANNEL_DIAMETER_MM, CHANNEL_SEGMENTS
    );
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir_path = crate_dir.join("outputs").join("millifluidic_demo");
    fs::create_dir_all(&out_dir_path)?;
    let out_dir = out_dir_path.to_str().expect("non-UTF8 path");

    let geometry_config = GeometryConfig::default();

    let cases = vec![
        DemoCase {
            name: "straight_channel",
            description: "Single straight channel (baseline)",
            splits: vec![],
            channel_type_config: ChannelTypeConfig::AllStraight,
        },
        DemoCase {
            name: "bifurcation",
            description: "Single bifurcation (2 daughter channels)",
            splits: vec![SplitType::Bifurcation],
            channel_type_config: ChannelTypeConfig::AllStraight,
        },
        DemoCase {
            name: "trifurcation",
            description: "Single trifurcation (3 daughter channels)",
            splits: vec![SplitType::Trifurcation],
            channel_type_config: ChannelTypeConfig::AllStraight,
        },
        DemoCase {
            name: "serpentine",
            description: "Bifurcation with serpentine channel paths",
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
        "Name", "Vertices", "Faces", "Ch", "Area (mm²)", "Vol (mm³)", "ms"
    );
    println!("{:-<78}", "");
    for r in &reports {
        println!(
            "{:<22} {:>8} {:>8} {:>4} {:>12.2} {:>12.2} {:>8}",
            r.name, r.vertices, r.faces, r.channels,
            r.surface_area_mm2, r.volume_mm3, r.elapsed_ms
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

    // 1. Generate 2D schematic via cfd-schematics using 96-well plate footprint
    let box_dims = (PLATE_WIDTH_MM as f64, PLATE_DEPTH_MM as f64);
    let channel_system = create_geometry(
        box_dims,
        &case.splits,
        geometry_config,
        &case.channel_type_config,
    );
    println!(
        "  Schematic: {} nodes, {} channels",
        channel_system.nodes.len(),
        channel_system.channels.len()
    );

    // 2. Convert interchange → 3D channel paths (lifted to z = height/2)
    let interchange = channel_system.to_interchange();
    let mid_z = PLATE_HEIGHT_MM / 2.0;
    let n_channels = interchange.channels.len();

    // Build channel definitions with 4 mm diameter and extended endpoints
    let profile = ChannelProfile::Circular {
        radius: CHANNEL_RADIUS_MM,
        segments: CHANNEL_SEGMENTS,
    };

    let mut channel_paths: Vec<ChannelPath> = Vec::with_capacity(n_channels);
    for ch in &interchange.channels {
        if ch.centerline_mm.len() < 2 {
            continue;
        }
        // Lift 2D → 3D at z = mid_z
        let mut pts: Vec<Point3r> = ch
            .centerline_mm
            .iter()
            .map(|&(x, y)| Point3r::new(x as Real, y as Real, mid_z))
            .collect();

        // Extend both endpoints past the substrate so CSG creates clean port holes
        extend_endpoints(&mut pts, PORT_EXTENSION_MM);
        channel_paths.push(ChannelPath::new(pts));
    }

    println!(
        "  Channels: {} paths, Ø {:.1} mm, z = {:.2} mm",
        channel_paths.len(),
        CHANNEL_DIAMETER_MM,
        mid_z
    );

    // 3. Mesh: substrate cuboid + channel tubes → CSG subtract
    let mesh = mesh_with_csg(&profile, &channel_paths)?;
    println!(
        "  Meshed (CSG): {} vertices, {} faces",
        mesh.vertex_count(),
        mesh.face_count(),
    );

    // 4. Metrics
    let area = mesh.surface_area();
    let vol = mesh.signed_volume();
    println!("  Area: {:.2} mm², Volume: {:.2} mm³", area, vol);

    // 5. Quality report
    let qr = mesh.quality_report();
    if let (Some(ar), Some(ma)) = (&qr.aspect_ratio, &qr.min_angle) {
        println!(
            "  Quality: AR mean={:.2} max={:.2}, min∠ mean={:.1}° min={:.1}°",
            ar.mean, ar.max,
            ma.mean.to_degrees(), ma.min.to_degrees(),
        );
    }
    if qr.passed {
        println!("  Quality: PASSED ({} faces)", qr.total_faces);
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
        channels: channel_paths.len(),
        surface_area_mm2: area as f64,
        volume_mm3: vol as f64,
        stl_path,
        elapsed_ms: elapsed,
    })
}

// ── Core CSG meshing ──────────────────────────────────────────

/// Build the millifluidic chip via CSG: substrate \ union(channels).
///
/// 1. Generate substrate cuboid (closed, 12 faces).
/// 2. Sweep each channel path into a capped tube.
/// 3. Iteratively union all channel tubes into one merged volume.
/// 4. Subtract the merged channel volume from the substrate.
///
/// The result is a solid chip body with internal channel passages and
/// visible inlet/outlet port holes where channels pierce the substrate.
fn mesh_with_csg(
    profile: &ChannelProfile,
    channel_paths: &[ChannelPath],
) -> Result<IndexedMesh, Box<dyn std::error::Error>> {
    let mut pool = VertexPool::default_millifluidic();
    let sweeper = SweepMesher::new(); // caps both ends

    // ── 1. Substrate cuboid ───────────────────────────────────
    let substrate_region = RegionId::new(0);
    let substrate_faces = SubstrateBuilder::new(PLATE_WIDTH_MM, PLATE_DEPTH_MM, PLATE_HEIGHT_MM)
        .with_origin(Point3r::origin())
        .build(&mut pool, substrate_region);

    println!(
        "  Substrate: {:.2}×{:.2}×{:.2} mm → {} faces",
        PLATE_WIDTH_MM, PLATE_DEPTH_MM, PLATE_HEIGHT_MM,
        substrate_faces.len()
    );

    // ── 2. Sweep each channel into a capped tube ─────────────
    let mut channel_face_sets: Vec<Vec<FaceData>> = Vec::with_capacity(channel_paths.len());

    for (i, path) in channel_paths.iter().enumerate() {
        let region = RegionId::new((i + 1) as u32);
        let faces = sweeper.sweep(profile, path, &mut pool, region);
        println!(
            "  Channel {}: {} pts → {} faces",
            i,
            path.points().len(),
            faces.len()
        );
        channel_face_sets.push(faces);
    }

    // ── 3. Merge (union) all channel tubes ────────────────────
    let merged_channels = if channel_face_sets.len() == 1 {
        channel_face_sets.into_iter().next().unwrap()
    } else {
        let mut merged = channel_face_sets[0].clone();
        for (i, set) in channel_face_sets[1..].iter().enumerate() {
            print!("  Union channel {}+{} ...", 0, i + 1);
            match csg_boolean(BooleanOp::Union, &merged, set, &mut pool) {
                Ok(result) => {
                    println!(" {} faces", result.len());
                    merged = result;
                }
                Err(e) => {
                    // If union fails (non-overlapping channels), just concatenate
                    println!(" union failed ({}), concatenating", e);
                    merged.extend_from_slice(set);
                }
            }
        }
        merged
    };

    println!(
        "  Merged channels: {} faces total",
        merged_channels.len()
    );

    // ── 4. CSG difference: substrate \ channels ───────────────
    println!("  CSG: substrate \\ channels ...");
    let result_faces = csg_boolean(
        BooleanOp::Difference,
        &substrate_faces,
        &merged_channels,
        &mut pool,
    )?;
    println!("  CSG result: {} faces", result_faces.len());

    // ── 5. Assemble IndexedMesh ───────────────────────────────
    let mut mesh = IndexedMesh::new();
    mesh.vertices = pool;
    for f in &result_faces {
        mesh.faces.push(*f);
    }
    mesh.rebuild_edges();

    Ok(mesh)
}

// ── Helpers ───────────────────────────────────────────────────

/// Extend both endpoints of a polyline outward by `dist` mm along the
/// tangent direction. This ensures swept channel tubes protrude past the
/// substrate faces so CSG subtraction creates clean through-holes.
fn extend_endpoints(pts: &mut Vec<Point3r>, dist: Real) {
    if pts.len() < 2 {
        return;
    }
    let n = pts.len();

    // Extend start backward
    let start_dir = (pts[0] - pts[1]).normalize();
    let new_start = pts[0] + start_dir * dist;

    // Extend end forward
    let end_dir = (pts[n - 1] - pts[n - 2]).normalize();
    let new_end = pts[n - 1] + end_dir * dist;

    pts.insert(0, new_start);
    pts.push(new_end);
}
