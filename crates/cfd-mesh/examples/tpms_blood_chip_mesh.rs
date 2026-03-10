//! TPMS Millifluidic Chip — 3D Mesh Assembly
//!
//! Builds a complete millifluidic chip device containing a Gyroid TPMS lattice
//! inside a rectangular cavity with inlet/outlet port connections. The output
//! is a recognisable physical device — not just a bare mathematical surface.
//!
//! # Chip Design
//!
//! ```text
//!   ┌─────────────────────────────────────┐
//!   │         Top plate (lid)             │ Z = 10..12
//!   ├────┬───────────────────────┬────────┤
//!   │wall│  Gyroid TPMS lattice  │  wall  │ Z = 2..10
//!  ═╡port│  inside 20×20×8 mm   │  port╞═
//!   │    │       cavity          │        │
//!   ├────┴───────────────────────┴────────┤
//!   │        Bottom plate                 │ Z = 0..2
//!   └─────────────────────────────────────┘
//!     X=0                              X=24
//! ```
//!
//! # Outputs
//!
//! - `chip_assembly.stl` — open-top chip body with port holes + TPMS lattice
//! - `chip_lid.stl`      — top plate (load alongside assembly for closed chip)
//! - `tpms_lattice.stl`  — standalone Gyroid lattice (detail inspection)
//!
//! # Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example tpms_blood_chip_mesh --release
//! ```

use std::collections::HashMap;
use std::fs;
use std::io::BufWriter;
use std::path::Path;

use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::domain::core::index::VertexId;
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::cube::Cube;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::domain::geometry::tpms::{build_tpms_box, Gyroid, TpmsBoxParams};
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::IndexedMesh;

/// Copy all vertices and faces from `source` into `target`.
///
/// Vertex IDs are remapped via `VertexPool` deduplication; shared boundary
/// vertices are automatically welded.
fn merge_into(target: &mut IndexedMesh, source: &IndexedMesh) {
    let mut remap: HashMap<VertexId, VertexId> = HashMap::new();
    for (old_id, _) in source.vertices.iter() {
        let pos = *source.vertices.position(old_id);
        let nrm = *source.vertices.normal(old_id);
        let new_id = target.add_vertex(pos, nrm);
        remap.insert(old_id, new_id);
    }
    for face in source.faces.iter() {
        target.add_face(
            remap[&face.vertices[0]],
            remap[&face.vertices[1]],
            remap[&face.vertices[2]],
        );
    }
}

fn build_cube(ox: f64, oy: f64, oz: f64, w: f64, h: f64, d: f64) -> IndexedMesh {
    Cube {
        origin: Point3r::new(ox, oy, oz),
        width: w,
        height: h,
        depth: d,
    }
    .build()
    .expect("valid cube dimensions")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════╗");
    println!("║  TPMS Millifluidic Chip — 3D Mesh Assembly             ║");
    println!("║  cfd-schematics → cfd-mesh pipeline                    ║");
    println!("╚══════════════════════════════════════════════════════════╝");
    println!();

    // ── 1. Chip geometry ─────────────────────────────────────────────────────
    //
    // Compact demonstrator chip:
    //   Outer : 24 × 24 × 12 mm
    //   Cavity: 20 × 20 × 8 mm  (2 mm walls on all six sides)
    //   Ports : 2 × 2 mm rectangular channels cut through left/right walls
    //   Gyroid λ = 5 mm → 4 periods across cavity width, ~1.6 vertically
    let wall = 2.0_f64;
    let cavity_w = 20.0_f64; // X extent of cavity
    let cavity_h = 20.0_f64; // Y extent of cavity
    let cavity_z = 8.0_f64; // Z extent of cavity
    let outer_w = cavity_w + 2.0 * wall; // 24 mm
    let outer_h = cavity_h + 2.0 * wall; // 24 mm
    let chip_z = cavity_z + 2.0 * wall; // 12 mm

    let port_size = 2.0_f64; // port cross-section (square)
    let port_y = (outer_h - port_size) / 2.0; // centred in Y
    let port_z = (chip_z - port_size) / 2.0; // centred in Z

    println!("Chip geometry:");
    println!("  Outer     : {outer_w:.1} × {outer_h:.1} × {chip_z:.1} mm");
    println!("  Cavity    : {cavity_w:.1} × {cavity_h:.1} × {cavity_z:.1} mm");
    println!("  Wall      : {wall:.1} mm all sides");
    println!(
        "  Inlet port: {port_size:.0}×{port_size:.0} mm hole through left wall (X = 0..{wall:.0})"
    );
    println!(
        "  Outlet port: {port_size:.0}×{port_size:.0} mm hole through right wall (X = {:.0}..{outer_w:.0})",
        outer_w - wall
    );
    println!();

    // ── 2. Build chip body components ────────────────────────────────────────
    //
    // Open-top chip = bottom plate + 4 side walls.
    // Inlet/outlet are rectangular channels cut through the left/right walls
    // via CSG Boolean Difference — fluid enters/exits through the holes.
    println!("Building chip body (open-top with port holes) ...");

    // Bottom plate: full footprint, Z = 0..wall
    let bottom = build_cube(0.0, 0.0, 0.0, outer_w, outer_h, wall);

    // Left wall:  X = 0..wall, Y = 0..outer_h, Z = wall..wall+cavity_z
    let left_solid = build_cube(0.0, 0.0, wall, wall, outer_h, cavity_z);

    // Right wall: X = outer_w-wall..outer_w, same Y/Z as left
    let right_solid = build_cube(outer_w - wall, 0.0, wall, wall, outer_h, cavity_z);

    // Front wall: X = wall..outer_w-wall, Y = 0..wall, Z = wall..wall+cavity_z
    let front = build_cube(wall, 0.0, wall, cavity_w, wall, cavity_z);

    // Back wall:  X = wall..outer_w-wall, Y = outer_h-wall..outer_h
    let back = build_cube(wall, outer_h - wall, wall, cavity_w, wall, cavity_z);

    // Port cutter geometry — cubes that extend beyond wall faces for clean CSG cuts.
    // Inlet cutter passes through the left wall (X = -0.1 .. wall+0.1).
    let inlet_cutter = build_cube(-0.1, port_y, port_z, wall + 0.2, port_size, port_size);
    // Outlet cutter passes through the right wall.
    let outlet_cutter = build_cube(
        outer_w - wall - 0.1,
        port_y,
        port_z,
        wall + 0.2,
        port_size,
        port_size,
    );

    // CSG Difference: wall − cutter → wall with rectangular port hole
    println!("  CSG: cutting inlet port through left wall ...");
    let left = csg_boolean_indexed(BooleanOp::Difference, &left_solid, &inlet_cutter)
        .map_err(|e| format!("Inlet CSG error: {e:?}"))?;
    println!(
        "    Left wall: {} faces → {} faces (hole cut)",
        left_solid.face_count(),
        left.face_count()
    );

    println!("  CSG: cutting outlet port through right wall ...");
    let right = csg_boolean_indexed(BooleanOp::Difference, &right_solid, &outlet_cutter)
        .map_err(|e| format!("Outlet CSG error: {e:?}"))?;
    println!(
        "    Right wall: {} faces → {} faces (hole cut)",
        right_solid.face_count(),
        right.face_count()
    );

    // Top plate (lid) — exported separately
    let lid = build_cube(0.0, 0.0, wall + cavity_z, outer_w, outer_h, wall);

    let body_faces: usize = [&bottom, &left, &right, &front, &back]
        .iter()
        .map(|m| m.face_count())
        .sum();
    println!("  Body parts: {body_faces} faces (bottom + 4 walls with port holes)");

    // ── 3. Build TPMS Gyroid lattice inside the cavity ───────────────────────
    //
    // Gyroid level-set:
    //   f(x,y,z) = sin(kx)·cos(ky) + sin(ky)·cos(kz) + sin(kz)·cos(kx)
    //   k = 2π/λ, iso = 0 → minimal surface mid-sheet
    //
    // `build_tpms_box` uses marching cubes with SDF boundary capping to produce
    // a watertight closed surface within the given AABB bounds.
    let lambda = 5.0;
    let resolution = 64;

    println!("Building Gyroid TPMS lattice ...");
    println!("  λ = {lambda:.1} mm, resolution = {resolution}, iso = 0.0");

    let bounds = [
        wall,
        wall,
        wall,
        wall + cavity_w,
        wall + cavity_h,
        wall + cavity_z,
    ];

    let params = TpmsBoxParams {
        bounds,
        period: lambda,
        resolution,
        iso_value: 0.0,
    };

    let tpms_mesh =
        build_tpms_box(&Gyroid, &params).map_err(|e| format!("TPMS build error: {e:?}"))?;

    let tpms_vol = tpms_mesh.signed_volume();
    println!(
        "  TPMS mesh : {:>7} verts, {:>7} faces, vol = {:>10.3} mm³",
        tpms_mesh.vertex_count(),
        tpms_mesh.face_count(),
        tpms_vol,
    );

    // ── 4. Assemble chip: merge body + lattice into one mesh ─────────────────
    println!("Assembling chip ...");
    let mut assembly = IndexedMesh::new();
    merge_into(&mut assembly, &bottom);
    merge_into(&mut assembly, &left);
    merge_into(&mut assembly, &right);
    merge_into(&mut assembly, &front);
    merge_into(&mut assembly, &back);
    merge_into(&mut assembly, &tpms_mesh);

    println!(
        "  Assembly  : {:>7} verts, {:>7} faces",
        assembly.vertex_count(),
        assembly.face_count(),
    );
    println!(
        "  Lid       : {:>7} verts, {:>7} faces",
        lid.vertex_count(),
        lid.face_count(),
    );

    // ── 5. Export STL ────────────────────────────────────────────────────────
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("tpms_blood_chip");
    fs::create_dir_all(&out_dir)?;

    // 5a. Main assembly — open-top chip with TPMS visible inside
    let asm_path = out_dir.join("chip_assembly.stl");
    stl::write_stl_binary(
        &mut BufWriter::new(fs::File::create(&asm_path)?),
        &assembly,
    )?;

    // 5b. Lid (top plate) — load alongside assembly for closed chip view
    let lid_path = out_dir.join("chip_lid.stl");
    stl::write_stl_binary(&mut BufWriter::new(fs::File::create(&lid_path)?), &lid)?;

    // 5c. Standalone TPMS lattice — for detailed lattice inspection
    let lattice_path = out_dir.join("tpms_lattice.stl");
    stl::write_stl_binary(
        &mut BufWriter::new(fs::File::create(&lattice_path)?),
        &tpms_mesh,
    )?;

    // ── 6. Physical summary ──────────────────────────────────────────────────
    let cavity_vol = cavity_w * cavity_h * cavity_z;
    let substrate_vol = outer_w * outer_h * chip_z;
    let porosity = tpms_vol.abs() / cavity_vol * 100.0;
    println!();
    println!("Physical summary:");
    println!("  Substrate vol   : {substrate_vol:.1} mm³");
    println!("  Cavity vol      : {cavity_vol:.1} mm³");
    println!("  TPMS fluid vol  : {:.1} mm³", tpms_vol.abs());
    println!("  Cavity porosity : {porosity:.1}%");
    println!();
    println!("STL files exported to: {}", out_dir.display());
    println!("  chip_assembly.stl — open-top chip body with port holes + Gyroid lattice");
    println!("  chip_lid.stl      — top plate (load with assembly for closed view)");
    println!("  tpms_lattice.stl  — standalone Gyroid lattice (detail inspection)");
    println!();
    println!(
        "View chip_assembly.stl to see the millifluidic chip with TPMS lattice"
    );
    println!("visible inside the cavity, with inlet/outlet port holes through the walls.");

    Ok(())
}
