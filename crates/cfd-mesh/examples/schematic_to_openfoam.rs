//! Canonical JSON schematic → STL + OpenFOAM pipeline.
//!
//! Demonstrates the complete path from `cfd-schematics`-produced `NetworkBlueprint`
//! presets to watertight 3-D surface meshes ready for:
//!
//! - **Manufacturing**: binary STL files (`*_fluid.stl`, `*_chip.stl`) for 3-D
//!   printing or CNC machining.
//!
//! - **3-D CFD simulation**: OpenFOAM `constant/polyMesh/` directories with
//!   named boundary patches (`inlet`, `outlet`, `walls`) accepted directly by
//!   `simpleFoam`, `icoFoam`, and `snappyHexMesh`.
//!
//! ## Pipeline
//!
//! ```text
//! NetworkBlueprint (cfd-schematics)
//!   └─▶ BlueprintMeshPipeline::run()
//!         ├─ fluid_mesh  —  channel interior (IndexedMesh, boundary-labelled)
//!         └─ chip_mesh   —  PDMS substrate minus channel voids
//!               │
//!               ├─▶ write_stl_binary()         → *_fluid.stl / *_chip.stl
//!               └─▶ reassign_regions_from_labels()
//!                     └─▶ write_openfoam_polymesh() → constant/polyMesh/
//! ```
//!
//! ## Boundary label → OpenFOAM patch mapping
//!
//! `BlueprintMeshPipeline` marks faces via `IndexedMesh::mark_boundary()`:
//!
//! | Label string | RegionId | PatchType  | OpenFOAM `type` |
//! |---|---|---|---|
//! | `"inlet"`  | 0 | `Inlet`  | `patch` (physicalType inlet) |
//! | `"outlet"` | 1 | `Outlet` | `patch` (physicalType outlet) |
//! | `"wall"`   | 2 | `Wall`   | `wall` |
//! | unlabelled | 2 | `Wall`   | `wall` (defaultFaces) |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example schematic_to_openfoam
//! ```

use std::fs;
use std::io::BufWriter;
use std::path::Path;

use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_mesh::domain::core::index::RegionId;
use cfd_mesh::domain::mesh::IndexedMesh;
use cfd_mesh::domain::topology::halfedge::PatchType;
use cfd_mesh::infrastructure::io::openfoam::write_openfoam_polymesh;
use cfd_mesh::infrastructure::io::stl::write_stl_binary;

use cfd_schematics::interface::presets::{
    serpentine_chain, serpentine_rect, symmetric_bifurcation, symmetric_trifurcation,
    venturi_chain, venturi_rect,
};

// ── Region ID constants ───────────────────────────────────────────────────────

const REGION_INLET: u32 = 0;
const REGION_OUTLET: u32 = 1;
const REGION_WALL: u32 = 2;

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════╗");
    println!("║  Schematic → 3D Mesh → STL + OpenFOAM Pipeline  ║");
    println!("╚══════════════════════════════════════════════════╝");

    let out_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("schematic_to_openfoam");
    fs::create_dir_all(&out_root)?;

    // Six validated millifluidic therapy designs
    let designs: Vec<(&str, cfd_schematics::NetworkBlueprint)> = vec![
        ("venturi_chain", venturi_chain("d1", 0.030, 0.004, 0.002)),
        (
            "symmetric_bifurcation",
            symmetric_bifurcation("d2", 0.010, 0.010, 0.004, 0.003),
        ),
        (
            "symmetric_trifurcation",
            symmetric_trifurcation("d3", 0.010, 0.008, 0.004, 0.004),
        ),
        ("serpentine_chain", serpentine_chain("d4", 3, 0.010, 0.004)),
        (
            "venturi_rect",
            venturi_rect("d5", 0.004, 0.002, 0.004, 0.005),
        ),
        (
            "serpentine_rect",
            serpentine_rect("d6", 3, 0.010, 0.004, 0.004),
        ),
    ];

    let config = PipelineConfig::default(); // chip_height=10mm, circular_segments=16
    let n = designs.len();
    let mut stl_count = 0usize;
    let mut of_count = 0usize;

    println!();
    println!(
        "  Running {} designs with chip_height = {} mm",
        n, config.chip_height_mm
    );
    println!("  Output root: {}", out_root.display());
    println!("{}", "─".repeat(60));

    for (i, (name, bp)) in designs.iter().enumerate() {
        print!("  [{}/{}] {:.<40}", i + 1, n, format!("{name} "));

        let mut out = BlueprintMeshPipeline::run(bp, &config)
            .map_err(|e| format!("{name}: pipeline failed — {e}"))?;

        // Invariant checks
        assert!(
            out.fluid_mesh.is_watertight(),
            "{name}: fluid mesh must be watertight"
        );
        assert!(
            out.fluid_mesh.signed_volume() > 0.0,
            "{name}: fluid mesh must have positive volume"
        );

        let topo = format!("{:?}", out.topology_class);
        println!(" {topo:.<22} {} faces", out.fluid_mesh.face_count());

        let design_dir = out_root.join(name);
        fs::create_dir_all(&design_dir)?;

        // ── STL: fluid mesh ───────────────────────────────────────────────────
        let fluid_path = design_dir.join(format!("{name}_fluid.stl"));
        let mut f = BufWriter::new(fs::File::create(&fluid_path)?);
        write_stl_binary(&mut f, &out.fluid_mesh)?;
        println!(
            "         fluid STL  → {} faces",
            out.fluid_mesh.face_count()
        );
        stl_count += 1;

        // ── STL: chip body ────────────────────────────────────────────────────
        if let Some(chip) = out.chip_mesh.as_mut() {
            assert!(chip.is_watertight(), "{name}: chip mesh must be watertight");
            assert!(
                chip.signed_volume() > 0.0,
                "{name}: chip mesh volume must be positive"
            );
            let chip_path = design_dir.join(format!("{name}_chip.stl"));
            let mut f = BufWriter::new(fs::File::create(&chip_path)?);
            write_stl_binary(&mut f, chip)?;
            println!("         chip  STL  → {} faces", chip.face_count());
            stl_count += 1;
        }

        // ── OpenFOAM: fluid mesh with boundary patches ────────────────────────
        //
        // BlueprintMeshPipeline labels faces via `mesh.mark_boundary(fid, label)`:
        //   "inlet"  / "outlet" / "wall"
        //
        // `write_openfoam_polymesh` partitions faces by `face.region: RegionId`.
        // We must re-assign these `RegionId`s from the string labels because CSG
        // operations reset region IDs to their sweep-time defaults.
        let mut fluid_for_of = out.fluid_mesh.clone();
        reassign_regions_from_labels(&mut fluid_for_of);

        let of_dir = design_dir.join("constant/polyMesh");
        write_openfoam_polymesh(
            &fluid_for_of,
            &of_dir,
            &[
                (RegionId::new(REGION_INLET), "inlet", PatchType::Inlet),
                (RegionId::new(REGION_OUTLET), "outlet", PatchType::Outlet),
                (RegionId::new(REGION_WALL), "walls", PatchType::Wall),
            ],
        )?;
        println!("         OpenFOAM   → {}/", of_dir.display());
        of_count += 1;

        // ── OpenFOAM: chip body (wall-only surface for snappyHexMesh) ─────────
        if let Some(chip) = out.chip_mesh.as_ref() {
            let chip_of_dir = design_dir.join("chip_polyMesh");
            write_openfoam_polymesh(
                chip,
                &chip_of_dir,
                &[(RegionId::new(0), "walls", PatchType::Wall)],
            )?;
            println!("         chip OF    → {}/", chip_of_dir.display());
            of_count += 1;
        }
    }

    println!("{}", "─".repeat(60));
    println!("  {stl_count} STL files + {of_count} OpenFOAM polyMesh dirs written");
    println!("  Output: {}", out_root.display());
    println!();

    Ok(())
}

// ── Boundary label → RegionId remapping ──────────────────────────────────────

/// Re-assign every face's `RegionId` in `mesh` from its `boundary_labels` entry.
///
/// This bridges [`BlueprintMeshPipeline`]'s string boundary labels
/// (`"inlet"` / `"outlet"` / `"wall"`) into the `RegionId` integers that
/// `write_openfoam_polymesh` uses to partition faces into named boundary patches.
///
/// CSG union/difference operations overwrite the face `region` field with sweep-time
/// bookkeeping IDs; this function restores the semantically correct IDs from the
/// separately stored `boundary_labels` map.
///
/// # Region assignment
///
/// | Label string | `RegionId` constant |
/// |---|---|
/// | `"inlet"`  | `REGION_INLET  (0)` |
/// | `"outlet"` | `REGION_OUTLET (1)` |
/// | all others | `REGION_WALL   (2)` |
fn reassign_regions_from_labels(mesh: &mut IndexedMesh) {
    use cfd_mesh::domain::core::index::FaceId;
    use std::collections::HashMap;

    // Build FaceId → RegionId lookup from the immutable label map first,
    // so the later mutable face iteration has no aliased borrows.
    let region_for: HashMap<FaceId, RegionId> = mesh
        .boundary_labels
        .iter()
        .map(|(&fid, label)| {
            let region = match label.as_str() {
                "inlet" => RegionId::new(REGION_INLET),
                "outlet" => RegionId::new(REGION_OUTLET),
                _ => RegionId::new(REGION_WALL),
            };
            (fid, region)
        })
        .collect();

    // Single mutable pass: default all faces to wall, then apply label-derived regions.
    for (fid, face) in mesh.faces.iter_mut_enumerated() {
        face.region = *region_for.get(&fid).unwrap_or(&RegionId::new(REGION_WALL));
    }
}
