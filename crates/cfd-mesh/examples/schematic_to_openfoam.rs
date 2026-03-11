//! Canonical JSON schematic → STL + OpenFOAM pipeline.
//!
//! Demonstrates the complete path from `cfd-schematics` designs to watertight
//! 3-D surface meshes ready for:
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
//! NetworkBlueprint / ChannelSystem (cfd-schematics)
//!   └─▶ BlueprintMeshPipeline::run()      (default for non-serpentine presets)
//!      or mesh_output_from_blueprint() (serpentine presets; centerline-driven)
//!         ├─ fluid_mesh  — channel interior (IndexedMesh, boundary-labelled)
//!         └─ chip_mesh   — PDMS substrate minus channel voids
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

use cfd_mesh::application::channel::path::ChannelPath;
use cfd_mesh::application::channel::substrate::SubstrateBuilder;
use cfd_mesh::application::channel::sweep::SweepMesher;
use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::application::pipeline::PipelineOutput;
use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig, PipelineVolumeTrace};
use cfd_mesh::domain::core::index::RegionId;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::mesh::IndexedMesh;
use cfd_mesh::domain::topology::halfedge::PatchType;
use cfd_mesh::infrastructure::io::openfoam::write_openfoam_polymesh;
use cfd_mesh::infrastructure::io::scheme;
use cfd_mesh::infrastructure::io::stl::write_stl_binary;

use cfd_schematics::config::{ChannelTypeConfig, FrustumConfig, GeometryConfig, SerpentineConfig};
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;
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
        if *name == "serpentine_chain" || *name == "serpentine_rect" {
            let blueprint = blueprint_for(name);
            let force_circular = *name == "serpentine_chain";
            out = mesh_output_from_blueprint(&blueprint, &config, force_circular)
                .map_err(|e| format!("{name}: schematic mesh failed — {e}"))?;
        }

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
        let fluid_vol = out.fluid_mesh.signed_volume();
        let fluid_path = design_dir.join(format!("{name}_fluid.stl"));
        let mut f = BufWriter::new(fs::File::create(&fluid_path)?);
        write_stl_binary(&mut f, &out.fluid_mesh)?;
        println!(
            "         fluid STL  → {:>7} faces, vol = {:>10.3} mm³",
            out.fluid_mesh.face_count(),
            fluid_vol
        );
        stl_count += 1;

        // ── STL: chip body ────────────────────────────────────────────────────
        if let Some(chip) = out.chip_mesh.as_mut() {
            assert!(chip.is_watertight(), "{name}: chip mesh must be watertight");
            let chip_vol = chip.signed_volume();
            assert!(chip_vol > 0.0, "{name}: chip mesh volume must be positive");
            let chip_path = design_dir.join(format!("{name}_chip.stl"));
            let mut f = BufWriter::new(fs::File::create(&chip_path)?);
            write_stl_binary(&mut f, chip)?;
            println!(
                "         chip  STL  → {:>7} faces, vol = {:>10.3} mm³",
                chip.face_count(),
                chip_vol
            );
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

fn blueprint_for(name: &str) -> cfd_schematics::NetworkBlueprint {
    const CHIP_W_MM: f64 = 127.76;
    const CHIP_D_MM: f64 = 85.47;
    let box_dims = (CHIP_W_MM, CHIP_D_MM);

    let geom = GeometryConfig {
        channel_width: 4.0,
        channel_height: 4.0,
        ..GeometryConfig::default()
    };

    match name {
        "venturi_chain" | "venturi_rect" => {
            let frustum = FrustumConfig {
                inlet_width: 4.0,
                throat_width: 2.0,
                outlet_width: 4.0,
                ..FrustumConfig::default()
            };
            create_geometry(
                box_dims,
                &[],
                &geom,
                &ChannelTypeConfig::AllFrustum(frustum),
            )
        }
        "symmetric_bifurcation" => create_geometry(
            box_dims,
            &[SplitType::Bifurcation],
            &geom,
            &ChannelTypeConfig::AllStraight,
        ),
        "symmetric_trifurcation" => create_geometry(
            box_dims,
            &[SplitType::Trifurcation],
            &geom,
            &ChannelTypeConfig::AllStraight,
        ),
        _ => create_geometry(
            box_dims,
            &[],
            &geom,
            &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
        ),
    }
}

fn mesh_output_from_blueprint(
    system: &cfd_schematics::NetworkBlueprint,
    config: &PipelineConfig,
    force_circular: bool,
) -> Result<PipelineOutput, Box<dyn std::error::Error>> {
    let schematic3d = scheme::from_blueprint(
        system,
        config.chip_height_mm as Real,
        config.circular_segments,
    )?;

    let mesher = SweepMesher::new();
    let mut all_channels: Option<IndexedMesh> = None;
    let mut all_void_channels: Option<IndexedMesh> = None;

    for channel_def in &schematic3d.channels {
        let profile = if force_circular {
            circularized_profile(&channel_def.profile, config.circular_segments)
        } else {
            channel_def.profile.clone()
        };

        let mut current = sweep_channel(
            &mesher,
            &profile,
            &channel_def.path,
            channel_def.width_scales.as_deref(),
        );
        current.rebuild_edges();

        all_channels = Some(if let Some(existing) = all_channels.take() {
            csg_boolean(BooleanOp::Union, &existing, &current)?
        } else {
            current
        });

        if config.include_chip_body {
            let extension_mm = profile_radius_mm(&profile).max(0.25);
            let extended_path = extend_path_ends(&channel_def.path, extension_mm);
            let mut void_current = sweep_channel(
                &mesher,
                &profile,
                &extended_path,
                channel_def.width_scales.as_deref(),
            );
            void_current.rebuild_edges();

            all_void_channels = Some(if let Some(existing) = all_void_channels.take() {
                csg_boolean(BooleanOp::Union, &existing, &void_current)?
            } else {
                void_current
            });
        }
    }

    let mut fluid_mesh = all_channels.ok_or("no channels produced from schematic")?;
    fluid_mesh.orient_outward();
    fluid_mesh.retain_largest_component();
    fluid_mesh.rebuild_edges();

    label_inlet_outlet_from_system(&mut fluid_mesh, &schematic3d)?;

    let chip_mesh = if config.include_chip_body {
        let substrate = SubstrateBuilder::well_plate_96(config.chip_height_mm).build_indexed()?;
        let void_mesh = all_void_channels.as_ref().unwrap_or(&fluid_mesh);
        let mut chip = csg_boolean(BooleanOp::Difference, &substrate, void_mesh)?;
        chip.orient_outward();
        chip.retain_largest_component();
        chip.rebuild_edges();

        if !chip.is_watertight() {
            chip = csg_boolean(BooleanOp::Difference, &substrate, void_mesh)?;
            chip.orient_outward();
            chip.retain_largest_component();
            chip.rebuild_edges();
        }
        if !chip.is_watertight() {
            return Err("schematic-driven chip body is not watertight".into());
        }
        Some(chip)
    } else {
        None
    };
    let schematic_summary = system.fluid_volume_summary();
    let fluid_mesh_volume_mm3 = fluid_mesh.signed_volume().abs();
    let chip_mesh_volume_mm3 = chip_mesh.as_ref().map(|mesh| mesh.signed_volume().abs());
    let fluid_mesh_volume_error_mm3 =
        fluid_mesh_volume_mm3 - schematic_summary.total_fluid_volume_mm3;
    let fluid_mesh_volume_error_pct = if schematic_summary.total_fluid_volume_mm3.abs() <= 1e-18 {
        0.0
    } else {
        fluid_mesh_volume_error_mm3.abs() / schematic_summary.total_fluid_volume_mm3.abs() * 100.0
    };
    let volume_trace = PipelineVolumeTrace {
        schematic_summary,
        channel_traces: Vec::new(),
        pre_csg_channel_volume_mm3: fluid_mesh_volume_mm3,
        synthetic_connector_volume_mm3: 0.0,
        fluid_mesh_volume_mm3,
        chip_mesh_volume_mm3,
        fluid_mesh_volume_error_mm3,
        fluid_mesh_volume_error_pct,
        csg_overlap_delta_mm3: 0.0,
    };

    Ok(PipelineOutput {
        fluid_mesh,
        chip_mesh,
        topology_class: cfd_mesh::application::pipeline::TopologyClass::LinearChain {
            n_segments: system.channels.len(),
        },
        segment_count: system.channels.len(),
        layout_segments: Vec::new(),
        volume_trace,
    })
}

fn label_inlet_outlet_from_system(
    mesh: &mut IndexedMesh,
    schematic3d: &scheme::Schematic,
) -> Result<(), Box<dyn std::error::Error>> {
    let ch = schematic3d
        .channels
        .first()
        .ok_or("cannot label boundaries: no channels in schematic")?;
    let pts = ch.path.points();
    let inlet = *pts.first().ok_or("channel path missing inlet point")?;
    let outlet = *pts.last().ok_or("channel path missing outlet point")?;
    let eps = profile_radius_mm(&ch.profile) * 2.0;

    let face_ids: Vec<_> = mesh.faces.iter_enumerated().map(|(id, _)| id).collect();
    let mut labels: Vec<_> = Vec::new();
    for fid in face_ids {
        let face = mesh.faces.get(fid);
        let p0 = mesh.vertices.position(face.vertices[0]);
        let p1 = mesh.vertices.position(face.vertices[1]);
        let p2 = mesh.vertices.position(face.vertices[2]);
        let c = Point3r::new(
            (p0.x + p1.x + p2.x) / 3.0,
            (p0.y + p1.y + p2.y) / 3.0,
            (p0.z + p1.z + p2.z) / 3.0,
        );

        if (c - inlet).norm() <= eps {
            labels.push((fid, "inlet"));
        } else if (c - outlet).norm() <= eps {
            labels.push((fid, "outlet"));
        }
    }

    mesh.boundary_labels.clear();
    for (fid, label) in labels {
        mesh.mark_boundary(fid, label);
    }
    Ok(())
}

fn sweep_channel(
    mesher: &SweepMesher,
    profile: &cfd_mesh::application::channel::profile::ChannelProfile,
    path: &ChannelPath,
    scales: Option<&[Real]>,
) -> IndexedMesh {
    let mut mesh = IndexedMesh::new();
    let faces = if let Some(scales) = scales {
        mesher.sweep_variable(profile, path, scales, &mut mesh.vertices, RegionId::new(0))
    } else {
        mesher.sweep(profile, path, &mut mesh.vertices, RegionId::new(0))
    };
    for face in faces {
        mesh.faces.push(face);
    }
    mesh
}

fn profile_radius_mm(profile: &cfd_mesh::application::channel::profile::ChannelProfile) -> Real {
    use cfd_mesh::application::channel::profile::ChannelProfile;
    match profile {
        ChannelProfile::Circular { radius, .. } => *radius,
        ChannelProfile::Rectangular { width, height } => 0.5 * width.min(*height),
        ChannelProfile::RoundedRectangular { width, height, .. } => 0.5 * width.min(*height),
    }
}

fn extend_path_ends(path: &ChannelPath, extension_mm: Real) -> ChannelPath {
    let pts = path.points();
    if pts.len() < 2 || extension_mm <= 0.0 {
        return path.clone();
    }

    let first_dir = (pts[1] - pts[0]).normalize();
    let last_dir = (pts[pts.len() - 1] - pts[pts.len() - 2]).normalize();

    let mut out = pts.to_vec();
    out[0] = pts[0] - first_dir * extension_mm;
    let n = out.len() - 1;
    out[n] = pts[n] + last_dir * extension_mm;
    ChannelPath::new(out)
}

fn circularized_profile(
    profile: &cfd_mesh::application::channel::profile::ChannelProfile,
    segments: usize,
) -> cfd_mesh::application::channel::profile::ChannelProfile {
    use cfd_mesh::application::channel::profile::ChannelProfile;
    use std::f64::consts::PI;

    match profile {
        ChannelProfile::Circular { radius, .. } => ChannelProfile::Circular {
            radius: *radius,
            segments,
        },
        ChannelProfile::Rectangular { width, height } => {
            let d_h = if *width > 0.0 && *height > 0.0 {
                2.0 * *width * *height / (*width + *height)
            } else {
                width.min(*height).max(0.0)
            };
            ChannelProfile::Circular {
                radius: 0.5 * d_h,
                segments,
            }
        }
        ChannelProfile::RoundedRectangular {
            width,
            height,
            corner_radius,
            ..
        } => {
            let area = (*width * *height - (4.0 - PI) * *corner_radius * *corner_radius).max(0.0);
            let d_eq = if area > 0.0 {
                2.0 * (area / PI).sqrt()
            } else {
                width.min(*height).max(0.0)
            };
            ChannelProfile::Circular {
                radius: 0.5 * d_eq,
                segments,
            }
        }
    }
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
