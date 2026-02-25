//! Export: All 6 millifluidic therapy designs → `outputs/millifluidic/`
//!
//! Runs each of the 6 validated therapy design presets through
//! [`BlueprintMeshPipeline`], asserts watertightness, and writes binary STL
//! files plus schematics produced by `cfd-schematics` for:
//!
//! - **Fluid-domain mesh** (`*_fluid.stl`) — the channel interior used by CFD
//!   solvers.
//! - **Chip body mesh** (`*_chip.stl`) — the PDMS substrate minus the channel
//!   void, suitable for 3-D printing or CAM toolpath generation.
//! - **Schematic JSON** (`*_schematic.json`) — interchange format from
//!   `cfd-schematics` with explicit centrelines, profiles and dimensions.
//! - **Schematic SVG** (`*_schematic.svg`) — 2-D top-down channel layout
//!   rendered by `cfd-schematics`.
//! - **Schematic PNG** (`*_schematic.png`) — raster version of the same layout.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --features scheme-io --example millifluidic_designs
//! ```

use std::fs;
use std::io::BufWriter;
use std::path::Path;

use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_mesh::infrastructure::io::stl;

use cfd_schematics::config::{ChannelTypeConfig, FrustumConfig, GeometryConfig, SerpentineConfig};
use cfd_schematics::geometry::{create_geometry, ChannelSystem, SplitType};
use cfd_schematics::interface::presets::{
    serpentine_chain, serpentine_rect, symmetric_bifurcation, symmetric_trifurcation,
    venturi_chain, venturi_rect,
};
use cfd_schematics::plot_geometry;

// ── Chip dimensions (SBS 96-well plate footprint, mm) ────────────────────────
const CHIP_W_MM: f64 = 127.76;
const CHIP_D_MM: f64 = 85.47;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("millifluidic");
    fs::create_dir_all(&out_dir)?;

    let config = PipelineConfig::default(); // include_chip_body = true

    // Six therapy designs — parameters match the integration tests in
    // tests/blueprint_pipeline.rs so all are known-good.
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
        ("venturi_rect", venturi_rect("d5", 0.004, 0.002, 0.004, 0.005)),
        (
            "serpentine_rect",
            serpentine_rect("d6", 3, 0.010, 0.004, 0.004),
        ),
    ];

    let n = designs.len();
    println!("=================================================================");
    println!("  Millifluidic Therapy Design STL Export");
    println!("  Output : {}", out_dir.display());
    println!("=================================================================");

    let mut stl_count = 0_usize;
    let mut schema_count = 0_usize;

    for (i, (name, bp)) in designs.iter().enumerate() {
        println!();
        let mut out = BlueprintMeshPipeline::run(bp, &config)
            .map_err(|e| format!("{name}: pipeline failed — {e}"))?;

        // ── Fluid mesh ────────────────────────────────────────────────────────
        assert!(
            out.fluid_mesh.is_watertight(),
            "{name}: fluid mesh is not watertight"
        );
        assert!(
            out.fluid_mesh.signed_volume() > 0.0,
            "{name}: fluid mesh has non-positive volume"
        );

        let topo_str = format!("{:?}", out.topology_class);

        println!(
            "  [{}/{}] {}  ({}, {} segments)",
            i + 1,
            n,
            name,
            topo_str,
            out.segment_count
        );
        println!(
            "    Fluid  : {:>6} faces,  vol = {:>10.3} mm³  [watertight]",
            out.fluid_mesh.face_count(),
            out.fluid_mesh.signed_volume()
        );

        let fluid_path = out_dir.join(format!("{name}_fluid.stl"));
        let file = fs::File::create(&fluid_path)?;
        stl::write_binary_stl(
            &mut BufWriter::new(file),
            &out.fluid_mesh.vertices,
            &out.fluid_mesh.faces,
        )?;
        println!("    -> {}", fluid_path.file_name().unwrap().to_string_lossy());
        stl_count += 1;

        // ── Chip body ─────────────────────────────────────────────────────────
        if let Some(chip) = out.chip_mesh.as_mut() {
            assert!(
                chip.is_watertight(),
                "{name}: chip mesh is not watertight"
            );
            assert!(
                chip.signed_volume() > 0.0,
                "{name}: chip mesh has non-positive volume"
            );
            println!(
                "    Chip   : {:>6} faces,  vol = {:>10.3} mm³  [watertight]",
                chip.face_count(),
                chip.signed_volume()
            );

            let chip_path = out_dir.join(format!("{name}_chip.stl"));
            let file = fs::File::create(&chip_path)?;
            stl::write_binary_stl(
                &mut BufWriter::new(file),
                &chip.vertices,
                &chip.faces,
            )?;
            println!("    -> {}", chip_path.file_name().unwrap().to_string_lossy());
            stl_count += 1;
        }

        // ── 2-D schematics via cfd-schematics ─────────────────────────────────
        // JSON interchange format (centrelines + profiles), SVG vector layout,
        // and PNG raster layout — all produced by cfd-schematics.
        let system = channel_system_for(name);

        let json_path = out_dir.join(format!("{name}_schematic.json"));
        fs::write(&json_path, system.to_interchange_json()?)?;
        println!("    -> {}", json_path.file_name().unwrap().to_string_lossy());
        schema_count += 1;

        let svg_path = out_dir.join(format!("{name}_schematic.svg"));
        plot_geometry(&system, svg_path.to_str().unwrap())?;
        println!("    -> {}", svg_path.file_name().unwrap().to_string_lossy());
        schema_count += 1;

        let png_path = out_dir.join(format!("{name}_schematic.png"));
        plot_geometry(&system, png_path.to_str().unwrap())?;
        println!("    -> {}", png_path.file_name().unwrap().to_string_lossy());
        schema_count += 1;
    }

    println!();
    println!("=================================================================");
    println!("  {stl_count} STL files written to {}", out_dir.display());
    println!("  {schema_count} schematic files (JSON + SVG + PNG) written");
    println!("=================================================================");
    Ok(())
}

// ── 2-D schematic geometry ────────────────────────────────────────────────────

/// Build a [`ChannelSystem`] for the named therapy design.
///
/// All dimensions use millimetres on the SBS 96-well plate footprint
/// (127.76 × 85.47 mm).  The geometry is fed to `cfd-schematics` to produce
/// the interchange JSON, SVG, and PNG schematic files.
fn channel_system_for(name: &str) -> ChannelSystem {
    let box_dims = (CHIP_W_MM, CHIP_D_MM);

    // 4 mm circular cross-section for all designs (main channel diameter).
    let geom = GeometryConfig {
        channel_width: 4.0,
        channel_height: 4.0,
        ..GeometryConfig::default()
    };

    match name {
        // ── Venturi designs: tapered (frustum) channel, no branching ──────────
        "venturi_chain" | "venturi_rect" => {
            let frustum = FrustumConfig {
                inlet_width: 4.0,   // mm — matches inlet_diameter_m = 0.004 m
                throat_width: 2.0,  // mm — matches throat_diameter_m = 0.002 m
                outlet_width: 4.0,
                ..FrustumConfig::default()
            };
            create_geometry(box_dims, &[], &geom, &ChannelTypeConfig::AllFrustum(frustum))
        }

        // ── Bifurcation: one symmetric Y-split ────────────────────────────────
        "symmetric_bifurcation" => create_geometry(
            box_dims,
            &[SplitType::Bifurcation],
            &geom,
            &ChannelTypeConfig::AllStraight,
        ),

        // ── Trifurcation: one symmetric T-split ───────────────────────────────
        "symmetric_trifurcation" => create_geometry(
            box_dims,
            &[SplitType::Trifurcation],
            &geom,
            &ChannelTypeConfig::AllStraight,
        ),

        // ── Serpentine designs: winding single channel, no branching ──────────
        _ => create_geometry(
            box_dims,
            &[],
            &geom,
            &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
        ),
    }
}
