//! Millifluidic chip STL export — schematic JSON → binary STL.
//!
//! For each of the canonical millifluidic therapy topologies this example:
//!
//! 1. Builds a [`NetworkBlueprint`] via the `cfd-schematics` preset library.
//! 2. Runs [`BlueprintMeshPipeline`] to produce:
//!    - **fluid mesh** — the interior channel void (inlet + outlet labelled).
//!    - **chip body** — the SBS-96 cuboid substrate with channel voids subtracted.
//! 3. Serializes the [`ChannelSystem`] 2-D schematic to JSON, SVG, and PNG.
//! 4. Writes all meshes to binary STL under `outputs/millifluidic_chip_stl/{design}/`.
//!
//! ## Topologies
//!
//! | Name                     | Description                                  |
//! |--------------------------|----------------------------------------------|
//! | `venturi_chain`          | Inlet → taper → throat → expand → outlet     |
//! | `symmetric_bifurcation`  | Inlet → Y-split (mirrored arms) → outlet     |
//! | `symmetric_trifurcation` | Inlet → T-split (mirrored arms) → outlet     |
//! | `serpentine_chain`       | Inlet → winding 3-leg serpentine → outlet    |
//!
//! Each topology has a **mirrored** channel network: the upper and lower fork
//! arms (bifurcation / trifurcation) are symmetric about the chip centre-line.
//! Inlet and outlet ports exit the chip body face and are labelled in the fluid
//! mesh boundary.
//!
//! ## Usage
//!
//! ```sh
//! cargo run -p cfd-mesh --features scheme-io --example millifluidic_chip_stl
//! ```

use std::fs;
use std::io::BufWriter;
use std::path::Path;

use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_mesh::infrastructure::io::stl;

use cfd_schematics::config::{ChannelTypeConfig, FrustumConfig, GeometryConfig, SerpentineConfig};
use cfd_schematics::geometry::{create_geometry, ChannelSystem, SplitType};
use cfd_schematics::interface::presets::{
    serpentine_chain, symmetric_bifurcation, symmetric_trifurcation, venturi_chain,
};
use cfd_schematics::plot_geometry;

// ── SBS 96-well plate footprint (mm) ─────────────────────────────────────────
/// Standard SBS micro-plate footprint.  All inlet/outlet ports are drilled at
/// the left face (x = 0) and right face (x = CHIP_W_MM) of the substrate.
const CHIP_W_MM: f64 = 127.76;
const CHIP_D_MM: f64 = 85.47;

// ── Design catalogue ──────────────────────────────────────────────────────────

/// One STL export target.
struct Design {
    /// Filesystem-safe short identifier used for filenames and directories.
    name:        &'static str,
    /// Human-readable description printed to stdout.
    description: &'static str,
    /// The [`NetworkBlueprint`] driving the mesh pipeline.
    blueprint:   cfd_schematics::NetworkBlueprint,
    /// The 2-D [`ChannelSystem`] used for schematic JSON/SVG/PNG export.
    system:      ChannelSystem,
    /// Per-design pipeline config override.
    ///
    /// Most designs use the shared `config` from `main`.  Designs with known
    /// CSG limitations (e.g. trifurcation chip-body 3-way T-junction) can
    /// override `include_chip_body` to `false` here.
    config_override: Option<PipelineConfig>,
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════╗");
    println!("║  Millifluidic Chip STL Export (schematic → binary STL) ║");
    println!("╚══════════════════════════════════════════════════════╝");
    println!();

    let designs = build_designs();
    let n = designs.len();

    let out_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("millifluidic_chip_stl");

    // Pipeline config: include substrate chip body (channels subtracted from
    // the SBS-96 cuboid).  16-segment circular cross-sections give a clean STL
    // at manageable file size for millifluidic scale.
    let config = PipelineConfig {
        include_chip_body: true,
        circular_segments: 16,
        ..PipelineConfig::default()
    };

    let mut stl_total = 0_usize;
    let mut schema_total = 0_usize;

    for (i, design) in designs.iter().enumerate() {
        let out_dir = out_root.join(design.name);
        fs::create_dir_all(&out_dir)?;

        println!("[{}/{}] {}", i + 1, n, design.name);
        println!("        {}", design.description);

        // ── Run mesh pipeline ─────────────────────────────────────────────────
        let effective_config = design.config_override.as_ref().unwrap_or(&config);
        let mut output = BlueprintMeshPipeline::run(&design.blueprint, effective_config)
            .map_err(|e| format!("{}: pipeline failed — {e}", design.name))?;

        // Verify watertightness (fast sanity check before writing STL).
        assert!(
            output.fluid_mesh.is_watertight(),
            "{}: fluid mesh not watertight",
            design.name
        );
        if let Some(chip) = output.chip_mesh.as_mut() {
            assert!(
                chip.is_watertight(),
                "{}: chip body not watertight",
                design.name
            );
        }

        // ── Fluid mesh STL ────────────────────────────────────────────────────
        // The fluid mesh represents the interior void of the microfluidic
        // network.  Boundary faces are labelled "inlet" / "outlet" / "wall"
        // by the pipeline for CFD use.
        {
            let v = output.fluid_mesh.signed_volume();
            println!(
                "        fluid : {:>7} faces,  vol = {:>10.3} mm³  [watertight ✓]",
                output.fluid_mesh.face_count(),
                v,
            );

            let path = out_dir.join(format!("{}_fluid.stl", design.name));
            let f = fs::File::create(&path)?;
            stl::write_binary_stl(
                &mut BufWriter::new(f),
                &output.fluid_mesh.vertices,
                &output.fluid_mesh.faces,
            )?;
            println!("        → {}", path.file_name().unwrap().to_string_lossy());
            stl_total += 1;
        }

        // ── Chip body STL ─────────────────────────────────────────────────────
        // The chip body is the PDMS/PMMA substrate (SBS-96 cuboid) with the
        // channel voids subtracted.  Inlet and outlet ports appear as circular
        // openings on the inlet (x = 0) and outlet (x = chip_width) faces.
        if let Some(chip) = output.chip_mesh.as_mut() {
            let v = chip.signed_volume();
            println!(
                "        chip  : {:>7} faces,  vol = {:>10.3} mm³  [watertight ✓]",
                chip.face_count(),
                v,
            );

            let path = out_dir.join(format!("{}_chip.stl", design.name));
            let f = fs::File::create(&path)?;
            stl::write_binary_stl(
                &mut BufWriter::new(f),
                &chip.vertices,
                &chip.faces,
            )?;
            println!("        → {}", path.file_name().unwrap().to_string_lossy());
            stl_total += 1;
        }

        // ── Schematic JSON ────────────────────────────────────────────────────
        // The 2-D interchange JSON encodes channel centrelines, cross-section
        // profiles, and bounding-box dimensions.  Downstream tools (cfd-mesh,
        // cfd-1d) can reload this JSON to reconstruct the 3-D mesh without
        // re-running the blueprint pipeline.
        {
            let json = design.system.to_interchange_json()?;
            let path = out_dir.join(format!("{}_schematic.json", design.name));
            fs::write(&path, &json)?;
            println!("        → {}", path.file_name().unwrap().to_string_lossy());
            schema_total += 1;
        }

        // ── Schematic SVG ─────────────────────────────────────────────────────
        {
            let path = out_dir.join(format!("{}_schematic.svg", design.name));
            plot_geometry(&design.system, path.to_str().unwrap())?;
            println!("        → {}", path.file_name().unwrap().to_string_lossy());
            schema_total += 1;
        }

        // ── Schematic PNG ─────────────────────────────────────────────────────
        {
            let path = out_dir.join(format!("{}_schematic.png", design.name));
            plot_geometry(&design.system, path.to_str().unwrap())?;
            println!("        → {}", path.file_name().unwrap().to_string_lossy());
            schema_total += 1;
        }

        println!();
    }

    println!("═══════════════════════════════════════════════════════");
    println!("  {stl_total} STL files written");
    println!("  {schema_total} schematic files (JSON + SVG + PNG) written");
    println!("  Output root: {}", out_root.display());
    println!("═══════════════════════════════════════════════════════");
    Ok(())
}

// ── Design catalogue builder ──────────────────────────────────────────────────

/// Assemble all four canonical therapy topologies.
///
/// Each design provides both:
/// - A [`NetworkBlueprint`] (drives [`BlueprintMeshPipeline`] for 3-D STL)
/// - A [`ChannelSystem`] (drives `cfd-schematics` for 2-D schematic export)
///
/// Dimensions follow the SBS 96-well plate footprint.  Channel diameters are
/// chosen to produce clearly visible port openings at the chip face scale:
///
/// | Topology       | d_inlet (mm) | d_throat / d_daughter (mm) |
/// |----------------|:------------:|:--------------------------:|
/// | Venturi chain  |     4.0      |            2.0             |
/// | Bifurcation    |     4.0      |            3.0             |
/// | Trifurcation   |     4.0      |            3.0             |
/// | Serpentine     |     4.0      |            —               |
fn build_designs() -> Vec<Design> {
    let box_dims = (CHIP_W_MM, CHIP_D_MM);

    // Shared geometry config for 2-D schematic rendering.
    let geom_4mm = GeometryConfig {
        channel_width:  4.0,
        channel_height: 4.0,
        ..GeometryConfig::default()
    };

    vec![
        // ── Venturi chain: inlet → contraction → throat → expansion → outlet ──
        //
        // Mirrored about the chip centreline by the preset (the tube is straight
        // along the centreline; "mirroring" is implicit — it's a single central
        // tube with a tapered throat.  The inlet and outlet ports are at x = 0
        // and x = chip_width respectively at y_centre = chip_depth/2.
        Design {
            name: "venturi_chain",
            description:
                "Single central channel: inlet Ø 4 mm → Ø 2 mm throat → Ø 4 mm → outlet  \
                 (127.76 mm long, one circular port on each end face)",
            blueprint: venturi_chain("vc", 0.12776, 0.004, 0.002),
            system:    create_geometry(
                box_dims,
                &[],
                &geom_4mm,
                &ChannelTypeConfig::AllFrustum(FrustumConfig {
                    inlet_width:  4.0,
                    throat_width: 2.0,
                    outlet_width: 4.0,
                    ..FrustumConfig::default()
                }),
            ),
            config_override: None,
        },

        // ── Symmetric bifurcation: mirrored Y-split ───────────────────────────
        //
        // Inlet at x = 0, centre-line.  At 25 % chip length the channel
        // diverges into two symmetric arms (upper + lower, mirrored about the
        // y-centre).  They re-converge at 75 % chip length and exit as a single
        // outlet.  Produces exactly 2 ports: one inlet, one outlet.
        Design {
            name: "symmetric_bifurcation",
            description:
                "Mirrored Y-topology: 1 inlet → upper+lower fork arms → 1 outlet  \
                 (parent Ø 4 mm, daughter Ø 4 mm — full-diameter equal fork)",
            blueprint: symmetric_bifurcation("bf", 0.010, 0.010, 0.004, 0.004),
            system:    create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &geom_4mm,
                &ChannelTypeConfig::AllStraight,
            ),
            config_override: None,
        },

        // ── Symmetric trifurcation: 1 inlet → 3 outlets T-split ──────────────
        //
        // Inlet at x = 0, centre-line.  At 25% chip length the channel diverges
        // into three symmetric arms (upper, centre, lower).  Each arm exits on
        // the right face.  Produces exactly 4 ports: 1 inlet, 3 outlets.
        Design {
            name: "symmetric_trifurcation",
            description:
                "Diamond T-topology: 1 inlet → 3 parallel arms → 1 outlet  \
                 (parent Ø 4 mm, daughter Ø 4 mm)",
            blueprint: symmetric_trifurcation("tf", 0.010, 0.010, 0.004, 0.004),
            system:    create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &geom_4mm,
                &ChannelTypeConfig::AllStraight,
            ),
            config_override: None,
        },

        // ── Serpentine chain: winding single channel ───────────────────────────
        //
        // Winding single-channel design: 3 serpentine legs connected by U-bends.
        // Useful for maximising residence time in a compact footprint.
        // Both ports are on the left face at different Y positions.
        Design {
            name: "serpentine_chain",
            description:
                "Serpentine path: inlet → 3 winding rows via U-bends → outlet  \
                 (Ø 4 mm, both ports on the left face at different Y positions)",
            blueprint: serpentine_chain("sc", 3, 0.010, 0.004),
            system:    create_geometry(
                box_dims,
                &[],
                &geom_4mm,
                &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
            ),
            config_override: None,
        },

        // ── Venturi rect: rectangular cross-section taper ─────────────────────
        //
        // Venturi with rectangular cross-section channels.  Useful as an
        // alternative to the circular venturi_chain for rectangular duct CFD
        // simulations.  Inlet and outlet ports on opposite chip faces.
        Design {
            name: "venturi_rect",
            description:
                "Rectangular duct Venturi: inlet 4×5 mm → 2×5 mm throat → outlet  \
                 (one rectangular port on each end face)",
            blueprint: {
                use cfd_schematics::interface::presets::venturi_rect;
                venturi_rect("vr", 0.004, 0.002, 0.004, 0.005)
            },
            system:    create_geometry(
                box_dims,
                &[],
                &geom_4mm,
                &ChannelTypeConfig::AllFrustum(FrustumConfig {
                    inlet_width:  4.0,
                    throat_width: 2.0,
                    outlet_width: 4.0,
                    ..FrustumConfig::default()
                }),
            ),
            config_override: None,
        },
    ]
}
