//! Millifluidic Shell TPMS export — schematic JSON → binary STL.
//!
//! Demonstrates the creation of a millifluidic chip modelled as a shell
//! enclosing a cavity filled with a TPMS lattice (Gyroid in this case).
//! Uses the `ShellMeshPipeline` to produce:
//! - **fluid mesh** — the watertight internal cavity filled with TPMS + ports.
//! - **chip body** — the SBS-96 cuboid substrate with the fluid void subtracted.
//!
//! Outputs STLs to `outputs/millifluidic_chip_stl/shell_tpms_gyroid/`.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --features scheme-io --example millifluidic_shell_tpms
//! ```

use std::fs;
use std::io::BufWriter;
use std::path::Path;

use cfd_mesh::application::pipeline::{ShellMeshPipeline, ShellPipelineConfig};
use cfd_mesh::infrastructure::io::stl;

use cfd_schematics::geometry::types::{
    InterchangeShellCuboid, InterchangeShellPort, TpmsFillSpec, TpmsSurfaceKind,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════╗");
    println!("║       Millifluidic Shell TPMS Example               ║");
    println!("╚══════════════════════════════════════════════════════╝");
    println!();

    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("millifluidic_chip_stl")
        .join("shell_tpms_gyroid");

    fs::create_dir_all(&out_dir)?;

    let name = "shell_tpms_gyroid";

    // Build the schematic defining the shell, cavity, and ports.
    let shell = InterchangeShellCuboid {
        schema_version: "1.0.0".to_string(),
        producer: "cfd-mesh-example".to_string(),
        length_units: "mm".to_string(),
        // SBS-96 microplate footprint
        outer_dims_mm: (127.76, 85.47),
        // 2mm walls
        inner_dims_mm: (123.76, 81.47),
        shell_thickness_mm: 2.0,
        ports: vec![
            InterchangeShellPort {
                label: "inlet".to_string(),
                outer_point_mm: (0.0, 42.735),
                inner_point_mm: (2.0, 42.735),
            },
            InterchangeShellPort {
                label: "outlet".to_string(),
                outer_point_mm: (127.76, 42.735),
                inner_point_mm: (125.76, 42.735),
            },
        ],
        // Fill the cavity with a coarse Gyroid
        tpms_fill: Some(TpmsFillSpec {
            surface: TpmsSurfaceKind::Gyroid,
            period_mm: 6.0,
            iso_value: 0.0,
            resolution: 48,
            gradient: None,
        }),
    };

    println!("Building {} ...", name);

    let config = ShellPipelineConfig::default();
    
    // Run the pipeline
    let mut output = ShellMeshPipeline::run(&shell, &config)
        .map_err(|e| format!("Pipeline failed: {}", e))?;

    assert!(output.fluid_mesh.is_watertight(), "Fluid mesh not watertight!");
    assert!(output.chip_body_mesh.is_watertight(), "Chip body mesh not watertight!");

    println!(
        "  fluid : {:>7} faces, vol = {:>10.3} mm³ [watertight ✓]",
        output.fluid_mesh.face_count(),
        output.fluid_mesh.signed_volume()
    );
    println!(
        "  chip  : {:>7} faces, vol = {:>10.3} mm³ [watertight ✓]",
        output.chip_body_mesh.face_count(),
        output.chip_body_mesh.signed_volume()
    );

    // Write STL
    let fluid_path = out_dir.join(format!("{}_fluid.stl", name));
    let chip_path = out_dir.join(format!("{}_chip.stl", name));

    stl::write_binary_stl(
        &mut BufWriter::new(fs::File::create(&fluid_path)?),
        &output.fluid_mesh.vertices,
        &output.fluid_mesh.faces,
    )?;
    println!("  → {}", fluid_path.display());

    stl::write_binary_stl(
        &mut BufWriter::new(fs::File::create(&chip_path)?),
        &output.chip_body_mesh.vertices,
        &output.chip_body_mesh.faces,
    )?;
    println!("  → {}", chip_path.display());

    println!("Done.");
    Ok(())
}
