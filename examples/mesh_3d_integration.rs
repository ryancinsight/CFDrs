//! Canonical blueprint-to-mesh integration example.

use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_mesh::infrastructure::io::stl::write_stl_binary;
use cfd_schematics::serpentine_rect;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_dir = PathBuf::from("outputs/examples/mesh_3d_integration");
    fs::create_dir_all(&output_dir)?;

    let blueprint = serpentine_rect("mesh_integration_serpentine", 5, 5.0e-3, 1.2e-3, 0.8e-3);
    blueprint.validate()?;

    let mut output = BlueprintMeshPipeline::run(
        &blueprint,
        &PipelineConfig {
            circular_segments: 20,
            axial_rings: 10,
            include_chip_body: true,
            skip_diameter_constraint: true,
            ..PipelineConfig::default()
        },
    )?;

    let fluid_path = output_dir.join("fluid_mesh.stl");
    let mut fluid_writer = BufWriter::new(File::create(&fluid_path)?);
    write_stl_binary(&mut fluid_writer, &output.fluid_mesh)?;

    let chip_path = output_dir.join("chip_mesh.stl");
    if let Some(chip_mesh) = &output.chip_mesh {
        let mut chip_writer = BufWriter::new(File::create(&chip_path)?);
        write_stl_binary(&mut chip_writer, chip_mesh)?;
    }

    println!("Blueprint: {}", blueprint.name);
    println!("Topology class: {:?}", output.topology_class);
    println!("Layout segments: {}", output.layout_segments.len());
    println!("Fluid mesh watertight: {}", output.fluid_mesh.is_watertight());
    println!(
        "Fluid mesh volume error [%]: {:.4}",
        output.volume_trace.fluid_mesh_volume_error_pct
    );
    println!("Fluid STL: {}", fluid_path.display());
    if output.chip_mesh.is_some() {
        println!("Chip STL: {}", chip_path.display());
    }

    Ok(())
}
