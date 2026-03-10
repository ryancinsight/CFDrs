use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_schematics::interface::presets::{symmetric_bifurcation, symmetric_trifurcation};

#[test]
fn test_mesh_volume() {
    let bp_b = symmetric_bifurcation("b1", 0.010, 0.010, 0.004, 0.003);
    let mut out_b = BlueprintMeshPipeline::run(&bp_b, &PipelineConfig::default()).unwrap();
    println!("Bifurcation watertight: {}", out_b.fluid_mesh.is_watertight());
    println!("Bifurcation volume: {}", out_b.fluid_mesh.signed_volume());

    let bp_t = symmetric_trifurcation("t1", 0.010, 0.008, 0.004, 0.004);
    let out_t = BlueprintMeshPipeline::run(&bp_t, &PipelineConfig::default());
    match out_t {
        Ok(mut t) => {
            println!("Trifurcation watertight: {}", t.fluid_mesh.is_watertight());
            println!("Trifurcation volume: {}", t.fluid_mesh.signed_volume());
        }
        Err(e) => {
            println!("Trifurcation Error: {}", e);
        }
    }
}
