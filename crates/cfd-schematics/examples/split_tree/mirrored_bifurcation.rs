use cfd_schematics::interface::presets::symmetric_bifurcation;

#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    tracing::info!("Mirrored bifurcation demo");
    tracing::info!("------------------------");

    let system = symmetric_bifurcation("mirrored_bifurcation", 30.0e-3, 50.0e-3, 2.0e-3, 2.0e-3);
    shared::save_example_output(&system, "mirrored_bifurcation");

    Ok(())
}
