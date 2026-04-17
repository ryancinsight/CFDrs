//! Blueprint-native mesh plus reduced-order venturi screening example.
//!
//! The historical example name is retained for compatibility. The executable
//! workflow now uses the canonical `cfd-schematics -> gaia -> cfd-1d` path
//! instead of an unimplemented ad hoc CSG stub.

use cfd_1d::{assess_venturi_screening, evaluate_venturi_screening, VenturiScreeningInput};
use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_schematics::venturi_rect;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inlet_width_m = 2.0e-3;
    let throat_width_m = 0.7e-3;
    let height_m = 1.0e-3;
    let throat_length_m = 2.4e-3;

    let blueprint = venturi_rect(
        "csg_cfd_simulation_reference",
        inlet_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
    );
    blueprint.validate()?;

    let mut mesh = BlueprintMeshPipeline::run(
        &blueprint,
        &PipelineConfig {
            circular_segments: 24,
            axial_rings: 12,
            include_chip_body: false,
            skip_diameter_constraint: true,
            ..PipelineConfig::default()
        },
    )?;

    let screening = evaluate_venturi_screening(VenturiScreeningInput {
        upstream_pressure_pa: 160_000.0,
        upstream_velocity_m_s: 0.35,
        throat_velocity_m_s: 2.85,
        throat_hydraulic_diameter_m: 2.0 * throat_width_m * height_m / (throat_width_m + height_m),
        throat_length_m,
        density_kg_m3: 1_025.0,
        viscosity_pa_s: 1.2e-3,
        vapor_pressure_pa: 3_170.0,
        vena_contracta_coeff: 0.91,
        diffuser_recovery_coeff: 0.62,
        upstream_nuclei_fraction: 0.015,
        selective_cavitation: None,
    })?;
    let assessment = assess_venturi_screening(&screening);

    println!("Blueprint: {}", blueprint.name);
    println!("Fluid mesh watertight: {}", mesh.fluid_mesh.is_watertight());
    println!("Fluid mesh faces: {}", mesh.fluid_mesh.face_count());
    println!("Fluid mesh signed volume [mm^3]: {:.6}", mesh.fluid_mesh.signed_volume());
    println!(
        "Mesh volume contract error [%]: {:.4}",
        mesh.volume_trace.fluid_mesh_volume_error_pct
    );
    println!(
        "Venturi throat pressure [Pa]: {:.3}",
        screening.throat_static_pressure_pa
    );
    println!(
        "Effective vapor threshold [Pa]: {:.3}",
        screening.effective_vapor_pressure_pa
    );
    println!("Cavitation margin [Pa]: {:.3}", assessment.cavitation_margin_pa);
    println!("Screening risk: {}", assessment.risk);

    Ok(())
}
