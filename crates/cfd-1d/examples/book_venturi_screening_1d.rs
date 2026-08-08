//! Book example: compact 1D Venturi selective-screening workflow.

use cfd_1d::{
    assess_venturi_screening, evaluate_venturi_screening, venturi_taper_length_m,
    VenturiScreeningInput,
};
use aequitas::systems::si::quantities::{
    DynamicViscosity, Length, MassDensity, Pressure, Velocity,
};
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("book_venturi_screening_1d");
    fs::create_dir_all(&output_dir)?;

    let inlet_width_m = 2.0e-3;
    let throat_width_m = 8.0e-4;
    let convergent_half_angle_deg = 12.0;
    let taper_length_m =
        venturi_taper_length_m(inlet_width_m, throat_width_m, convergent_half_angle_deg)?;

    let input = VenturiScreeningInput {
        upstream_pressure_pa: Pressure::from_base(140_000.0),
        upstream_velocity_m_s: Velocity::from_base(0.28),
        throat_velocity_m_s: Velocity::from_base(0.28 * (inlet_width_m / throat_width_m).powi(2)),
        throat_hydraulic_diameter_m: Length::from_base(throat_width_m),
        throat_length_m: Length::from_base(taper_length_m),
        density_kg_m3: MassDensity::from_base(1_000.0),
        viscosity_pa_s: DynamicViscosity::from_base(1.0e-3),
        vapor_pressure_pa: Pressure::from_base(2_300.0),
        vena_contracta_coeff: 0.97,
        diffuser_recovery_coeff: 0.60,
        upstream_nuclei_fraction: 0.004,
        selective_cavitation: None,
    };

    let upstream_velocity = input.upstream_velocity_m_s.into_base();
    let result = evaluate_venturi_screening(input)?;
    let assessment = assess_venturi_screening(&result);

    assert!(result.cavitation_number.is_finite(), "cavitation number must be finite");
    assert!(result.effective_throat_velocity_m_s.into_base() > upstream_velocity);
    assert!(result.bernoulli_drop_pa.into_base() > 0.0, "Bernoulli drop must be positive");
    assert!(result.throat_static_pressure_pa.into_base() > 0.0, "throat pressure must stay physical");

    let summary = format!(
        "sigma={:.6}\nvelocity_up={:.6}\nvelocity_throat={:.6}\npressure_drop_pa={:.6}\nregime={:?}\nrisk={:?}\n",
        result.cavitation_number,
        upstream_velocity,
        result.effective_throat_velocity_m_s.into_base(),
        result.bernoulli_drop_pa.into_base(),
        result.screening_regime,
        assessment.risk
    );
    fs::write(output_dir.join("summary.txt"), summary)?;

    Ok(())
}
