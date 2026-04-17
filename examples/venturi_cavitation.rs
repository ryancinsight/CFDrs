//! Selective venturi cavitation screening example.

use cfd_1d::{assess_venturi_screening, evaluate_venturi_screening, VenturiScreeningInput};
use cfd_core::physics::cavitation::{
    CellMechanicalState, CellPopulationIdentity, PopulationNucleationState,
    SelectiveCavitationInput, SelectiveCavitationPopulation,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let selective = SelectiveCavitationInput {
        base_vapor_pressure_pa: 3_170.0,
        ambient_pressure_pa: 101_325.0,
        density_kg_m3: 1_025.0,
        populations: vec![
            SelectiveCavitationPopulation {
                identity: CellPopulationIdentity::HealthyRbc,
                label: "Healthy RBC".to_string(),
                mechanical_state: CellMechanicalState {
                    membrane_stiffness_pa: 72_000.0,
                    interfacial_tension_n_m: 0.052,
                    particle_radius_m: 3.3e-6,
                    deformability_factor: 1.0,
                },
                nucleation_state: PopulationNucleationState {
                    volume_fraction: 0.88,
                    upstream_nuclei_fraction: 0.01,
                    seed_density_factor: 0.10,
                    inception_weight: 1.0,
                },
            },
            SelectiveCavitationPopulation {
                identity: CellPopulationIdentity::CirculatingTumorCell,
                label: "CTC".to_string(),
                mechanical_state: CellMechanicalState {
                    membrane_stiffness_pa: 28_000.0,
                    interfacial_tension_n_m: 0.036,
                    particle_radius_m: 8.0e-6,
                    deformability_factor: 1.35,
                },
                nucleation_state: PopulationNucleationState {
                    volume_fraction: 0.02,
                    upstream_nuclei_fraction: 0.02,
                    seed_density_factor: 0.35,
                    inception_weight: 1.0,
                },
            },
            SelectiveCavitationPopulation {
                identity: CellPopulationIdentity::HealthyWbc,
                label: "Healthy WBC".to_string(),
                mechanical_state: CellMechanicalState {
                    membrane_stiffness_pa: 45_000.0,
                    interfacial_tension_n_m: 0.044,
                    particle_radius_m: 5.0e-6,
                    deformability_factor: 1.1,
                },
                nucleation_state: PopulationNucleationState {
                    volume_fraction: 0.10,
                    upstream_nuclei_fraction: 0.015,
                    seed_density_factor: 0.20,
                    inception_weight: 1.0,
                },
            },
        ],
    };

    let screening = evaluate_venturi_screening(VenturiScreeningInput {
        upstream_pressure_pa: 150_000.0,
        upstream_velocity_m_s: 0.30,
        throat_velocity_m_s: 3.10,
        throat_hydraulic_diameter_m: 7.5e-4,
        throat_length_m: 2.5e-3,
        density_kg_m3: 1_025.0,
        viscosity_pa_s: 1.2e-3,
        vapor_pressure_pa: 3_170.0,
        vena_contracta_coeff: 0.92,
        diffuser_recovery_coeff: 0.64,
        upstream_nuclei_fraction: 0.02,
        selective_cavitation: Some(selective),
    })?;
    let assessment = assess_venturi_screening(&screening);

    println!("Throat pressure [Pa]: {:.3}", screening.throat_static_pressure_pa);
    println!(
        "Mixture inception threshold [Pa]: {:.3}",
        screening.mixture_inception_threshold_pa
    );
    println!("Selectivity margin [Pa]: {:.3}", screening.selectivity_margin_pa);
    println!("Hydrodynamic risk: {}", assessment.risk);
    println!("Selective regime: {:?}", screening.screening_regime);
    if let Some(label) = &screening.dominant_selective_label {
        println!("Dominant population: {label}");
    }
    for threshold in &screening.population_thresholds_pa {
        println!(
            "  {} -> effective threshold {:.3} Pa",
            threshold.label, threshold.effective_threshold_pressure_pa
        );
    }

    Ok(())
}
