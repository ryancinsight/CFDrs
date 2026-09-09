//! Cavitation physics models for CFD simulations
//!
//! This module implements hydrodynamic cavitation models based on established
//! literature including Brennen (1995), Franc & Michel (2004), and Rayleigh-Plesset
//! dynamics for bubble growth and collapse.

pub mod bio_damage;
pub mod constants;
pub mod damage;
pub mod heterogeneous_nucleation;
pub mod models;
pub mod nuclei_transport;
pub mod number;
pub mod rayleigh_plesset;
pub mod regimes;
pub mod venturi;

// Re-export main types
pub use bio_damage::{CellularInjuryProfile, CellularMembraneMechanics};
pub use damage::CavitationDamage;
pub use heterogeneous_nucleation::{
    evaluate_selective_cavitation_thresholds, heterogeneous_inception_threshold_pa,
    rank_population_selectivity, validate_selective_cavitation_input, CellMechanicalState,
    CellPopulationIdentity, CellularPopulation, PopulationCavitationThreshold,
    PopulationNucleationState, SelectiveCavitationInput, SelectiveCavitationPopulation,
    SelectiveCavitationResult,
};
pub use models::CavitationModel;
pub use nuclei_transport::{NucleiTransport, NucleiTransportConfig};
pub use number::CavitationNumber;
pub use rayleigh_plesset::RayleighPlesset;
pub use regimes::{CavitationRegime, CavitationRegimeAnalysis, CavitationRegimeClassifier};
pub use venturi::VenturiCavitation;
