//! Turbulence modeling for 2D CFD simulations
//!
//! Implements various turbulence models including:
//! - k-ε model
//! - k-ω SST model
//! - Spalart-Allmaras model
//! - Wall functions for near-wall treatment

pub mod boundary_conditions;
pub mod constants;
pub mod k_epsilon;
pub mod k_omega_sst;
pub mod spalart_allmaras;
pub mod traits;
pub mod validation;
pub mod wall_functions;

pub use boundary_conditions::{
    TurbulenceBoundaryCondition, TurbulenceBoundaryManager,
};
pub use constants::*;
pub use k_epsilon::KEpsilonModel;
pub use k_omega_sst::KOmegaSSTModel;
pub use spalart_allmaras::SpalartAllmaras;
pub use traits::TurbulenceModel;
pub use validation::{run_turbulence_validation, TurbulenceValidator, ValidationResult};
pub use wall_functions::{WallFunction, WallTreatment};

// Literature-based validation tests
#[cfg(test)]
mod literature_validation_tests;

// Comprehensive k-ω SST tests
#[cfg(test)]
#[path = "k_omega_sst_tests.rs"]
mod k_omega_sst_tests;
