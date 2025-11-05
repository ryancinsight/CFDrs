//! Turbulence modeling for 2D CFD simulations
//!
//! Implements various turbulence models including:
//! - k-ε model (RANS)
//! - k-ω SST model (RANS)
//! - Reynolds Stress Transport Model (RSTM) (RANS)
//! - Spalart-Allmaras model (RANS)
//! - Smagorinsky LES model (LES)
//! - Detached Eddy Simulation (DES)
//! - Wall functions for near-wall treatment

pub mod boundary_conditions;
pub mod constants;
pub mod constants_validation;
pub mod des;
pub mod k_epsilon;
pub mod k_omega_sst;
pub mod les_smagorinsky;
pub mod reynolds_stress;
pub mod spalart_allmaras;
pub mod traits;
pub mod validation;
pub mod wall_functions;

pub use boundary_conditions::{
    TurbulenceBoundaryCondition, TurbulenceBoundaryManager,
};
pub use constants::*;
pub use constants_validation::{run_turbulence_constants_validation, TurbulenceConstantsValidator, ConstantsValidationResult};
pub use des::DetachedEddySimulation;
pub use k_epsilon::KEpsilonModel;
pub use k_omega_sst::KOmegaSSTModel;
pub use les_smagorinsky::SmagorinskyLES;
pub use reynolds_stress::{ReynoldsStressModel, ReynoldsStressTensor, PressureStrainModel};
pub use spalart_allmaras::SpalartAllmaras;
pub use traits::{TurbulenceModel, LESTurbulenceModel};
pub use validation::{run_turbulence_validation, run_rans_benchmark_suite, run_les_benchmark_suite, TurbulenceValidator, ValidationResult};
pub use wall_functions::{WallFunction, WallTreatment};

// Literature-based validation tests
#[cfg(test)]
mod literature_validation_tests;

// Comprehensive k-ω SST tests
#[cfg(test)]
#[path = "k_omega_sst_tests.rs"]
mod k_omega_sst_tests;
