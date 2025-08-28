//! Material properties domain module
//!
//! Provides fluid and material property models following Domain-Driven Design

pub mod database;
pub mod fluids;
pub mod interfaces;
pub mod non_newtonian;
pub mod property_calculators;
pub mod service;
pub mod solids;
pub mod traits;

// Re-export main types
pub use database::MaterialDatabase;
pub use fluids::NewtonianFluid;
pub use interfaces::{FluidSolidInterface, WettingProperties};
pub use non_newtonian::{BinghamFluid, PowerLawFluid};
pub use property_calculators::{
    KinematicViscosityCalculator, PrandtlNumberCalculator, PropertyCalculator,
    ReynoldsNumberCalculator,
};
pub use service::MaterialPropertiesService;
pub use solids::ElasticSolid;
pub use traits::{FluidProperties, InterfaceProperties, SolidProperties};
