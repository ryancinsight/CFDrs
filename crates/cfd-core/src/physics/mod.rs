//! Physics module consolidating fluid properties, boundary conditions, and constants.
//!
//! This module provides a unified interface to the physical models and properties
//! used in CFD simulations, following a deep vertical module structure.

pub mod boundary;
pub mod cavitation;
pub mod constants;
pub mod fluid;
pub mod fluid_dynamics;
pub mod hemolysis;
pub mod material;
pub mod values;

// Re-exports for convenience

// Boundary Conditions
pub use boundary::{
    BoundaryCondition, BoundaryConditionManager, BoundaryConditionSet, WallType,
};

// Fluid Properties
pub use fluid::{ConstantPropertyFluid, Fluid, FluidProperties};

// Fluid Dynamics (Fields & Models)
pub use fluid_dynamics::{
    FlowClassifier, FlowField, FlowOperations, FlowRegime, FluidDynamicsService, PressureField,
    RANSModel, RhieChowInterpolation, ScalarField, TurbulenceModel, VelocityField,
};

// Biological Models
pub use hemolysis::{
    BloodTrauma, BloodTraumaSeverity, HemolysisCalculator, HemolysisModel, PlateletActivation,
};

// Cavitation Models
pub use cavitation::{
    CavitationDamage, CavitationModel, CavitationNumber, CavitationRegime, RayleighPlesset,
    VenturiCavitation,
};

// Material Properties
pub use material::{MaterialDatabase, SolidProperties};

// Physical Values
pub use values::{Pressure, ReynoldsNumber, Temperature, Velocity};

// Constants
pub use constants::{mathematical, physics};
