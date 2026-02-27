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

// Re-export core types for simplified access
pub use self::boundary::{
    BoundaryCondition, BoundaryConditionManager, BoundaryConditionSet, WallType,
};
pub use self::cavitation::{
    CavitationDamage, CavitationModel, CavitationNumber, CavitationRegime,
    CavitationRegimeAnalysis, CavitationRegimeClassifier, RayleighPlesset, VenturiCavitation,
};
pub use self::constants::{mathematical, physics};
pub use self::fluid::{ConstantPropertyFluid, Fluid, FluidProperties};
pub use self::fluid_dynamics::{
    FlowClassifier, FlowField, FlowOperations, FlowRegime, FluidDynamicsService, PressureField,
    RANSModel, RhieChowInterpolation, ScalarField, TurbulenceModel, VelocityField,
};
pub use self::hemolysis::{
    BloodTrauma, BloodTraumaSeverity, HemolysisCalculator, HemolysisModel, PlateletActivation,
};
pub use self::material::{MaterialDatabase, SolidProperties};
pub use self::values::{Pressure, ReynoldsNumber, Temperature, Velocity};
