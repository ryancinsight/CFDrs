//! Domain-specific modules organized by CFD domain concerns.
//!
//! This module organizes CFD functionality by domain following Domain-Driven Design principles:
//! - Fluid Dynamics: Core fluid mechanics concepts
//! - Numerical Methods: Mathematical algorithms and solvers
//! - Mesh Operations: Geometry and discretization
//! - Boundary Conditions: Physical constraints and conditions
//! - Material Properties: Physical properties and constitutive relations

pub mod boundary_conditions;
pub mod fluid_dynamics;
pub mod material_properties;
pub mod mesh_operations;
pub mod numerical_methods;
// Re-export domain-specific functionality
pub use boundary_conditions::BoundaryConditionApplicator;
pub use fluid_dynamics::{
    FlowClassifier, FlowField, FlowOperations, FlowRegime, KEpsilonConstants, KEpsilonModel,
    MixingLengthModel, PressureField, RANSModel, ScalarField, SmagorinskyModel, TurbulenceModel,
    VelocityField,
};
pub use material_properties::{FluidProperties, InterfaceProperties, SolidProperties};
pub use mesh_operations::{MeshGeneration, MeshQuality, MeshRefinement};
pub use numerical_methods::{DiscretizationScheme, LinearSystemSolver, TimeIntegrationScheme};
