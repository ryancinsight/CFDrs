//! Domain-specific modules organized by CFD domain concerns.
//!
//! This module organizes CFD functionality by domain following Domain-Driven Design principles:
//! - Fluid Dynamics: Core fluid mechanics concepts
//! - Numerical Methods: Mathematical algorithms and solvers
//! - Mesh Operations: Geometry and discretization
//! - Boundary Conditions: Physical constraints and conditions
//! - Material Properties: Physical properties and constitutive relations

pub mod fluid_dynamics;
pub mod numerical_methods;
pub mod mesh_operations;
pub mod boundary_conditions;
pub mod material_properties;

// Re-export domain-specific functionality
pub use fluid_dynamics::{FlowField, VelocityField, PressureField, TurbulenceModel};
pub use numerical_methods::{DiscretizationScheme, TimeIntegrationScheme, LinearSystemSolver};
pub use mesh_operations::{MeshGeneration, MeshRefinement, MeshQuality};
pub use boundary_conditions::{BoundaryConditionType, BoundaryConditionApplicator};
pub use material_properties::{FluidProperties, SolidProperties, InterfaceProperties};
