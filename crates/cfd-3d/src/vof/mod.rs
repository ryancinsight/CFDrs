//! Volume of Fluid (VOF) method for interface tracking in 3D multiphase flows
//!
//! The VOF method tracks interfaces by advecting volume fractions,
//! providing excellent mass conservation properties.

mod config;
mod solver;
mod reconstruction;
mod advection;
mod initialization;

pub use config::{VofConfig, constants};
pub use solver::VofSolver;

// Re-export key types for convenience
pub use reconstruction::InterfaceReconstruction;
pub use advection::AdvectionMethod;