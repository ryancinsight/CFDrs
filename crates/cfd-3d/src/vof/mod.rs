//! Volume of Fluid (VOF) method for interface tracking in 3D multiphase flows
//!
//! The VOF method tracks interfaces by advecting volume fractions,
//! providing excellent mass conservation properties.

mod advection;
mod config;
mod initialization;
mod reconstruction;
mod solver;
pub use config::{constants, VofConfig};
pub use solver::VofSolver;
// Re-export key types for convenience
pub use advection::AdvectionMethod;
pub use reconstruction::InterfaceReconstruction;
