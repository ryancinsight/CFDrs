//! 2D Navier-Stokes FVM solver hierarchy.
//!
//! ## Structure
//! ```text
//! ns_fvm/
//!   mod.rs        — canonical re-exports
//!   grid.rs       — StaggeredGrid2D
//!   field.rs      — FlowField2D
//!   boundary.rs   — BCType, BoundaryCondition, BloodModel
//!   config.rs     — SIMPLEConfig, SolveResult
//!   solver.rs     — NavierStokesSolver2D + SIMPLE implementation
//! ```

pub mod boundary;
pub mod config;
pub mod field;
pub mod grid;
pub mod solver;

pub use boundary::{BCType, BloodModel, BoundaryCondition};
pub use config::{SIMPLEConfig, SolveResult};
pub use field::FlowField2D;
pub use grid::StaggeredGrid2D;
pub use solver::NavierStokesSolver2D;
