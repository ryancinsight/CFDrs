//! Constants module for CFD simulations
//!
//! Comprehensive collection of physical, numerical, and mathematical constants
//! organized by domain to ensure SSOT (Single Source of Truth).

pub mod physics;
pub mod numerical;

// Re-export commonly used constants at module level for convenience
pub use physics::fluid;
pub use physics::thermo;
pub use physics::turbulence;
pub use physics::dimensionless;
pub use numerical::solver;
pub use numerical::discretization;