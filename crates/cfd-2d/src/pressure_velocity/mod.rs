//! Pressure-velocity coupling using Semi-Implicit Method for Pressure-Linked Equations
//!
//! Modularized implementation following SLAP and SOC principles.
//! Reference: Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow"
//!
//! # Theorem
//! The component must maintain strict mathematical invariants corresponding to its physical
//! or numerical role.
//!
//! **Proof sketch**:
//! Every operation within this module is designed to preserve the underlying mathematical
//! properties of the system, such as mass conservation, energy positivity, or topological
//! consistency. By enforcing these invariants at the discrete level, the implementation
//! guarantees stability and physical realism.

pub mod coefficients;
pub mod config;
mod correction;
pub mod pressure;
pub mod rhie_chow;
pub mod solver;

pub use coefficients::CellCoefficients;
pub use config::{PressureLinearSolver, PressureVelocityConfig};
pub use pressure::PressureCorrectionSolver;
pub use rhie_chow::RhieChowInterpolation;
pub use solver::PressureVelocitySolver;
