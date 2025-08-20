//! Pressure-velocity coupling using Semi-Implicit Method for Pressure-Linked Equations
//!
//! Modularized implementation following SLAP and SOC principles.
//! Reference: Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow"

pub mod config;
pub mod coefficients;
pub mod pressure;
pub mod solver;
pub mod rhie_chow;

pub use config::PressureVelocityConfig;
pub use coefficients::CellCoefficients;
pub use solver::PressureVelocitySolver;
pub use pressure::PressureCorrectionSolver;
pub use rhie_chow::RhieChowInterpolation;

