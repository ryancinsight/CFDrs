//! SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
//!
//! Modularized implementation following SLAP and SOC principles.
//! Reference: Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow"

pub mod config;
pub mod coefficients;
pub mod momentum;
pub mod pressure;
pub mod solver;
pub mod rhie_chow;

pub use config::SIMPLEConfig;
pub use coefficients::CellCoefficients;
pub use solver::SIMPLESolver;
pub use momentum::MomentumSolver;
pub use pressure::PressureCorrectionSolver;
pub use rhie_chow::RhieChowInterpolation;

// Re-export for backward compatibility (temporary)
pub use solver::SIMPLESolver as PressureVelocityCouplerSolver;
pub use config::SIMPLEConfig as PressureVelocityCouplingConfig;