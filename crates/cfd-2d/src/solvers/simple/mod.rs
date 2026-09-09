//! SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
//!
//! Implements the SIMPLE algorithm for pressure-velocity coupling in
//! incompressible CFD on collocated grids with Rhie-Chow interpolation.
//!
//! # Theorem 1 — SIMPLE Convergence (Patankar 1980)
//! (See mathematical proofs in component files)

/// Coordinated iterative algorithm loop.
pub mod algorithm;
/// Predictor momentum definitions.
pub mod momentum;
/// Predictor pressure-correction definitions.
pub mod pressure;

pub use algorithm::SimpleAlgorithm;

/// Minimum A_P threshold below which a cell is treated as stagnant (D = 0,
/// no pressure-velocity correction applied).
pub const STAGNANT_CELL_AP_THRESHOLD: f64 = 1e-10;
