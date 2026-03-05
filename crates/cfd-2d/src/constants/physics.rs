//! Physical constants for 2D CFD simulations
//!
//! # Invariant
//!
//! All constants are defined as exact values from their respective physical or
//! numerical derivations (e.g., $c_s^2 = 1/3$ for D2Q9 LBM follows from the
//! second-order isotropy condition of the lattice).

/// Speed of sound squared for lattice Boltzmann (cs^2)
pub const LATTICE_SOUND_SPEED_SQUARED: f64 = 1.0 / 3.0;

/// Central difference coefficient
pub const CENTRAL_DIFF_COEFF: f64 = 2.0;
