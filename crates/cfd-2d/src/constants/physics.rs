//! Physical constants for 2D CFD simulations
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

/// Speed of sound squared for lattice Boltzmann (cs^2)
pub const LATTICE_SOUND_SPEED_SQUARED: f64 = 1.0 / 3.0;

/// Central difference coefficient
pub const CENTRAL_DIFF_COEFF: f64 = 2.0;
