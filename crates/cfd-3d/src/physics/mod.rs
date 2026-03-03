//! Physics models for 3D CFD
//!
//! # Theorem — Eddy Viscosity Positivity
//!
//! For closure models implemented in this module, turbulent viscosity is
//! constrained non-negative ($\nu_t \ge 0$), ensuring added dissipation and
//! preventing anti-diffusive instability in momentum equations.

#[path = "turbulence/mod.rs"]
pub mod turbulence;
