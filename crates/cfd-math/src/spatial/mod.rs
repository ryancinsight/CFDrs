//! High-order spatial discretization schemes
//!
//! This module provides advanced spatial discretization methods for CFD,
//! including shock-capturing schemes and high-order accurate reconstructions.

pub mod weno;

// Re-export main types
pub use weno::{Weno5, Weno7, WenoConfig, WenoReconstruction};
