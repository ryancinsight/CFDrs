//! High-order spatial discretization schemes
//!
//! This module provides advanced spatial discretization methods for CFD,
//! including shock-capturing schemes and high-order accurate reconstructions.

pub use crate::high_order::weno::{WENO5, WENO5 as Weno5, WENO7, WenoReconstruction};
