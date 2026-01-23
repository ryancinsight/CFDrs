//! High-order spatial discretization schemes
//!
//! This module provides advanced spatial discretization methods for CFD,
//! including shock-capturing schemes and high-order accurate reconstructions.

// TODO(MEDIUM): Eliminate WENO duplication - consolidate high_order/weno.rs and spatial re-export into single module
// Dependencies: Architectural restructuring (REST-001 to REST-005)
// Reference: Shu (1997) WENO schemes - ensure single implementation for consistency

pub use crate::high_order::weno::{WenoReconstruction, WENO5, WENO5 as Weno5, WENO7};
