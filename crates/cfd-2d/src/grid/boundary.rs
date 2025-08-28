//! Boundary types and handling for grids

use serde::{Deserialize, Serialize};

/// Boundary types for grid cells
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryType {
    /// Wall boundary (no-slip)
    Wall,
    /// Inlet boundary
    Inlet,
    /// Outlet boundary
    Outlet,
    /// Symmetry boundary
    Symmetry,
    /// Periodic boundary
    Periodic,
}
