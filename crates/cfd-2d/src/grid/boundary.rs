//! Boundary types and handling for grids
//!
//! # Theorem
//! The grid topology must form a valid, non-overlapping partition of the computational domain.
//!
//! **Proof sketch**:
//! For a finite volume discretization to be conservative, the control volumes $\Omega_i$
//! must satisfy $\cup_i \Omega_i = \Omega$ and $\Omega_i \cap \Omega_j = \emptyset$ for $i \neq j$.
//! The grid data structures enforce this by maintaining strict adjacency invariants
//! and ensuring that the sum of face area vectors for any closed cell is exactly zero:
//! $\sum_f \mathbf{A}_f = \mathbf{0}$.

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
