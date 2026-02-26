//! 2D grid structures for CFD simulations.
//!
//! This module provides various grid types for 2D CFD simulations, including
//! structured and unstructured grids with support for adaptive refinement.
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

pub mod boundary;
pub mod refinement;
pub mod structured;
pub mod traits;
pub mod unstructured;

pub use boundary::BoundaryType;
pub use refinement::{AdaptiveGrid2D, RefinementCriterion};
pub use structured::StructuredGrid2D;
pub use traits::Grid2D;
pub use unstructured::UnstructuredGrid2D;
