//! 2D grid structures for CFD simulations.
//!
//! This module provides various grid types for 2D CFD simulations, including
//! structured and unstructured grids with support for adaptive refinement.

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
