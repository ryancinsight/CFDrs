//! Constrained Delaunay Triangulation (CDT) extension.
//!
//! After constructing a Delaunay triangulation, constraint segments from the
//! PSLG are enforced by flipping non-Delaunay edges that cross them.

pub mod cavity;
pub mod enforce;

pub use enforce::Cdt;
