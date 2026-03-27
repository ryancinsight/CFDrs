//! Core Delaunay triangulation kernel.
//!
//! The triangulation is stored as a flat array of [`Triangle`] structures with
//! explicit triangle-triangle adjacency, following the indexing conventions of
//! Guibas & Stolfi (1985) and Shewchuk (1996).

pub mod adjacency;
pub mod bowyer_watson;
pub mod locate;
pub(crate) mod ordering;
pub mod triangle;
pub(crate) mod vertex_star;

pub use adjacency::Adjacency;
pub use bowyer_watson::DelaunayTriangulation;
pub use triangle::{Triangle, TriangleId, GHOST_TRIANGLE};
