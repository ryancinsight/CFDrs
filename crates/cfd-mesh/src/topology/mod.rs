//! Mesh topology components: vertices, edges, faces, and cells
//!
//! This module provides the fundamental building blocks for mesh representation:
//!
//! - [`Vertex`]: A point in 3D space
//! - [`Edge`]: A connection between two vertices
//! - [`Face`]: A 2D polygon (2-simplex) defined by vertices
//! - [`Cell`]: A 3D volume defined by faces
//!
//! # Topological Invariants
//!
//! - Edges require exactly 2 distinct vertices
//! - Faces require at least 3 distinct vertices
//! - Cells are defined by face indices
//!
//! # Distributed Mesh Support
//!
//! [`Vertex`], [`Face`], and [`Cell`] all support distributed mesh properties
//! (`global_id`, `partition_id`) for MPI-based mesh partitioning.

mod cell;
mod edge;
mod face;
mod vertex;

pub use cell::Cell;
pub use edge::Edge;
pub use face::{Face, MIN_FACE_VERTICES};
pub use vertex::Vertex;

// Re-export ElementType for convenience
pub use cfd_core::domain::ElementType;
