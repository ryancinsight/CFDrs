//! PSLG vertex, segment, and graph types.
//!
//! A Planar Straight-Line Graph (PSLG) is the standard input representation
//! for constrained Delaunay triangulation.

pub mod graph;
pub mod segment;
pub mod vertex;

pub use graph::Pslg;
pub use segment::PslgSegmentId;
pub use vertex::PslgVertexId;
