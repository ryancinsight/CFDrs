//! Vertex welding and snapping.
//!
//! Provides spatial hash-based vertex deduplication (already in `VertexPool`)
//! and additional snap-to-grid / snap-to-vertex utilities for post-processing.

pub mod spatial_hash;
pub mod welder;
pub mod snap;

pub use welder::MeshWelder;
pub use snap::SnapConfig;
