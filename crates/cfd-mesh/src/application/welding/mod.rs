//! Vertex welding and snapping.
//!
//! Provides spatial hash-based vertex deduplication (also inside `VertexPool`)
//! and additional snap-to-grid / snap-to-vertex utilities for post-processing.
//!
//! ## Primary API
//!
//! - **[`SnappingGrid`]** — unified 26-neighbor grid for coordinate snapping
//!   AND vertex deduplication in one structure.  Integrates with
//!   `HalfEdgeMesh<'id>` via [`SnappingGrid::insert_or_weld_he`].
//!
//! ## Legacy API
//!
//! - [`SpatialHashGrid`] — point-set dedup only (no coordinate snapping).
//! - [`SnapConfig`] — coordinate snapping only (no dedup); shim in `snap.rs`.
//! - [`MeshWelder`] — batch vertex welding for existing `IndexedMesh`es.

pub mod snap;
pub mod spatial_hash;
pub mod welder;

pub use snap::{GridCell, SnapConfig, SnappingGrid};
pub use spatial_hash::SpatialHashGrid;
pub use welder::MeshWelder;
